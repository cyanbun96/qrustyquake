// The renderer's worker threads (r_threads). The screen is cut into bands of
// rows. The main thread records the drawing of a frame as commands, each for
// the bands its rows touch, and the workers run the commands of a band in the
// order they were recorded, one thread per band at a time. So every pixel goes
// through the same steps in the same order as when one thread draws it all.
// The commands only write the pixels of their band and draw state kept in
// Q_TLS variables, which the main thread only reads right after writing it.
// A worker that finds nothing to do for a moment sleeps until more is
// published, and r_threads -1 picks the count while playing (R_TuneFrame).
#ifdef __linux__
#define _GNU_SOURCE // sched_getaffinity
#include <sched.h>
#endif
#include "quakedef.h"
#include <fenv.h>

#define R_BANDSPERWORKER 4
#define R_STACKSIZE (4 << 20) // D_PolysetDraw keeps a span per screen row
#define R_SPINNS 50000 // how long a worker looks for work before it sleeps
#define R_AUTOMAX 8 // the most workers r_threads -1 tries, on 16 CPUs the
	// last few added 1-3% each

typedef struct rcmd_s {
	struct rcmd_s *next;
	rcmdfn_t fn;
	void *arg;
} rcmd_t;
typedef union { // a cache line pair of its own, workers poll every band
	struct { // the counts only ever go up, so a stale read is harmless
		SDL_AtomicInt pushed; // commands published by the main thread
		SDL_AtomicInt done; // commands that have run
		SDL_AtomicInt lease; // held by the thread running the commands
		rcmd_t *cur; // the last command that ran, the lease holder's
		rcmd_t *tail; // the last command added, the main thread's
		s32 added; // added but not published yet, the main thread's
		s32 y0, y1;
		rcmd_t head;
	} b;
	u8 pad[128];
} rband_t;
typedef struct {
	SDL_Thread *thread;
	SDL_Semaphore *wake;
	SDL_AtomicInt sleeping; // waiting for wake, see R_Sleep
} rworker_t;

Q_TLS s32 d_bandy0 = 0, d_bandy1 = MAXHEIGHT; // the rows being drawn
#ifdef __ANDROID__ // Bionic keeps the 64 bytes after the thread pointer for
	// itself on arm and won't start a program whose TLS is aligned less
static Q_TLS u8 r_tlsalign __attribute__((aligned(64), used));
#endif
bool r_threaded;
s32 r_nbands, r_bandrows, r_recordband;
static rband_t r_bands[R_MAXBANDS];
static rworker_t r_workers[R_MAXWORKERS];
static s32 r_nworkers, r_wanted; // the threads running and asked for
static s32 r_nactive; // the ones of those drawing this frame
static SDL_AtomicInt r_pollactive; // r_nactive for the workers
static SDL_AtomicInt r_pollbands; // r_nbands for the workers
static SDL_AtomicInt r_active; // a frame is being drawn
static SDL_AtomicInt r_quit;
static union { SDL_AtomicInt n; u8 pad[128]; } r_sleepers; // active ones
	// waiting to be woken, read at every R_BandPublish
static u64 r_framestart; // when the frame began, for r_threads -1
static u8 *r_arena; // the commands and their data, for one frame
static size_t r_arenasize, r_arenaused;
static void **r_oldarenas; // outgrown ones, freed when the frame is over
static s32 r_numold, r_capold;
static fenv_t r_fenv; // the rounding R_RenderView draws with, every thread
static s32 r_fenvframe; // counts the times r_fenv was taken

static bool R_BandPending(rband_t *band)
{
	return SDL_GetAtomicInt(&band->b.pushed)
		!= SDL_GetAtomicInt(&band->b.done);
}

// Runs what the band has, unless another thread is running it already
static bool R_RunBand(rband_t *band, bool worker)
{
	if(!R_BandPending(band)
		|| !SDL_CompareAndSwapAtomicInt(&band->b.lease, 0, 1))
		return 0;
	// The last holder may have run it in between, and then the main thread
	// can be setting the band up for the next frame already
	if(!R_BandPending(band)){
		SDL_SetAtomicInt(&band->b.lease, 0);
		return 0;
	}
	static Q_TLS s32 envframe = -1; // the workers keep what they have
	fenv_t mainenv;
	if(!worker)
		fegetenv(&mainenv);
	if(!worker || envframe != r_fenvframe){
		fesetenv(&r_fenv);
		envframe = worker ? r_fenvframe : -1;
	}
	u32 n = SDL_GetAtomicInt(&band->b.pushed);
	u32 done = SDL_GetAtomicInt(&band->b.done);
	d_bandy0 = band->b.y0;
	d_bandy1 = band->b.y1;
	rcmd_t *c = band->b.cur;
	for(; done != n; done++){
		c = c->next;
		c->fn(c->arg);
	}
	band->b.cur = c;
	if(!worker)
		fesetenv(&mainenv);
	// SDL_SetAtomicInt is only an acquire barrier with GCC on ARM and
	// RISC-V. What the commands read and drew, and cur, come before done
	// and the lease, or the next thread finds an old cur and skips some.
	SDL_MemoryBarrierRelease();
	SDL_SetAtomicInt(&band->b.done, done);
	SDL_SetAtomicInt(&band->b.lease, 0);
	return 1;
}

// Runs bands until none has anything left, w is the worker or NULL for the
// main thread, which only waits without help. A worker goes on looking for as
// long as the frame lasts, unless it found nothing for R_SPINNS or it isn't
// drawing anymore.
static void R_Work(rworker_t *w, bool help)
{
	u64 idle = 0;
	for(u32 spins = 0;; spins++){
		bool any = 0, ran = 0;
		s32 n = SDL_GetAtomicInt(&r_pollbands);
		for(s32 i = 0; i < n && !ran; i++){
			any |= R_BandPending(&r_bands[i]);
			if(help)
				ran = R_RunBand(&r_bands[i], w != NULL);
		}
		if(ran){
			idle = 0;
			continue;
		}
		if(!any && (!w || !SDL_GetAtomicInt(&r_active)
				|| SDL_GetAtomicInt(&r_quit)))
			break;
		if(w && !(spins & 31)){ // the clock costs, now and then
			u64 now = SDL_GetTicksNS();
			if(!idle)
				idle = now;
			else if(now - idle > R_SPINNS || w - r_workers
					>= SDL_GetAtomicInt(&r_pollactive))
				break;
		}
		SDL_CPUPauseInstruction();
	}
	d_bandy0 = 0;
	d_bandy1 = MAXHEIGHT;
}

// Is there a band with commands that no thread is running?
static bool R_Claimable()
{
	s32 n = SDL_GetAtomicInt(&r_pollbands);
	for(s32 i = 0; i < n; i++)
		if(R_BandPending(&r_bands[i])
				&& !SDL_GetAtomicInt(&r_bands[i].b.lease))
			return 1;
	return 0;
}

// Waits to be woken. A drawing worker counts itself in r_sleepers and looks
// for work once more after its flag is up, while R_BandPublish adds commands
// before it looks at the count, so one of them sees what the other did. The
// lease holder of a band looks at it again before it stops.
static void R_Sleep(rworker_t *w)
{
	SDL_SetAtomicInt(&w->sleeping, 1);
	bool active = w - r_workers < SDL_GetAtomicInt(&r_pollactive);
	if(active){
		SDL_AddAtomicInt(&r_sleepers.n, 1);
		if(R_Claimable()
			&& SDL_CompareAndSwapAtomicInt(&w->sleeping, 1, 0)){
			SDL_AddAtomicInt(&r_sleepers.n, -1);
			return;
		}
	}
	SDL_WaitSemaphore(w->wake);
	if(active)
		SDL_AddAtomicInt(&r_sleepers.n, -1);
}

static s32 SDLCALL R_Worker(void *data)
{
	rworker_t *w = data;
	for(;;){
		R_Sleep(w);
		if(SDL_GetAtomicInt(&r_quit))
			return 0;
		R_Work(w, 1);
	}
}

static void R_WakeWorkers(s32 from, s32 to)
{
	for(s32 i = from; i < to; i++)
		if(SDL_CompareAndSwapAtomicInt(&r_workers[i].sleeping, 1, 0))
			SDL_SignalSemaphore(r_workers[i].wake);
}

static void R_WakeOne()
{ // for commands just published, if no worker may be looking
	for(s32 i = 0; i < r_nactive; i++)
		if(SDL_CompareAndSwapAtomicInt(&r_workers[i].sleeping, 1, 0)){
			SDL_SignalSemaphore(r_workers[i].wake);
			return;
		}
}

static void R_StopWorkers()
{
	SDL_SetAtomicInt(&r_active, 0);
	SDL_SetAtomicInt(&r_quit, 1);
	for(s32 i = 0; i < r_nworkers; i++)
		SDL_SignalSemaphore(r_workers[i].wake);
	for(s32 i = 0; i < r_nworkers; i++){
		SDL_WaitThread(r_workers[i].thread, NULL);
		SDL_DestroySemaphore(r_workers[i].wake);
	}
	r_nworkers = r_nactive = 0;
	SDL_SetAtomicInt(&r_pollactive, 0);
	SDL_SetAtomicInt(&r_sleepers.n, 0);
	SDL_SetAtomicInt(&r_quit, 0);
}

static void R_StartWorkers(s32 n, bool quiet)
{
#ifdef Q_NOTHREADS
	if(n && !quiet)
		Con_Printf("r_threads: this build has no thread support\n");
#else
	for(; r_nworkers < n; r_nworkers++){
		rworker_t *w = &r_workers[r_nworkers];
		SDL_SetAtomicInt(&w->sleeping, 0);
		if(!(w->wake = SDL_CreateSemaphore(0)))
			break;
		SDL_PropertiesID p = SDL_CreateProperties();
		SDL_SetPointerProperty(p,
			SDL_PROP_THREAD_CREATE_ENTRY_FUNCTION_POINTER,
			(void *)R_Worker);
		SDL_SetStringProperty(p, SDL_PROP_THREAD_CREATE_NAME_STRING,
			"qq_render");
		SDL_SetPointerProperty(p,
			SDL_PROP_THREAD_CREATE_USERDATA_POINTER, w);
		SDL_SetNumberProperty(p,
			SDL_PROP_THREAD_CREATE_STACKSIZE_NUMBER, R_STACKSIZE);
		w->thread = SDL_CreateThreadWithProperties(p);
		SDL_DestroyProperties(p);
		if(!w->thread){
			SDL_DestroySemaphore(w->wake);
			break;
		}
	}
	if(r_nworkers < n && !quiet)
		Con_Printf("r_threads: %d of %d threads started: %s\n",
				r_nworkers, n, SDL_GetError());
#endif
}

static void R_GrowArena(size_t size)
{ // the old one stays until the frame is over, commands point into it
	if(r_arena && r_numold == r_capold){
		r_capold = r_capold ? r_capold * 2 : 16;
		r_oldarenas = realloc(r_oldarenas,
				r_capold * sizeof(*r_oldarenas));
		if(!r_oldarenas)
			Sys_Error("R_FrameAlloc: out of memory");
	}
	if(r_arena)
		r_oldarenas[r_numold++] = r_arena;
	r_arenasize = q_max(r_arenasize * 2, size + (1 << 20));
	if(!(r_arena = malloc(r_arenasize)))
		Sys_Error("R_FrameAlloc: out of memory");
	r_arenaused = 0;
}

void *R_FrameAlloc(size_t size)
{
	size = (size + 15) & ~(size_t)15;
	if(r_arenaused + size > r_arenasize)
		R_GrowArena(size);
	void *p = r_arena + r_arenaused;
	r_arenaused += size;
	return p;
}

void R_BandAdd(s32 band, rcmdfn_t fn, void *arg)
{
	rband_t *b = &r_bands[band];
	rcmd_t *c = R_FrameAlloc(sizeof(*c));
	c->next = NULL;
	c->fn = fn;
	c->arg = arg;
	b->b.tail->next = c;
	b->b.tail = c;
	b->b.added++;
}

void R_BandPublish(s32 band)
{
	rband_t *b = &r_bands[band];
	if(!b->b.added)
		return;
	SDL_AddAtomicInt(&b->b.pushed, b->b.added);
	b->b.added = 0;
	if(SDL_GetAtomicInt(&r_sleepers.n))
		R_WakeOne();
}

// Waits until everything recorded so far has been drawn, the commands that
// haven't been published yet too. Anywhere: the main thread only waits, as
// running commands would change the draw state it may be setting up.
void R_ThreadsSync()
{
	if(!r_nworkers)
		return;
	for(s32 b = 0; b < r_nbands; b++)
		R_BandPublish(b);
	R_Work(NULL, 0);
}

// R_ThreadsSync with the main thread drawing too, where its draw state is free
void R_ThreadsHelp()
{
	if(!r_nworkers)
		return;
	for(s32 b = 0; b < r_nbands; b++)
		R_BandPublish(b);
	R_Work(NULL, 1);
}

void R_RowsAdd(s32 y0, s32 y1, rcmdfn_t fn, void *arg)
{ // rows y0 to y1 - 1
	y0 = q_max(y0, 0);
	y1 = q_min(y1, (s32)vid.height);
	if(y0 >= y1)
		return;
	for(s32 b = y0 / r_bandrows; b <= (y1 - 1) / r_bandrows; b++){
		R_BandAdd(b, fn, arg);
		R_BandPublish(b);
	}
}

static void R_SetupBands()
{ // nothing is pending, so the workers leave the bands alone
	s32 h = vid.height, n = R_BANDSPERWORKER * r_nactive;
	r_bandrows = q_max((h + n - 1) / n, 2);
	r_nbands = q_min((h + r_bandrows - 1) / r_bandrows, R_MAXBANDS);
	r_bandrows = (h + r_nbands - 1) / r_nbands;
	for(s32 i = 0; i < r_nbands; i++){
		rband_t *b = &r_bands[i];
		b->b.head.next = NULL;
		b->b.cur = b->b.tail = &b->b.head;
		b->b.added = 0;
		b->b.y0 = i * r_bandrows;
		b->b.y1 = q_min((i + 1) * r_bandrows, h);
	}
	SDL_SetAtomicInt(&r_pollbands, r_nbands);
}

// The logical CPUs the game may run on
s32 R_ThreadsCPUs()
{
	static s32 n;
	if(n)
		return n;
	n = q_max(SDL_GetNumLogicalCPUCores(), 1);
#if defined(__linux__) && defined(CPU_COUNT) // taskset, containers
	cpu_set_t set;
	if(!sched_getaffinity(0, sizeof(set), &set) && CPU_COUNT(&set) > 0)
		n = q_min(n, CPU_COUNT(&set));
#endif
	return n;
}

// r_threads -1 starts with a worker for every CPU but the main thread's, at
// most 6, and tunes up to R_AUTOMAX. On every CPU set tried, cores or
// hyperthreads, that was the best count or close to it.
s32 R_ThreadsAutoMax(){ return CLAMP(0, R_ThreadsCPUs() - 1, R_AUTOMAX); }
static s32 R_AutoStart(){ return q_min(R_ThreadsAutoMax(), 6); }

// The tuning of r_threads -1. Every so often windows of frames with one
// worker more or less take turns with windows of the current count, and each
// is compared with the two around it, so a scene that changes slowly takes
// nothing from either. The other count stays if it's 2% faster and won in all
// but one window, a load on the CPU from something else rarely does that. A
// count that's clearly slower is dropped after two. A new map or screen size
// makes it look again soon.
#define TUNE_FRAMES 12 // a window has at least these, the first two after a
#define TUNE_NS 30000000 // change don't count, and lasts at least this long
#define TUNE_MAXFRAMES 64
#define TUNE_PAIRS 6 // trial windows, between TUNE_PAIRS + 1 current ones
#define TUNE_SOON 30 // frames until the next trial after a change
#define TUNE_MAXWAIT 3600 // and the longest, when neither neighbour helps
static struct {
	bool on;
	s32 cur, trial, dir; // trial < 0: none running
	s32 window, frame;
	u64 start; // of the window
	f32 med[2 * TUNE_PAIRS + 1]; // median frame time of each window
	f32 t[TUNE_MAXFRAMES];
	s32 wait, backoff, fails;
	u32 w, h;
	void *world;
} r_tune;

static f32 R_Median(f32 *v, s32 n)
{
	for(s32 i = 1; i < n; i++) // they're few
		for(s32 j = i; j > 0 && v[j - 1] > v[j]; j--){
			f32 t = v[j]; v[j] = v[j - 1]; v[j - 1] = t;
		}
	return n & 1 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}

static void R_TuneSoon(s32 wait)
{ // no trial running, the next one after wait frames
	r_tune.trial = -1;
	r_tune.wait = wait;
	r_tune.backoff = TUNE_SOON * 4;
	r_tune.fails = 0;
}

static void R_TuneReset()
{
	r_tune.on = r_threads.value < 0;
	r_tune.cur = q_min(R_AutoStart(), r_nworkers);
	r_tune.dir = 1;
	R_TuneSoon(TUNE_SOON * 2);
}

static s32 R_TuneCount()
{ return r_tune.trial >= 0 && r_tune.window & 1 ? r_tune.trial : r_tune.cur; }

// Trial window i over the current ones around it
static f32 R_TunePair(s32 i)
{ return r_tune.med[2*i + 1] * 2 / (r_tune.med[2*i] + r_tune.med[2*i + 2]); }

static void R_TuneFail()
{ // the other way next, and when neither helps, longer and longer waits
	r_tune.dir = -r_tune.dir;
	r_tune.wait = TUNE_SOON;
	if(++r_tune.fails >= 2){
		r_tune.wait = r_tune.backoff;
		r_tune.backoff = q_min(r_tune.backoff * 2, TUNE_MAXWAIT);
	}
}

static void R_TuneDecide(s32 pairs)
{
	f32 r[TUNE_PAIRS];
	s32 won = 0;
	for(s32 i = 0; i < pairs; i++){
		r[i] = R_TunePair(i);
		won += r[i] < 1;
	}
	f32 ratio = R_Median(r, pairs);
	s32 from = r_tune.cur, to = r_tune.trial;
	r_tune.trial = -1;
	if(pairs == TUNE_PAIRS && ratio < 0.98f && won >= TUNE_PAIRS - 1){
		r_tune.cur = to; // and on the same way, soon
		R_TuneSoon(TUNE_SOON);
	}else
		R_TuneFail();
	Con_DPrintf("r_threads -1: %d workers %.1f%% the time of %d, "
		"using %d\n", to, ratio * 100, from, r_tune.cur);
}

static void R_TuneFrame(u64 now, f32 ms)
{
	if(vid.width != r_tune.w || vid.height != r_tune.h
			|| (void *)cl.worldmodel != r_tune.world){
		r_tune.w = vid.width;
		r_tune.h = vid.height;
		r_tune.world = cl.worldmodel;
		R_TuneSoon(TUNE_SOON);
		return;
	}
	if(r_tune.trial < 0){
		if(--r_tune.wait > 0)
			return;
		s32 to = r_tune.cur + r_tune.dir;
		if(to < 0 || to > r_nworkers){ // there's nothing that way
			R_TuneFail();
			return;
		}
		r_tune.trial = to;
		r_tune.window = r_tune.frame = 0;
		return;
	}
	if(!r_tune.frame)
		r_tune.start = now;
	if(r_tune.frame >= 2)
		r_tune.t[r_tune.frame - 2] = ms;
	if((++r_tune.frame < TUNE_FRAMES || now - r_tune.start < TUNE_NS)
			&& r_tune.frame < TUNE_MAXFRAMES + 2)
		return;
	r_tune.med[r_tune.window] = R_Median(r_tune.t, r_tune.frame - 2);
	r_tune.frame = 0;
	if(++r_tune.window == 2 * TUNE_PAIRS + 1)
		R_TuneDecide(TUNE_PAIRS);
	else if(r_tune.window == 5 && R_TunePair(0) > 1.05f
			&& R_TunePair(1) > 1.05f)
		R_TuneDecide(2);
}

// The workers drawing this frame
s32 R_ThreadsActive(){ return r_nactive; }

// Decides whether this frame gets drawn by the workers, and how many
void R_ThreadsBeginFrame()
{
	R_ThreadsSync(); // in case a frame was cut short
	r_threaded = 0;
	for(; r_numold > 0; r_numold--)
		free(r_oldarenas[r_numold - 1]);
	r_arenaused = 0;
	bool tune = r_threads.value < 0;
	s32 want = tune ? R_ThreadsAutoMax()
		: CLAMP(0, (s32)r_threads.value, R_MAXWORKERS);
	if(want != r_wanted){
		R_StopWorkers();
		R_StartWorkers(want, tune);
		r_wanted = want;
		R_TuneReset();
	}
	if(tune != r_tune.on)
		R_TuneReset();
	s32 active = tune ? R_TuneCount() : r_nworkers, was = r_nactive;
	r_nactive = active;
	SDL_SetAtomicInt(&r_pollactive, active);
	R_WakeWorkers(active, was); // those that stop sleep again uncounted
	r_framestart = tune ? SDL_GetTicksNS() : 0;
	if(!active || r_showtris.value)
		return; // the debug lines are drawn in the middle of models
	R_SetupBands();
	r_threaded = 1;
	SDL_SetAtomicInt(&r_active, 1);
	R_WakeWorkers(0, active);
}

// The floating point state the frame is drawn with, before any command
void R_ThreadsFPEnv()
{
	if(!r_threaded)
		return;
	fegetenv(&r_fenv);
	r_fenvframe++;
}

// The frame's drawing is done. The workers stay for VID_Update's composite.
void R_ThreadsEndFrame()
{
	if(!r_threaded)
		return;
	R_ThreadsHelp();
	r_threaded = 0;
}

// Are the workers and the bands of the frame still there for more?
bool R_ThreadsAwake()
{ return r_nworkers && SDL_GetAtomicInt(&r_active); }

void R_ThreadsIdle()
{ // the workers go to sleep once they find nothing more to do
	R_ThreadsSync();
	SDL_SetAtomicInt(&r_active, 0);
	r_framestart = 0;
}

// VID_Update has made the frame
void R_ThreadsFrameDone()
{
	u64 start = r_framestart;
	R_ThreadsIdle();
	u64 now = SDL_GetTicksNS();
	if(start && r_tune.on && !r_showtris.value)
		R_TuneFrame(now, (now - start) / 1e6f);
}

void R_ThreadsShutdown()
{
	R_ThreadsSync();
	R_StopWorkers();
}
