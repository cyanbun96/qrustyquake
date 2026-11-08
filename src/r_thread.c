// The renderer's worker threads (r_threads). The screen is cut into bands of
// rows. The main thread records the drawing of a frame as commands, each for
// the bands its rows touch, and the workers run the commands of a band in the
// order they were recorded, one thread per band at a time. So every pixel goes
// through the same steps in the same order as when one thread draws it all.
// The commands only write the pixels of their band and draw state kept in
// Q_TLS variables, which the main thread only reads right after writing it.
#include "quakedef.h"
#include <fenv.h>

#define R_MAXWORKERS 64
#define R_BANDSPERWORKER 4
#define R_STACKSIZE (4 << 20) // D_PolysetDraw keeps a span per screen row

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
	SDL_AtomicInt sleeping;
} rworker_t;

Q_TLS s32 d_bandy0 = 0, d_bandy1 = MAXHEIGHT; // the rows being drawn
bool r_threaded;
s32 r_nbands, r_bandrows, r_recordband;
static rband_t r_bands[R_MAXBANDS];
static rworker_t r_workers[R_MAXWORKERS];
static s32 r_nworkers, r_wanted;
static SDL_AtomicInt r_pollbands; // r_nbands for the workers
static SDL_AtomicInt r_active; // a frame is being drawn, workers stay awake
static SDL_AtomicInt r_quit;
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
	SDL_SetAtomicInt(&band->b.done, done);
	SDL_SetAtomicInt(&band->b.lease, 0);
	return 1;
}

// Runs bands until none has anything left; workers wait for more for as long
// as the frame lasts. Without help the main thread only waits.
static void R_Work(bool worker, bool help)
{
	for(;;){
		bool any = 0, ran = 0;
		s32 n = SDL_GetAtomicInt(&r_pollbands);
		for(s32 i = 0; i < n && !ran; i++){
			any |= R_BandPending(&r_bands[i]);
			if(help)
				ran = R_RunBand(&r_bands[i], worker);
		}
		if(ran)
			continue;
		if(!any && (!worker || !SDL_GetAtomicInt(&r_active)
				|| SDL_GetAtomicInt(&r_quit)))
			break;
		SDL_CPUPauseInstruction();
	}
	d_bandy0 = 0;
	d_bandy1 = MAXHEIGHT;
}

static s32 SDLCALL R_Worker(void *data)
{
	rworker_t *w = data;
	for(;;){
		SDL_SetAtomicInt(&w->sleeping, 1);
		if(!SDL_GetAtomicInt(&r_active) && !SDL_GetAtomicInt(&r_quit))
			SDL_WaitSemaphore(w->wake);
		SDL_SetAtomicInt(&w->sleeping, 0);
		if(SDL_GetAtomicInt(&r_quit))
			return 0;
		R_Work(1, 1);
	}
}

static void R_WakeWorkers()
{
	for(s32 i = 0; i < r_nworkers; i++)
		if(SDL_CompareAndSwapAtomicInt(&r_workers[i].sleeping, 1, 0))
			SDL_SignalSemaphore(r_workers[i].wake);
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
	r_nworkers = 0;
	SDL_SetAtomicInt(&r_quit, 0);
}

static void R_StartWorkers(s32 n)
{
#ifdef Q_NOTHREADS
	if(n)
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
	if(r_nworkers < n)
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
	if(b->b.added){
		SDL_AddAtomicInt(&b->b.pushed, b->b.added);
		b->b.added = 0;
	}
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
	R_Work(0, 0);
}

// R_ThreadsSync with the main thread drawing too, where its draw state is free
void R_ThreadsHelp()
{
	if(!r_nworkers)
		return;
	for(s32 b = 0; b < r_nbands; b++)
		R_BandPublish(b);
	R_Work(0, 1);
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
	s32 h = vid.height, n = R_BANDSPERWORKER * r_nworkers;
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

// Decides whether this frame gets drawn by the workers
void R_ThreadsBeginFrame()
{
	R_ThreadsSync(); // in case a frame was cut short
	r_threaded = 0;
	for(; r_numold > 0; r_numold--)
		free(r_oldarenas[r_numold - 1]);
	r_arenaused = 0;
	s32 want = CLAMP(0, (s32)r_threads.value, R_MAXWORKERS);
	if(want != r_wanted){
		R_StopWorkers();
		R_StartWorkers(want);
		r_wanted = want;
	}
	if(!r_nworkers || r_showtris.value)
		return; // the debug lines are drawn in the middle of models
	R_SetupBands();
	r_threaded = 1;
	SDL_SetAtomicInt(&r_active, 1);
	R_WakeWorkers();
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
}

void R_ThreadsShutdown()
{
	R_ThreadsSync();
	R_StopWorkers();
}
