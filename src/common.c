// Copyright(C) 1996-2001 Id Software, Inc.
// Copyright(C) 2002-2009 John Fitzgibbons and others
// Copyright(C) 2010-2014 QuakeSpasm developers
// GPLv3 See LICENSE for details.
// common.c -- misc functions used in client and server
#include "quakedef.h"
#include "miniz.h"

/* All of Quake's data access is through a hierchal file system, but the
contents of the file system can be transparently merged from several sources.

The "base directory" is the path to the directory holding the quake.exe and all
game directories. The sys_ files pass this to host_init in quakeparms_t->basedir
This can be overridden with the "-basedir" command line parm to allow code
debugging in a different directory. The base directory is only used during
filesystem initialization.

The "game directory" is the first tree on the search path and directory that all
generated files(savegames, screenshots, demos, config files) will be saved to.
This can be overridden with the "-game" command line parameter. The game
directory can never be changed while quake is executing. This is a precacution
against having a malicious server instruct clients to write files over areas
they shouldn't.

The "cache directory" is only used during development to save network bandwidth,
especially over ISDN / T1 lines. If there is a cache directory specified, when
a file is found by the normal search path, it will be mirrored into the cache
directory, then opened there.

FIXME:
The file "parms.txt" will be read out of the game directory and appended to the
current command line arguments to allow different games to initialize startup
parms differently. This could be used to add a "-sspeed 22050" for the high
quality sound edition. Because they are added at the end, they will not
override an explicit setting on the original command line. */

static bool com_modified; // set 1 if using non-id files
static s8 *largv[MAX_NUM_ARGVS + 1];
static s8 argvdummy[] = " ";
static void COM_Path_f();
static u8 *loadbuf;
static cache_user_t *loadcache;
static s32 loadsize;
static localization_t localization;
static bool fitzmode;
static s8 com_cmdline[CMDLINE_LENGTH];
static searchpath_t *com_base_searchpaths;
static searchpath_t *com_searchpaths;
static u16 pop[] =
{ // this graphic needs to be in the pak file to use registered features
	     0,     0,     0,     0,     0,     0,     0,     0,
	     0,     0,0x6600,     0,     0,     0,0x6600,     0,
	     0,0x0066,     0,     0,     0,     0,0x0067,     0,
	     0,0x6665,     0,     0,     0,     0,0x0065,0x6600,
	0x0063,0x6561,     0,     0,     0,     0,0x0061,0x6563,
	0x0064,0x6561,     0,     0,     0,     0,0x0061,0x6564,
	0x0064,0x6564,     0,0x6469,0x6969,0x6400,0x0064,0x6564,
	0x0063,0x6568,0x6200,0x0064,0x6864,     0,0x6268,0x6563,
	     0,0x6567,0x6963,0x0064,0x6764,0x0063,0x6967,0x6500,
	     0,0x6266,0x6769,0x6a68,0x6768,0x6a69,0x6766,0x6200,
	     0,0x0062,0x6566,0x6666,0x6666,0x6666,0x6562,     0,
	     0,     0,0x0062,0x6364,0x6664,0x6362,     0,     0,
	     0,     0,     0,0x0062,0x6662,     0,     0,     0,
	     0,     0,     0,0x0061,0x6661,     0,     0,     0,
	     0,     0,     0,     0,0x6500,     0,     0,     0,
	     0,     0,     0,     0,0x6400,     0,     0,     0
};

static s32 q_islower(s32 c){return(c >= 'a' && c <= 'z');}
static s32 q_isupper(s32 c){return(c >= 'A' && c <= 'Z');}
static s32 q_tolower(s32 c){return((q_isupper(c)) ? (c | ('a' - 'A')) : c);}
static s32 q_toupper(s32 c){return((q_islower(c)) ? (c & ~('a' - 'A')) : c);}
static s32 q_isdigit(s32 c){return(c >= '0' && c <= '9');}
static s32 q_isblank(s32 c){return(c == ' ' || c == '\t');}
static s32 q_isspace(s32 c){
	switch(c){
		case ' ': case '\t':
		case '\n': case '\r':
		case '\f': case '\v': return 1;
	}
	return 0;
}

void COM_CloseFile(s32 h);
void ClearLink(link_t *l){ l->prev = l->next = l; } // used for new headnodes
void RemoveLink(link_t *l){ l->next->prev = l->prev; l->prev->next = l->next;}

void InsertLinkBefore(link_t *l, link_t *before)
{
	l->next = before;
	l->prev = before->prev;
	l->prev->next = l;
	l->next->prev = l;
}

void InsertLinkAfter(link_t *l, link_t *after)
{
	l->next = after->next;
	l->prev = after;
	l->prev->next = l;
	l->next->prev = l;
}

void Vec_Grow(void **pvec, size_t element_size, size_t count)
{
	vec_header_t header;
	if(*pvec) header = VEC_HEADER(*pvec);
	else header.size = header.capacity = 0;
	if(header.size + count > header.capacity) {
		void *new_buffer;
		size_t total_size;
		header.capacity = header.size + count;
		header.capacity += header.capacity >> 1;
		if(header.capacity < 16) header.capacity = 16;
		total_size = sizeof(vec_header_t)+header.capacity*element_size;
		if(*pvec)
			new_buffer=realloc(((vec_header_t*)*pvec)-1,total_size);
		else new_buffer = malloc(total_size);
		if(!new_buffer)
			Sys_Error("Vec_Grow: failed to allocate %lu bytes\n",
					(u64)total_size);
		*pvec = 1 + (vec_header_t*)new_buffer;
		VEC_HEADER(*pvec) = header;
	}
}

void Vec_Append(void **pvec, size_t element_size, const void *data, size_t count)
{
	if(!count) return;
	Vec_Grow(pvec, element_size, count);
	memcpy((u8 *)*pvec + VEC_HEADER(*pvec).size * element_size, data, count * element_size);
	VEC_HEADER(*pvec).size += count;
}

void Vec_Clear(void **pvec) { if(*pvec) VEC_HEADER(*pvec).size = 0; }
void Vec_Free(void **pvec){if(*pvec) { free(&VEC_HEADER(*pvec)); *pvec = NULL;}}


s32 q_strlcpy(s8 *dst, const s8 *src, size_t siz) // $OpenBSD: q_strlcpy.c,v1.11
{ // Copyright(c) 1998 Todd C. Miller <Todd.Miller@courtesan.com>
	s8 *d = dst;
	const s8 *s = src;
	size_t n = siz;
	if(n != 0){ // Copy as many bytes as will fit
		while(--n != 0){ if((*d++ = *s++) == '\0') break; }
	}
	if(n == 0){// Not enough room in dst, add NUL and traverse rest of src
		if(siz != 0) *d = '\0'; // NUL-terminate dst
		while(*s++);
	}
	return(s - src - 1); // count does not include NUL
}


size_t q_strlcat(s8 *dst, const s8 *src, size_t siz)//$OpenBSD:q_strlcat.c,v1.13
{
	s8 *d = dst;
	const s8 *s = src;
	size_t n = siz;
	// Find the end of dst and adjust bytes left but don't go past end
	while(n-- != 0 && *d != '\0') d++;
	size_t dlen = d - dst;
	n = siz - dlen;
	if(n == 0) return(dlen + strlen(s));
	while(*s != '\0'){ if(n != 1){ *d++ = *s; n--; } s++; }
	*d = '\0';
	return(dlen + (s - src)); // count does not include NUL
}

s32 q_strcasecmp(const s8 * s1, const s8 * s2)
{
	const s8 * p1 = s1;
	const s8 * p2 = s2;
	s8 c1, c2;
	if(p1 == p2) return 0;
	do {
		c1 = q_tolower(*p1++);
		c2 = q_tolower(*p2++);
		if(c1 == '\0') break;
	} while(c1 == c2);
	return(s32)(c1 - c2);
}

s32 q_strncasecmp(const s8 *s1, const s8 *s2, size_t n)
{
	const s8 * p1 = s1;
	const s8 * p2 = s2;
	s8 c1, c2;
	if(p1 == p2 || n == 0) return 0;
	do {
		c1 = q_tolower(*p1++);
		c2 = q_tolower(*p2++);
		if(c1 == '\0' || c1 != c2)
			break;
	} while(--n > 0);
	return(s32)(c1 - c2);
}

s8 *q_strcasestr(const s8 *haystack, const s8 *needle)
{
	const size_t len = strlen(needle);
	while(*haystack) {
		if(!q_strncasecmp(haystack, needle, len))
			return(s8 *)haystack;
		++haystack;
	}
	return NULL;
}

s8 *q_strlwr(s8 *str)
{
	s8 *c = str;
	while(*c) { *c = q_tolower(*c); c++; }
	return str;
}

s8 *q_strupr(s8 *str)
{
	s8 *c = str;
	while(*c) { *c = q_toupper(*c); c++; }
	return str;
}

s8 *q_strdup(const s8 *str)
{
	size_t len = strlen(str) + 1;
	s8 *newstr = (s8 *)malloc(len);
	memcpy(newstr, str, len);
	return newstr;
}


s32 q_vsnprintf(s8 *str, size_t size, const s8 *format, va_list args)
{
	s32 ret = vsnprintf_func(str, size, format, args);
	if(ret < 0) ret = (s32)size;
	if(size == 0) return ret; // no buffer
	if((size_t)ret >= size) str[size - 1] = '\0';
	return ret;
}

s32 q_snprintf(s8 *str, size_t size, const s8 *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	s32 ret = q_vsnprintf(str, size, format, argptr);
	va_end(argptr);
	return ret;
}

void Q_memset(void *dest, s32 fill, size_t count)
{
	size_t i;
	if( (((uintptr_t)dest | count) & 3) == 0) {
		count >>= 2;
		fill = fill | (fill<<8) | (fill<<16) | (fill<<24);
		for(i = 0; i < count; i++)
			((s32 *)dest)[i] = fill;
	}
	else for(i = 0; i < count; i++)
		((u8 *)dest)[i] = fill;
}

void Q_memmove(void *dest, const void *src, size_t count)
{ memmove(dest, src, count); }

void Q_memcpy(void *dest, const void *src, size_t count)
{ memcpy(dest, src, count); }

s32 Q_memcmp(const void *m1, const void *m2, size_t count)
{
	while(count) {
		count--;
		if(((u8 *)m1)[count] != ((u8 *)m2)[count])
			return -1;
	}
	return 0;
}

void Q_strcpy(s8 *dest, const s8 *src)
{
	while(*src) *dest++ = *src++;
	*dest++ = 0;
}

void Q_strncpy(s8 *dest, const s8 *src, s32 count)
{
	while(*src && count--)
		*dest++ = *src++;
	if(count) *dest++ = 0;
}

s32 Q_strlen(const s8 *str)
{
	s32 count = 0;
	while(str[count]) count++;
	return count;
}

s8 *Q_strrchr(const s8 *s, s8 c)
{
	s32 len = Q_strlen(s);
	s += len;
	while(len--) if(*--s == c) return(s8 *)s;
	return NULL;
}

void Q_strcat(s8 *dest, const s8 *src)
{
	dest += Q_strlen(dest);
	Q_strcpy(dest, src);
}

s32 Q_strcmp(const s8 *s1, const s8 *s2)
{
	while(1) {
		if(*s1 != *s2) return -1; // strings not equal
		if(!*s1) return 0; // strings are equal
		s1++;
		s2++;
	}
	return -1;
}

s32 Q_strncmp(const s8 *s1, const s8 *s2, s32 count)
{
	while(1) {
		if(!count--) return 0;
		if(*s1 != *s2) return -1; // strings not equal
		if(!*s1) return 0; // strings are equal
		s1++;
		s2++;
	}
	return -1;
}

s32 Q_atoi(const s8 *str)
{
	while(q_isspace(*str)) ++str;
	s32 sign = 1;
	if(*str == '-') {
		sign = -1;
		str++;
	}
	s32 val = 0;
	if(str[0] == '0' && (str[1] == 'x' || str[1] == 'X') ) { //check for hex
		str += 2;
		while(1) {
			s32 c = *str++;
			if(c >= '0' && c <= '9')
				val = (val<<4) + c - '0';
			else if(c >= 'a' && c <= 'f')
				val = (val<<4) + c - 'a' + 10;
			else if(c >= 'A' && c <= 'F')
				val = (val<<4) + c - 'A' + 10;
			else
				return val*sign;
		}
	}
	if(str[0] == '\'') // check for character
		return sign * str[1];
	while(1) { // assume decimal
		s32 c = *str++;
		if(c <'0' || c > '9')
			return val*sign;
		val = val*10 + c - '0';
	}
	return 0;
}


f32 Q_atof(const s8 *str)
{
	while(q_isspace(*str)) ++str;
	s32 sign = 1;
	if(*str == '-') {
		sign = -1;
		str++;
	}
	f64 val = 0;
	if(str[0] == '0' && (str[1] == 'x' || str[1] == 'X') ) { //check for hex
		str += 2;
		while(1) {
			s32 c = *str++;
			if(c >= '0' && c <= '9')
				val = (val*16) + c - '0';
			else if(c >= 'a' && c <= 'f')
				val = (val*16) + c - 'a' + 10;
			else if(c >= 'A' && c <= 'F')
				val = (val*16) + c - 'A' + 10;
			else
				return val*sign;
		}
	}
	if(str[0] == '\'') // check for character
		return sign * str[1];
	s32 decimal = -1; // assume decimal
	s32 total = 0;
	while(1) {
		s32 c = *str++;
		if(c == '.') {
			decimal = total;
			continue;
		}
		if(c <'0' || c > '9')
			break;
		val = val*10 + c - '0';
		total++;
	}
	if(decimal == -1)
		return val*sign;
	while(total > decimal) {
		val /= 10;
		total--;
	}
	return val*sign;
}

s16 ShortSwap(s16 l)
{
	u8 b1 = l&255;
	u8 b2 = (l>>8)&255;
	return(b1<<8) + b2;
}

s16 ShortNoSwap(s16 l) { return l; }

s32 LongSwap(s32 l)
{
	u8 b1 = l&255;
	u8 b2 = (l>>8)&255;
	u8 b3 = (l>>16)&255;
	u8 b4 = (l>>24)&255;
	return((s32)b1<<24) + ((s32)b2<<16) + ((s32)b3<<8) + b4;
}

s32 LongNoSwap(s32 l) { return l; }

f32 FloatSwap(f32 f)
{
	union { f32 f; u8 b[4]; } dat1, dat2;
	dat1.f = f;
	dat2.b[0] = dat1.b[3];
	dat2.b[1] = dat1.b[2];
	dat2.b[2] = dat1.b[1];
	dat2.b[3] = dat1.b[0];
	return dat2.f;
}

f32 FloatNoSwap(f32 f) { return f; }

void MSG_WriteChar(sizebuf_t *sb, s32 c)
{ // Handles u8 ordering and avoids alignment errors
	u8 *buf = (u8 *) SZ_GetSpace(sb, 1);
	buf[0] = c;
}

void MSG_WriteByte(sizebuf_t *sb, s32 c)
{
	u8 *buf = (u8 *) SZ_GetSpace(sb, 1);
	buf[0] = c;
}

void MSG_WriteShort(sizebuf_t *sb, s32 c)
{
	u8 *buf = (u8 *) SZ_GetSpace(sb, 2);
	buf[0] = c&0xff;
	buf[1] = c>>8;
}

void MSG_WriteLong(sizebuf_t *sb, s32 c)
{
	u8 *buf = (u8 *) SZ_GetSpace(sb, 4);
	buf[0] = c&0xff;
	buf[1] = (c>>8)&0xff;
	buf[2] = (c>>16)&0xff;
	buf[3] = c>>24;
}

void MSG_WriteFloat(sizebuf_t *sb, f32 f)
{
	union { f32 f; s32 l; } dat;
	dat.f = f;
	dat.l = LittleLong(dat.l);
	SZ_Write(sb, &dat.l, 4);
}

void MSG_WriteString(sizebuf_t *sb, const s8 *s)
{
	if(!s) SZ_Write(sb, "", 1);
	else SZ_Write(sb, s, Q_strlen(s)+1);
}

void MSG_WriteCoord16(sizebuf_t *sb, f32 f)
{ //johnfitz -- original behavior, 13.3 fixed point coords, max range +-4096
	MSG_WriteShort(sb, Q_rint(f*8));
}

void MSG_WriteCoord24(sizebuf_t *sb, f32 f)
{ //johnfitz -- 16.8 fixed point coords, max range +-32768
	MSG_WriteShort(sb, f);
	MSG_WriteByte(sb, (s32)(f*255)%255);
}

void MSG_WriteCoord32f(sizebuf_t *sb, f32 f)
{ //johnfitz -- 32-bit f32 coords
	MSG_WriteFloat(sb, f);
}

void MSG_WriteCoord(sizebuf_t *sb, f32 f, u32 flags)
{
	if(flags & PRFL_FLOATCOORD)      MSG_WriteFloat(sb, f);
	else if(flags & PRFL_INT32COORD) MSG_WriteLong(sb, Q_rint(f * 16));
	else if(flags & PRFL_24BITCOORD) MSG_WriteCoord24(sb, f);
	else MSG_WriteCoord16(sb, f);
}

void MSG_WriteAngle(sizebuf_t *sb, f32 f, u32 flags)
{
	if(flags & PRFL_FLOATANGLE)
		MSG_WriteFloat(sb, f);
	else if(flags & PRFL_SHORTANGLE)
		MSG_WriteShort(sb, Q_rint(f * 65536.0 / 360.0) & 65535);
	else MSG_WriteByte(sb, Q_rint(f * 256.0 / 360.0) & 255);
}
void MSG_WriteAngle16(sizebuf_t *sb, f32 f, u32 flags)
{ //johnfitz -- for PROTOCOL_FITZQUAKE
	if(flags & PRFL_FLOATANGLE)
		MSG_WriteFloat(sb, f);
	else MSG_WriteShort(sb, Q_rint(f * 65536.0 / 360.0) & 65535);
}

void MSG_BeginReading() { msg_readcount = 0; msg_badread = 0; }


s32 MSG_ReadChar()
{ // returns -1 and sets msg_badread if no more characters are available
	if(msg_readcount+1 > net_message.cursize) {
		msg_badread = 1;
		return -1;
	}
	s32 c = (s8)net_message.data[msg_readcount];
	msg_readcount++;
	return c;
}

s32 MSG_ReadByte()
{
	if(msg_readcount+1 > net_message.cursize) {
		msg_badread = 1;
		return -1;
	}
	s32 c = (u8)net_message.data[msg_readcount];
	msg_readcount++;
	return c;
}

s32 MSG_ReadShort()
{
	if(msg_readcount+2 > net_message.cursize) {
		msg_badread = 1;
		return -1;
	}
	s32 c = (s16)(net_message.data[msg_readcount]
			+ (net_message.data[msg_readcount+1]<<8));
	msg_readcount += 2;
	return c;
}

s32 MSG_ReadLong()
{
	if(msg_readcount+4 > net_message.cursize) {
		msg_badread = 1;
		return -1;
	}
	s32 c = net_message.data[msg_readcount]
		+ (net_message.data[msg_readcount+1]<<8)
		+ (net_message.data[msg_readcount+2]<<16)
		+ (net_message.data[msg_readcount+3]<<24);
	msg_readcount += 4;
	return c;
}

f32 MSG_ReadFloat()
{
	union { u8 b[4]; f32 f; s32 l; } dat;
	dat.b[0] = net_message.data[msg_readcount];
	dat.b[1] = net_message.data[msg_readcount+1];
	dat.b[2] = net_message.data[msg_readcount+2];
	dat.b[3] = net_message.data[msg_readcount+3];
	msg_readcount += 4;
	dat.l = LittleLong(dat.l);
	return dat.f;
}

const s8 *MSG_ReadString()
{
	static s8 string[2048];
	size_t l = 0;
	do {
		s32 c = MSG_ReadByte();
		if(c == -1 || c == 0) break;
		string[l] = c;
		l++;
	} while(l < sizeof(string) - 1);
	string[l] = 0;
	return string;
}

f32 MSG_ReadCoord16()
{ //johnfitz -- original behavior, 13.3 fixed point coords, max range +-4096
	return MSG_ReadShort() * (1.0/8);
}

f32 MSG_ReadCoord24()
{ //johnfitz -- 16.8 fixed point coords, max range +-32768
	return MSG_ReadShort() + MSG_ReadByte() * (1.0/255);
}


f32 MSG_ReadCoord32f()
{ //johnfitz -- 32-bit f32 coords
	return MSG_ReadFloat();
}

f32 MSG_ReadCoord(u32 flags)
{
	if(flags & PRFL_FLOATCOORD)      return MSG_ReadFloat();
	else if(flags & PRFL_INT32COORD) return MSG_ReadLong() * (1.0 / 16.0);
	else if(flags & PRFL_24BITCOORD) return MSG_ReadCoord24();
	else return MSG_ReadCoord16();
}

f32 MSG_ReadAngle(u32 flags)
{
	if(flags & PRFL_FLOATANGLE)      return MSG_ReadFloat();
	else if(flags & PRFL_SHORTANGLE) return MSG_ReadShort() * (360.0/65536);
	else return MSG_ReadChar() * (360.0 / 256);
}


f32 MSG_ReadAngle16(u32 flags)
{ //johnfitz -- for PROTOCOL_FITZQUAKE
	if(flags & PRFL_FLOATANGLE)
		return MSG_ReadFloat(); // make sure
	else return MSG_ReadShort() * (360.0 / 65536);
}

void SZ_Alloc(sizebuf_t *buf, s32 startsize)
{
	if(startsize < 256) startsize = 256;
	buf->data = (u8 *) Hunk_AllocName(startsize, "sizebuf");
	buf->maxsize = startsize;
	buf->cursize = 0;
}

void SZ_Free(sizebuf_t *buf) { buf->cursize = 0; }
void SZ_Clear(sizebuf_t *buf) { buf->cursize = 0; }

void *SZ_GetSpace(sizebuf_t *buf, s32 length)
{
	if(buf->cursize + length > buf->maxsize) {
		if(!buf->allowoverflow)
		Host_Error("SZ_GetSpace: overflow without allowoverflow set");
		if(length > buf->maxsize)
		Sys_Error("SZ_GetSpace: %i is > full buffer size", length);
		buf->overflowed = 1;
		Con_Printf("SZ_GetSpace: overflow\n");
		SZ_Clear(buf);
	}
	void *data = buf->data + buf->cursize;
	buf->cursize += length;
	return data;
}

void SZ_Write(sizebuf_t *buf, const void *data, s32 length)
{ Q_memcpy(SZ_GetSpace(buf,length),data,length); }

void SZ_Print(sizebuf_t *buf, const s8 *data)
{
	s32 len = Q_strlen(data) + 1;
	if(buf->data[buf->cursize-1]) // no trailing 0
		Q_memcpy((u8 *)SZ_GetSpace(buf, len ) , data, len);
	else // write over trailing 0
		Q_memcpy((u8 *)SZ_GetSpace(buf, len-1)-1, data, len);
}

const s8 *COM_SkipPath(const s8 *pathname)
{
	const s8 *last = pathname;
	while(*pathname) {
		if(*pathname == '/')
			last = pathname + 1;
		pathname++;
	}
	return last;
}

void COM_StripExtension(const s8 *in, s8 *out, size_t outsize)
{
	if(!*in) {
		*out = '\0';
		return;
	}
	if(in != out) // copy when not in-place editing
		q_strlcpy(out, in, outsize);
	s32 length = (s32)strlen(out) - 1;
	while(length > 0 && out[length] != '.') {
		--length;
		if(out[length] == '/' || out[length] == '\\')
			return; // no extension
	}
	if(length > 0) out[length] = '\0';
}

const s8 *COM_FileGetExtension(const s8 *in)
{ // COM_FileGetExtension - doesn't return NULL
	const s8 *src;
	size_t len = strlen(in);
	if(len < 2) // nothing meaningful
		return "";
	src = in + len - 1;
	while(src != in && src[-1] != '.')
		src--;
	if(src == in || strchr(src, '/') != NULL || strchr(src, '\\') != NULL)
		return ""; // no extension, or parent directory has a dot 
	return src;
}

void COM_ExtractExtension(const s8 *in, s8 *out, size_t outsize)
{
	const s8 *ext = COM_FileGetExtension(in);
	if(! *ext) *out = '\0';
	else q_strlcpy(out, ext, outsize);
}

void COM_FileBase(const s8 *in, s8 *out, size_t outsize)
{ // take 'somedir/otherdir/filename.ext', write only 'filename' to the output
	const s8 *s = in;
	const s8 *slash = in;
	const s8 *dot = NULL;
	while(*s) {
		if(*s == '/' || *s == '\\')
			slash = s + 1;
		if(*s == '.')
			dot = s;
		s++;
	}
	if(dot == NULL) dot = s;
	if(dot - slash < 2) q_strlcpy(out, "?model?", outsize);
	else {
		size_t len = dot - slash;
		if(len >= outsize)
			len = outsize - 1;
		memcpy(out, slash, len);
		out[len] = '\0';
	}
}

void COM_AddExtension(s8 *path, const s8 *extension, size_t len)
{ // if path extension != .EXT, append it (extension should include leading '.')
	if(strcmp(COM_FileGetExtension(path), extension + 1) != 0)
		q_strlcat(path, extension, len);
}

// Parse a token out of a string
// The mode argument controls how overflow is handled:
// - CPE_NOTRUNC: return NULL(abort parsing)
// - CPE_ALLOWTRUNC: truncate com_token(ignore extra characters in this token)
const s8 *COM_ParseEx(const s8 *data, cpe_mode mode)
{
	s32 len = 0;
	com_token[0] = 0;
	if(!data) return NULL;
skipwhitespace:
	s32 c;
	while((c = *data) <= ' ') {
		if(c == 0) return NULL; // end of file
		data++;
	}
	if(c == '/' && data[1] == '/') { // skip // comments
		while(*data && *data != '\n') data++;
		goto skipwhitespace;
	}
	if(c == '/' && data[1] == '*') { // skip /*..*/ comments
		data += 2;
		while(*data && !(*data == '*' && data[1] == '/')) data++;
		if(*data) data += 2;
		goto skipwhitespace;
	}
	if(c == '\"') { // handle quoted strings specially
		data++;
		while(1) {
			if((c = *data) != 0) ++data;
			if(c == '\"' || !c) { com_token[len] = 0; return data; }
			if((u64)len<Q_COUNTOF(com_token)-1)com_token[len++] = c;
			else if(mode == CPE_NOTRUNC) return NULL;
		}
	} // parse single characters
	if(c == '{' || c == '}'|| c == '('|| c == ')' || c == '\'' || c == ':'){
		if((u64)len < Q_COUNTOF(com_token) - 1) com_token[len++] = c;
		else if(mode == CPE_NOTRUNC) return NULL;
		com_token[len] = 0;
		return data+1;
	}
	do { // parse a regular word
		if((u64)len < Q_COUNTOF(com_token) - 1) com_token[len++] = c;
		else if(mode == CPE_NOTRUNC) return NULL;
		data++;
		c = *data;
		if(c == '{'||c == '}'||c== '('||c == ')'||c == '\'') break;
	} while(c > 32);
	com_token[len] = 0;
	return data;
}

const s8 *COM_Parse(const s8 *data)
{ // Parse a token out of a string, return NULL in case of overflow
	return COM_ParseEx(data, CPE_NOTRUNC);
}


s32 COM_CheckParm(const s8 *parm) // Returns the position(1 to argc-1) in the
{ // program's argument list where given parameter apears, or 0 if not present
	for(s32 i = 1; i < com_argc; i++) {
		if(!com_argv[i]) continue;
		if(!Q_strcmp(parm,com_argv[i])) return i;
	}
	return 0;
}

// Looks for the pop.txt file and verifies it.
// Sets the "registered" cvar.
// Immediately exits out if an alternate game was attempted to be started
// without being registered.
static void COM_CheckRegistered()
{
	s32 h;
	u16 check[128];
	s32 i;
	COM_OpenFile("gfx/pop.lmp", &h, NULL);
	if(h == -1) {
		Cvar_Set("registered", "0");
		Con_Printf("Playing shareware version.\n");
		if(com_modified)
Sys_Error("You must have the registered version to use modified games.\n\n"
"Basedir is: %s\n\n"
"Check that this has an \"id1\" subdirectory containing pak0.pak and pak1.pak, "
"or use the -basedir command-line option to specify another directory.",
		com_basedir);
		return;
	}
	i = Sys_FileRead(h, check, sizeof(check));
	COM_CloseFile(h);
	if(i != (s32) sizeof(check)) goto corrupt;
	for(i = 0; i < 128; i++) {
		if(pop[i] != (u16)BigShort(check[i]))
		{ corrupt:
			Sys_Error("Corrupted data file.");
		}
	}
	for(i = 0; com_cmdline[i]; i++)
		if(com_cmdline[i]!= ' ')
			break;
	Cvar_Set("cmdline", &com_cmdline[i]);
	Cvar_Set("registered", "1");
	Con_Printf("Playing registered version.\n");
}


void COM_InitArgv(s32 argc, s8 **argv)
{
	// reconstitute the command line for the cmdline externally visible cvar
	s32 n = 0;
	for(s32 j = 0; (j<MAX_NUM_ARGVS) && (j< argc); j++)
	{
		s32 i = 0;
		while((n < (CMDLINE_LENGTH - 1)) && argv[j][i])
			com_cmdline[n++] = argv[j][i++];
		if(n < (CMDLINE_LENGTH - 1))
			com_cmdline[n++] = ' ';
		else break;
	}
	if(n > 0 && com_cmdline[n-1] == ' ') com_cmdline[n-1] = 0;
	Con_Printf("Command line: %s\n", com_cmdline);
	for(com_argc=0; (com_argc<MAX_NUM_ARGVS)&&(com_argc<argc); com_argc++) {
		largv[com_argc] = argv[com_argc];
		if(!Q_strcmp("-safe", argv[com_argc]))
			safemode = 1;
	}
	largv[com_argc] = argvdummy;
	com_argv = largv;
	if(COM_CheckParm("-rogue")) {
		rogue = 1;
		standard_quake = 0;
	}
	if(COM_CheckParm("-hipnotic") || COM_CheckParm("-quoth")) {
		hipnotic = 1;
		standard_quake = 0;
	}
}

void COM_Init()
{
	BigShort = ShortSwap;
	LittleShort = ShortNoSwap;
	BigLong = LongSwap;
	LittleLong = LongNoSwap;
	BigFloat = FloatSwap;
	LittleFloat = FloatNoSwap;
	fitzmode = COM_CheckParm("-fitz");
}

// Does a varargs printf into a temp buffer. Cycles between 4 different static
// buffers. The number of buffers cycled is defined in VA_NUM_BUFFS.
static s8 *get_va_buffer()
{
	static s8 va_buffers[VA_NUM_BUFFS][VA_BUFFERLEN];
	static s32 buffer_idx = 0;
	buffer_idx = (buffer_idx + 1) & (VA_NUM_BUFFS - 1);
	return va_buffers[buffer_idx];
}

s8 *va(const s8 *format, ...)
{
	va_list argptr;
	s8 *va_buf = get_va_buffer();
	va_start(argptr, format);
	q_vsnprintf(va_buf, VA_BUFFERLEN, format, argptr);
	va_end(argptr);
	return va_buf;
}

static void COM_Path_f()
{
	Con_Printf("Current search path:\n");
	for(searchpath_t *s = com_searchpaths; s; s = s->next) {
		if(s->pack) Con_Printf("%s(%i files)\n",
				s->pack->filename, s->pack->numfiles);
		else Con_Printf("%s\n", s->filename);
	}
}


void COM_WriteFile(const s8 *filename, const void *data, s32 len)
{ // The filename will be prefixed by the current game directory
	s32 handle;
	s8 name[MAX_OSPATH];
	Sys_mkdir(com_gamedir); // create if nonexistant gamedir to avoid crash
	q_snprintf(name, sizeof(name), "%s/%s", com_gamedir, filename);
	handle = Sys_FileOpenWrite(name);
	if(handle == -1) {
		Sys_Printf("COM_WriteFile: failed on %s\n", name);
		return;
	}
	Sys_Printf("COM_WriteFile: %s\n", name);
	Sys_FileWrite(handle, data, len);
	Sys_FileClose(handle);
}

void COM_CreatePath(s8 *path)
{
	for(s8 *ofs = path + 1; *ofs; ofs++) {
		if(*ofs == '/') { // create the directory
			*ofs = 0;
			Sys_mkdir(path);
			*ofs = '/';
		}
	}
}

s64 COM_filelength(FILE *f)
{
	s64 pos = ftell(f);
	fseek(f, 0, SEEK_END);
	s64 end = ftell(f);
	fseek(f, pos, SEEK_SET);
	return end;
}

// Finds the file in the search path. Sets com_filesize and one of handle or file If neither of file or handle is set, this can be used for detecting a file's presence.
static s32 COM_FindFile(const s8 *name, s32 *handle, FILE **file, u32 *path_id)
{
	s8 netpath[MAX_OSPATH];
	s32 i;
	if(file && handle) Sys_Error("COM_FindFile: both handle and file set");
	file_from_pak = 0;
	// search through the path, one element at a time
	for(searchpath_t *search=com_searchpaths; search; search=search->next) {
		if(search->pack){ // look through all the pak file elements
			pack_t *pak = search->pack;
			for(i = 0; i < pak->numfiles; i++) {
				if(strcmp(pak->files[i].name, name) != 0)
					continue;
				com_filesize = pak->files[i].filelen;//found it!
				file_from_pak = 1;
				if(path_id) *path_id = search->path_id;
				if(handle){
					*handle = pak->handle;
					Sys_FileSeek(pak->handle,
							pak->files[i].filepos);
					return com_filesize;
				}
				else if(file){ // open a new file on the pakfile
					*file = fopen(pak->filename, "rb");
					if(*file)
						fseek(*file,
						pak->files[i].filepos,SEEK_SET);
					return com_filesize;
				}
				else return com_filesize;
			}
		}
		else { // check a file in the directory tree
			if(!registered.value)
				if(strchr(name, '/') || strchr(name,'\\'))
					continue;
			q_snprintf(netpath, sizeof(netpath),
					"%s/%s",search->filename, name);
			if(!(Sys_FileType(netpath) & FS_ENT_FILE)) continue;
			if(path_id) *path_id = search->path_id;
			if(handle) {
				com_filesize = Sys_FileOpenRead(netpath, &i);
				*handle = i;
				return com_filesize;
			}
			else if(file) {
				*file = fopen(netpath, "rb");
				com_filesize = (*file == NULL) ?
					-1 : COM_filelength(*file);
				return com_filesize;
			}
			else return 0;
		}
	}
	if(strcmp(COM_FileGetExtension(name), "pcx") != 0
		&& strcmp(COM_FileGetExtension(name), "tga") != 0
		&& strcmp(COM_FileGetExtension(name), "lit") != 0
		&& strcmp(COM_FileGetExtension(name), "vis") != 0
		&& strcmp(COM_FileGetExtension(name), "ent") != 0)
		Con_DPrintf("FindFile: can't find %s\n", name);
	else Con_DPrintf("FindFile: can't find %s\n", name);
	if(handle) *handle = -1;
	if(file) *file = NULL;
	com_filesize = -1;
	return com_filesize;
}

bool COM_FileExists(const s8 *filename, u32 *path_id)
{ // Returns whether the file is found in the quake filesystem.
	s32 ret = COM_FindFile(filename, NULL, NULL, path_id);
	return(ret == -1) ? 0 : 1;
}

// filename never has a leading slash, but may contain directory walks.
// Returns a handle and a length. It may actually be inside a pak file.
s32 COM_OpenFile(const s8 *filename, s32 *handle, u32 *path_id)
{ return COM_FindFile(filename, handle, NULL, path_id); }

// If requested file is inside a pak, a new FILE * will be opened into the file.
s32 COM_FOpenFile(const s8 *filename, FILE **file, u32 *path_id)
{ return COM_FindFile(filename, NULL, file, path_id); }

void COM_CloseFile(s32 h)
{ // If it is a pak file handle, don't really close it
	for(searchpath_t *s = com_searchpaths; s; s = s->next)
		if(s->pack && s->pack->handle == h)
			return;
	Sys_FileClose(h);
}

u8 *COM_LoadFile(const s8 *path, s32 usehunk, u32 *path_id)
{ // Filename are reletive to the quake directory. Allways appends a 0.
	s32 h;
	u8 *buf = NULL;
	s8 base[32];
	s32 len = COM_OpenFile(path, &h, path_id); //look in filesystem or packs
	if(h == -1) return NULL;
	// extract the filename base name for hunk tag
	COM_FileBase(path, base, sizeof(base));
	switch(usehunk) {
	case LOADFILE_HUNK: buf = (u8 *) Hunk_AllocName(len+1, base); break;
	case LOADFILE_TEMPHUNK: buf = (u8 *) Hunk_TempAlloc(len+1); break;
	case LOADFILE_ZONE: buf = (u8 *) Z_Malloc(len+1); break;
	case LOADFILE_CACHE: buf = (u8*)Cache_Alloc(loadcache,len+1,base);break;
	case LOADFILE_STACK: if(len < loadsize) buf = loadbuf;
			else buf = (u8 *) Hunk_TempAlloc(len+1);
			break;
	case LOADFILE_MALLOC: buf = (u8 *) malloc(len+1); break;
	default: Sys_Error("COM_LoadFile: bad usehunk");
	}
	if(!buf) Sys_Error("COM_LoadFile: not enough space for %s", path);
	((u8 *)buf)[len] = 0;
	s32 nread = Sys_FileRead(h, buf, len);
	COM_CloseFile(h);
	if(nread != len) Sys_Error("COM_LoadFile: Error reading %s", path);
	return buf;
}

u8 *COM_LoadHunkFile(const s8 *path, u32 *path_id)
{ return COM_LoadFile(path, LOADFILE_HUNK, path_id); }

u8 *COM_LoadZoneFile(const s8 *path, u32 *path_id)
{ return COM_LoadFile(path, LOADFILE_ZONE, path_id); }

u8 *COM_LoadTempFile(const s8 *path, u32 *path_id)
{ return COM_LoadFile(path, LOADFILE_TEMPHUNK, path_id); }

void COM_LoadCacheFile(const s8 *path, struct cache_user_s *cu, u32 *path_id)
{
	loadcache = cu;
	COM_LoadFile(path, LOADFILE_CACHE, path_id);
}


u8 *COM_LoadStackFile(const s8 *path, void *buffer, s32 bufsize, u32 *path_id)
{ // uses temp hunk if larger than bufsize
	loadbuf = (u8 *)buffer;
	loadsize = bufsize;
	u8 *buf = COM_LoadFile(path, LOADFILE_STACK, path_id);
	return buf;
}


u8 *COM_LoadMallocFile(const s8 *path, u32 *path_id) // returns malloc'd memory
{ return COM_LoadFile(path, LOADFILE_MALLOC, path_id); }

u8 *COM_LoadMallocFile_TextMode_OSPath(const s8 *path, s64 *len_out)
{
// ericw -- this is used by Host_Loadgame_f. Translate CRLF to LF on load games,
// othewise multiline messages have a garbage character at the end of each line.
// TODO: could handle in a way that allows loading CRLF savegames on mac/linux
// without the junk characters appearing.
	FILE *f = fopen(path, "rt");
	if(f == NULL) return NULL;
	s64 len = COM_filelength(f);
	if(len < 0) { fclose(f); return NULL; }
	u8 *data = (u8 *) malloc(len + 1);
	if(data == NULL) { fclose(f); return NULL; }
	// (actuallen < len) if CRLF to LF translation was performed
	s64 actuallen = fread(data, 1, len, f);
	if(ferror(f)) {
		fclose(f);
		free(data);
		return NULL;
	}
	data[actuallen] = '\0';
	if(len_out != NULL) *len_out = actuallen;
	fclose(f);
	return data;
}

const s8 *COM_ParseIntNewline(const s8 *buffer, s32 *value)
{
	s32 consumed = 0;
	sscanf(buffer, "%i\n%n", value, &consumed);
	return buffer + consumed;
}

const s8 *COM_ParseFloatNewline(const s8 *buffer, f32 *value)
{
	s32 consumed = 0;
	sscanf(buffer, "%f\n%n", value, &consumed);
	return buffer + consumed;
}

const s8 *COM_ParseStringNewline(const s8 *buffer)
{
	s32 consumed = 0;
	com_token[0] = '\0';
	sscanf(buffer, "%1023s\n%n", com_token, &consumed);
	return buffer + consumed;
}

// Takes an explicit(not game tree related) path to a pak file.
// Loads the header and directory, adding the files at the beginning
// of the list so they override previous pack files.
static pack_t *COM_LoadPackFile(const s8 *packfile)
{
	dpackfile_t info[MAX_FILES_IN_PACK];
	s32 packhandle;
	dpackheader_t header;
	if(Sys_FileOpenRead(packfile, &packhandle) == -1) return NULL;
	if(Sys_FileRead(packhandle,&header,sizeof(header))!=(s32)sizeof(header)
			|| memcmp(header.id, "PACK", 4) != 0)
		Sys_Error("%s is not a packfile", packfile);
	header.dirofs = LittleLong(header.dirofs);
	header.dirlen = LittleLong(header.dirlen);
	s32 numpackfiles = header.dirlen / sizeof(dpackfile_t);
	if(header.dirlen < 0 || header.dirofs < 0)
		Sys_Error("Invalid packfile %s(dirlen: %i, dirofs: %i)",
				packfile, header.dirlen, header.dirofs);
	if(!numpackfiles) {
		Sys_Printf("WARNING: %s has no files, ignored\n", packfile);
		Sys_FileClose(packhandle);
		return NULL;
	}
	if(numpackfiles > MAX_FILES_IN_PACK)
		Sys_Error("%s has %i files", packfile, numpackfiles);
	if(numpackfiles != PAK0_COUNT)
		com_modified = 1; // not the original file
	packfile_t *newf=(packfile_t*)Z_Malloc(numpackfiles*sizeof(packfile_t));
	Sys_FileSeek(packhandle, header.dirofs);
	if(Sys_FileRead(packhandle, info, header.dirlen) != header.dirlen)
		Sys_Error("Error reading %s", packfile);
	u16 crc;
	CRC_Init(&crc); // crc the directory to check for modifications
	for(s32 i = 0; i < header.dirlen; i++)
		CRC_ProcessByte(&crc, ((u8 *)info)[i]);
	if(crc != PAK0_CRC_V106 && crc != PAK0_CRC_V101 && crc != PAK0_CRC_V100)
		com_modified = 1;
	for(s32 i = 0; i < numpackfiles; i++) { // parse the directory
		q_strlcpy(newf[i].name,info[i].name,sizeof(newf[i].name));
		newf[i].filepos = LittleLong(info[i].filepos);
		newf[i].filelen = LittleLong(info[i].filelen);
	}
	pack_t *pack = (pack_t *) Z_Malloc(sizeof(pack_t));
	q_strlcpy(pack->filename, packfile, sizeof(pack->filename));
	pack->handle = packhandle;
	pack->numfiles = numpackfiles;
	pack->files = newf;
	return pack;
}

static void COM_AddGameDirectory(const s8 *base, const s8 *dir)
{
	s8 pakfile[MAX_OSPATH+32];
	bool been_here = 0;
	q_strlcpy(com_gamedir, va("%s/%s", base, dir), sizeof(com_gamedir));
	u32 path_id = com_searchpaths ? com_searchpaths->path_id << 1 : 1U;
	searchpath_t *search = (searchpath_t *) Z_Malloc(sizeof(searchpath_t));
	search->path_id = path_id;
	q_strlcpy(search->filename, com_gamedir, sizeof(search->filename));
	search->next = com_searchpaths;
	com_searchpaths = search;
	for(s32 i=0;;i++){//add any pak files in the format pak0.pak pak1.pak...
		snprintf(pakfile,sizeof(pakfile),"%s/pak%i.pak",com_gamedir,i);
		pack_t *pak = COM_LoadPackFile(pakfile);
		pack_t *qspak;
		if(i != 0 || path_id != 1 || fitzmode) qspak = NULL;
		else {
			bool old = com_modified;
			if(been_here) base = host_parms.userdir;
			q_snprintf(pakfile, sizeof(pakfile),
					"%s/quakespasm.pak", base);
			qspak = COM_LoadPackFile(pakfile);
			com_modified = old;
		}
		if(pak){
			search = (searchpath_t*) Z_Malloc(sizeof(searchpath_t));
			search->path_id = path_id;
			search->pack = pak;
			search->next = com_searchpaths;
			com_searchpaths = search;
		}
		if(qspak){
			search = (searchpath_t*) Z_Malloc(sizeof(searchpath_t));
			search->path_id = path_id;
			search->pack = qspak;
			search->next = com_searchpaths;
			com_searchpaths = search;
		}
		if(!pak) break;
	}
}

void SetWorldPal(s8 *path, s8 *cmappath)
{
	u8 worldpalbuf[768];
	s8 ppath[MAX_OSPATH] = "gfx/";
	q_strlcat(ppath, path, MAX_OSPATH);
	q_strlcat(ppath, ".lmp", MAX_OSPATH);
	FILE *f;
	COM_FOpenFile(ppath, &f, NULL);
	if(!f){Con_Printf("Couldn't load %s\n", ppath); return;}
	if(fread(worldpalbuf, 768, 1, f) != 1)
		{Con_Printf("Failed reading %s\n", ppath); fclose(f); return;}
	fclose(f);
	u8 worldcmapbuf[256*64];
	s8 cmpath[MAX_OSPATH] = "gfx/";
	q_strlcat(cmpath, cmappath, MAX_OSPATH);
	q_strlcat(cmpath, ".lmp", MAX_OSPATH);
	COM_FOpenFile(cmpath, &f, NULL);
	if(!f){Con_Printf("Couldn't load %s\n", cmpath); return;}
	if(fread(worldcmapbuf, 256*64, 1, f) != 1)
		{Con_Printf("Failed reading %s\n", cmpath); fclose(f); return;}
	fclose(f);
	Con_DPrintf("Setting %s %s\n", ppath, cmpath);
	memcpy(worldpal, worldpalbuf, 768);
	q_strlcpy(worldpalname, path, MAX_OSPATH);
	VID_SetPalette(worldpal, screen);
	memcpy(worldcmap, worldcmapbuf, 64*256);
	q_strlcpy(worldcmapname, cmappath, MAX_OSPATH);
	vid.colormap = worldcmap;
	fog_lut_built = lit_lut_initialized = 0;
}

void SetUiPal(s8 *path)
{
	u8 uipalbuf[768];
	s8 ppath[MAX_OSPATH] = "gfx/";
	q_strlcat(ppath, path, MAX_OSPATH);
	q_strlcat(ppath, ".lmp", MAX_OSPATH);
	FILE *f;
	COM_FOpenFile(ppath, &f, NULL);
	if(!f){Con_Printf("Couldn't load %s\n", ppath); return;}
	if(fread(uipalbuf, 768, 1, f) != 1)
		{Con_Printf("Failed reading %s\n", ppath); fclose(f); return;}
	fclose(f);
	Con_DPrintf("Setting %s\n", ppath);
	memcpy(uipal, uipalbuf, 768);
	q_strlcpy(uipalname, path, MAX_OSPATH);
	VID_SetPalette(uipal, screenui);
}

static void COM_WorldPal_f()
{
	if(Cmd_Argc() <= 2){
		Con_Printf("Usage: worldpal [palette] [colormap]\n");
		if (worldpalname[0] && worldcmapname[0])
		    Con_Printf("Current: %s %s\n", worldpalname, worldcmapname);
		return;
	}
	SetWorldPal(Cmd_Argv(1), Cmd_Argv(2));
}

static void COM_UiPal_f()
{
	if(Cmd_Argc() <= 1){
		Con_Printf("Usage: uipal [palette]\n");
		if (uipalname[0]) Con_Printf("Current: %s\n", uipalname);
		return;
	}
	SetUiPal(Cmd_Argv(1));
}

static void COM_Game_f()
{
	if(Cmd_Argc() <= 1){
		Con_Printf("\"game\" is \"%s\"\n", COM_SkipPath(com_gamedir));
		return;
	}
	const s8 *p = Cmd_Argv(1);
	const s8 *p2 = Cmd_Argv(2);
	searchpath_t *search;
	if(!registered.value){ // disable shareware quake
Con_Printf("You must have the registered version to use modified games\n");
		return;
	}
	if(!*p || !strcmp(p, ".") || strstr(p, "..") || strstr(p, "/")
			|| strstr(p, "\\") || strstr(p, ":")) {
	Con_Printf("gamedir should be a single directory name, not a path\n");
		return;
	}

	if(*p2) {
		if(strcmp(p2,"-hipnotic") && strcmp(p2,"-rogue")
				&& strcmp(p2,"-quoth")){
		Con_Printf("invalid mission pack argument to \"game\"\n");
			return;
		}
		if(!q_strcasecmp(p, GAMENAME)){
		Con_Printf("no mission pack arguments to %s game\n", GAMENAME);
			return;
		}
	}
	if(Sys_FileType(va("%s/%s", com_basedir, p)) != FS_ENT_DIRECTORY) {
		if(host_parms.userdir == host_parms.basedir ||
				(Sys_FileType(va("%s/%s",host_parms.userdir,p))
				 != FS_ENT_DIRECTORY)) {
			Con_Printf("No such game directory \"%s\"\n", p);
			return;
		}
	}
	if(!q_strcasecmp(p, COM_SkipPath(com_gamedir))){ //no change
		if(com_searchpaths->path_id > 1){ //current game not id1
			if(*p2 && com_searchpaths->path_id == 2){
				// rely on QSpasm treating '-game missionpack'
				// as '-missionpack', otherwise would be a mess
				if(!q_strcasecmp(p, &p2[1])) goto _same;
	Con_Printf("reloading game \"%s\" with \"%s\" support\n", p, &p2[1]);
			}
			else if(!*p2 && com_searchpaths->path_id > 2)
	Con_Printf("reloading game \"%s\" without mission pack support\n", p);
			else goto _same;
		}
		else { _same:
	Con_Printf("\"game\" is already \"%s\"\n", COM_SkipPath(com_gamedir));
			return;
		}
	}
	com_modified = 1;
	CL_Disconnect(); // Kill the server
	Host_ShutdownServer(1);
	Host_WriteConfiguration();
	while(com_searchpaths!=com_base_searchpaths){//Kill extra game if loaded
		if(com_searchpaths->pack) {
			Sys_FileClose(com_searchpaths->pack->handle);
			Z_Free(com_searchpaths->pack->files);
			Z_Free(com_searchpaths->pack);
		}
		search = com_searchpaths->next;
		Z_Free(com_searchpaths);
		com_searchpaths = search;
	}
	hipnotic = 0;
	rogue = 0;
	standard_quake = 1;
	if(q_strcasecmp(p, GAMENAME)) { //game is not id1
		if(*p2){
			COM_AddGameDirectory(com_basedir, &p2[1]);
			standard_quake = 0;
			if(!strcmp(p2,"-hipnotic") || !strcmp(p2,"-quoth"))
				hipnotic = 1;
			else if(!strcmp(p2,"-rogue"))
				rogue = 1;
			if(q_strcasecmp(p, &p2[1])) //don't load twice
				COM_AddGameDirectory(com_basedir, p);
		} else {
			COM_AddGameDirectory(com_basedir, p);
			if(!q_strcasecmp(p,"hipnotic") // -game missionpack ==
				||!q_strcasecmp(p,"quoth")){ // -missionpack
				hipnotic = 1;
				standard_quake = 0;
			} else if(!q_strcasecmp(p,"rogue")){
				rogue = 1;
				standard_quake = 0;
			}
		}
	}
	else { // just update com_gamedir
		q_snprintf(com_gamedir, sizeof(com_gamedir), "%s/%s",
				(host_parms.userdir != host_parms.basedir)?
				(s8 *)host_parms.userdir : com_basedir,
				GAMENAME);
	}
	Cache_Flush();
	Mod_ResetAll();
	if(!isDedicated) {
		W_LoadWadFile();
		Draw_Init();
		SCR_Init();
		Sbar_Init();
	}
	Con_Printf("\"game\" changed to \"%s\"\n", COM_SkipPath(com_gamedir));
	Cbuf_AddText("exec quake.rc\n");
	Cbuf_AddText("vid_unlock\n");
}

void COM_InitFilesystem()
{
	Cvar_RegisterVariable(&registered);
	Cvar_RegisterVariable(&cmdline);
	Cmd_AddCommand("path", COM_Path_f);
	Cmd_AddCommand("game", COM_Game_f);
	Cmd_AddCommand("worldpal", COM_WorldPal_f);
	Cmd_AddCommand("uipal", COM_UiPal_f);
	s32 i = COM_CheckParm("-basedir");
	if(i && i < com_argc-1)
		q_strlcpy(com_basedir, com_argv[i + 1], sizeof(com_basedir));
	else q_strlcpy(com_basedir, host_parms.basedir, sizeof(com_basedir));
	s32 j = strlen(com_basedir);
	if(j < 1) Sys_Error("Bad argument to -basedir");
	if((com_basedir[j-1]=='\\')||(com_basedir[j-1]=='/'))com_basedir[j-1]=0;
	// start up with GAMENAME by default(id1)
	COM_AddGameDirectory(com_basedir, GAMENAME);
	// this is the end of our base searchpath: any set gamedirs, such as
	// those from -game command line arguments or by the 'game' console
	// command will be freed up to here upon a new game command.
	com_base_searchpaths = com_searchpaths;
	if(COM_CheckParm("-rogue"))
		COM_AddGameDirectory(com_basedir, "rogue");
	if(COM_CheckParm("-hipnotic"))
		COM_AddGameDirectory(com_basedir, "hipnotic");
	if(COM_CheckParm("-quoth"))
		COM_AddGameDirectory(com_basedir, "quoth");
	i = COM_CheckParm("-game");
	if(i && i < com_argc-1) {
		const s8 *p = com_argv[i + 1];
		if(!*p || !strcmp(p, ".") || strstr(p, "..") || strstr(p, "/")
				|| strstr(p, "\\") || strstr(p, ":"))
	Sys_Error("gamedir should be a single directory name, not a path\n");
		com_modified = 1;
		// don't load mission packs twice
		if(COM_CheckParm("-rogue") && !q_strcasecmp(p, "rogue"))
			p = NULL;
		if(p && COM_CheckParm("-hipnotic")&&!q_strcasecmp(p,"hipnotic"))
			p = NULL;
		if(p && COM_CheckParm("-quoth") && !q_strcasecmp(p, "quoth")) 
			p = NULL;
		if(p != NULL){
			COM_AddGameDirectory(com_basedir, p);
			// QSpasm extension: treat '-game mpack' as '-mpack'
			if(!q_strcasecmp(p,"rogue")){
				rogue = 1;
				standard_quake = 0;
			}
			if(!q_strcasecmp(p,"hipnotic") ||
					!q_strcasecmp(p,"quoth")){
				hipnotic = 1;
				standard_quake = 0;
			}
		}
	}
	COM_CheckRegistered();
}

// The following FS_*() stdio replacements are needed to perform non-sequential
// reads on files reopened on pak files because we need the bookkeeping about
// file start/end positions. Allocating and filling in the fshandle_t structure
// is the users' responsibility when the file is initially opened.
size_t FS_fread(void *ptr, size_t size, size_t nmemb, fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return 0; }
	if(!ptr){ errno = EFAULT; return 0; }
	if(!size || !nmemb){ // no error, just zero bytes wanted
		errno = 0;
		return 0;
	}
	s64 byte_size = nmemb * size;
	if(byte_size > fh->length - fh->pos) // just read to end
		byte_size = fh->length - fh->pos;
	s64 bytes_read = fread(ptr, 1, byte_size, fh->file);
	fh->pos += bytes_read;
	//fread() must return number of elements read, not total number of bytes
	s64 nmemb_read = bytes_read / size;
	// even if the last member is only read partially
	// it is counted as a whole in the return value.
	if(bytes_read % size) nmemb_read++;
	return nmemb_read;
}

s32 FS_fseek(fshandle_t *fh, s64 offset, s32 whence)
{
	if(!fh){ errno = EBADF; return -1; }
	switch(whence) {
	case SEEK_SET: break;
	case SEEK_CUR: offset += fh->pos; break;
	case SEEK_END: offset = fh->length + offset; break;
	default: errno = EINVAL; return -1;
	}
	if(offset < 0){ errno = EINVAL; return -1; }
	if(offset > fh->length) offset = fh->length; // just seek to end
	s32 ret = fseek(fh->file, fh->start + offset, SEEK_SET);
	if(ret < 0) return ret;
	fh->pos = offset;
	return 0;
}

s32 FS_fclose(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return -1; }
	return fclose(fh->file);
}

s64 FS_ftell(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return -1; }
	return fh->pos;
}

void FS_rewind(fshandle_t *fh)
{
	if(!fh) return;
	clearerr(fh->file);
	fseek(fh->file, fh->start, SEEK_SET);
	fh->pos = 0;
}

s32 FS_feof(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return -1; }
	if(fh->pos >= fh->length) return -1;
	return 0;
}

s32 FS_ferror(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return -1; }
	return ferror(fh->file);
}

s32 FS_fgetc(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return EOF; }
	if(fh->pos >= fh->length) return EOF;
	fh->pos += 1;
	return fgetc(fh->file);
}

s8 *FS_fgets(s8 *s, s32 size, fshandle_t *fh)
{
	if(FS_feof(fh)) return NULL;
	if(size > (fh->length - fh->pos) + 1)
		size = (fh->length - fh->pos) + 1;
	s8 *ret = fgets(s, size, fh->file);
	fh->pos = ftell(fh->file) - fh->start;
	return ret;
}

s64 FS_filelength(fshandle_t *fh)
{
	if(!fh){ errno = EBADF; return -1; }
	return fh->length;
}


unsigned COM_HashString(const s8 *str)
{ // Computes the FNV-1a hash of string str
	unsigned hash = 0x811c9dc5u;
	while(*str) { hash ^= *str++; hash *= 0x01000193u; }
	return hash;
}


unsigned COM_HashBlock(const void *data, size_t size)
{ // Computes the FNV-1a hash of a memory block
	const u8 *ptr = (const u8 *)data;
	unsigned hash = 0x811c9dc5u;
	while(size--) { hash ^= *ptr++; hash *= 0x01000193u; }
	return hash;
}

size_t mz_zip_file_read_func(void *opaque, mz_uint64 ofs, void *buf, size_t n)
{
	if(SDL_SeekIO((SDL_IOStream*)opaque, (Sint64)ofs, SDL_IO_SEEK_SET)<0) return 0;
	return SDL_ReadIO((SDL_IOStream*)opaque, buf, n);
}

bool LOC_LoadFile (const s8 *file)
{
	s8 path[1024];
	s32 i,lineno,warnings;
	s8 *cursor;
	SDL_IOStream *rw = NULL;
	Sint64 sz;
	mz_zip_archive archive;
	size_t size = 0;
	if (localization.text) { // clear existing data
		free(localization.text);
		localization.text = NULL;
	}
	localization.numentries = 0;
	localization.numindices = 0;
	if (!file || !*file) return 0;
	memset(&archive, 0, sizeof(archive));
	q_snprintf(path, sizeof(path), "%s", file);
	rw = SDL_IOFromFile(path, "rb");
	if (!rw) {
		q_snprintf(path, sizeof(path), "%s/QuakeEX.kpf", com_basedir);
		rw = SDL_IOFromFile(path, "rb");
		if (!rw) goto fail;
		sz = SDL_GetIOSize(rw);
		if (sz <= 0) goto fail;
		archive.m_pRead = mz_zip_file_read_func;
		archive.m_pIO_opaque = rw;
		if (!mz_zip_reader_init(&archive, sz, 0)) goto fail;
		localization.text = (s8 *)
		   mz_zip_reader_extract_file_to_heap(&archive, file, &size, 0);
		if (!localization.text) goto fail;
		mz_zip_reader_end(&archive);
		SDL_CloseIO(rw);
		localization.text = (s8 *) realloc(localization.text, size+1);
		localization.text[size] = 0;
	}
	else {
		sz = SDL_GetIOSize(rw);
		if (sz <= 0) goto fail;
		localization.text = (s8 *) calloc(1, sz+1);
		if (!localization.text) {
fail:			mz_zip_reader_end(&archive);
			if (rw) SDL_CloseIO(rw);
			Con_Printf("Couldn't load '%s'\n", file);
			return 0;
		}
		SDL_ReadIO(rw, localization.text, sz);
		SDL_CloseIO(rw);
	}
	cursor = localization.text;
	if((u8)(cursor[0])==0xEF&&(u8)(cursor[1])==0xBB&&(u8)(cursor[2])==0xBF)
		cursor += 3; // skip BOM
	warnings = 0;
	lineno = 0;
	while (*cursor) {
		s8 *line, *equals;
		lineno++;
		while (q_isblank(*cursor)) ++cursor; // skip leading whitespace
		line = cursor;
		equals = NULL;
		// find line end and first equals sign, if any
		while (*cursor && *cursor != '\n') {
			if (*cursor == '=' && !equals)
				equals = cursor;
			cursor++;
		}
		if (line[0] == '/') {
			if (line[1] != '/')
	Con_DPrintf("LOC_LoadFile: malformed comment on line %d\n", lineno);
		}
		else if (equals) {
			s8 *key_end = equals;
			bool leading_quote;
			bool trailing_quote;
			locentry_t *entry;
			s8 *value_src;
			s8 *value_dst;
			s8 *value;
			// trim whitespace before equals sign
			while (key_end != line && q_isspace(key_end[-1]))
				key_end--;
			*key_end = 0;
			value = equals + 1;
			// skip whitespace after equals sign
			while (value != cursor && q_isspace(*value))
				value++;
			leading_quote = (*value == '\"');
			trailing_quote = 0;
			value += leading_quote;
			// transform escape sequences in-place
			value_src = value;
			value_dst = value;
			while (value_src != cursor) {
				if(*value_src=='\\' && value_src+1!=cursor){
					s8 c = value_src[1];
					value_src += 2;
					switch (c) {
					case 'n': *value_dst++ = '\n'; break;
					case 't': *value_dst++ = '\t'; break;
					case 'v': *value_dst++ = '\v'; break;
					case 'b': *value_dst++ = '\b'; break;
					case 'f': *value_dst++ = '\f'; break;
					case '"':
					case '\'':
						*value_dst++ = c;
						break;
					default:
Con_Printf("LOC_LoadFile: unrecognized escape sequence \\%c on line %d\n",
		c, lineno);
						*value_dst++ = c;
						break;
					}
					continue;
				}
				if (*value_src == '\"') {
					trailing_quote = 1;
					*value_dst = 0;
					break;
				}
				*value_dst++ = *value_src++;
			}
			// if not a quoted string, trim trailing whitespace
			if (!trailing_quote) {
				while (value_dst != value
					&& q_isblank(value_dst[-1])) {
					*value_dst = 0;
					value_dst--;
				}
			}
			if(localization.numentries==localization.maxnumentries){
				// grow by 50%
				localization.maxnumentries +=
					localization.maxnumentries >> 1;
				localization.maxnumentries =
					q_max(localization.maxnumentries, 32);
				localization.entries = (locentry_t*)
					realloc(localization.entries,
					sizeof(*localization.entries)
					* localization.maxnumentries);
			}
			entry=&localization.entries[localization.numentries++];
			entry->key = line;
			entry->value = value;
		}
		if (*cursor) *cursor++ = 0; //terminate line and advance to next
	}
	if (developer.value >= 2.f && warnings > 0)
		Sys_Printf("%d strings with unprintable characters\n",warnings);
	// hash all entries
	localization.numindices = localization.numentries * 2; //50% load factor
	if (localization.numindices == 0) {
		Con_Printf("No localized strings in file '%s'\n", file);
		return 0;
	}
	localization.indices = (u32*) realloc(localization.indices,
		localization.numindices * sizeof(*localization.indices));
	memset(localization.indices, 0, localization.numindices
			* sizeof(*localization.indices));
	for (i = 0; i < localization.numentries; i++) {
		locentry_t *entry = &localization.entries[i];
		u32 pos = COM_HashString(entry->key)
			% localization.numindices, end = pos;
		for (;;) {
			if (!localization.indices[pos]) {
				localization.indices[pos] = i + 1;
				break;
			}
			++pos;
			if ((s32)pos == localization.numindices) pos = 0;
			if (pos == end) Sys_Error("LOC_LoadFile failed");
		}
	}
	Con_SafePrintf ("Loaded %d strings from '%s'\n",
			localization.numentries, path);
	return 1;
}

void LOC_Init() { LOC_LoadFile("localization/loc_english.txt"); }

void LOC_Shutdown()
{
	free(localization.indices);
	free(localization.entries);
	free(localization.text);
}


const s8* LOC_GetRawString(const s8 *key)
{ // Returns localized string if available, or NULL otherwise
	if(!localization.numindices || !key || !*key || *key != '$')return NULL;
	key++;
	u32 pos = COM_HashString(key) % localization.numindices;
	u32 end = pos;
	do {
		unsigned idx = localization.indices[pos];
		locentry_t *entry;
		if(!idx) return NULL;
		entry = &localization.entries[idx - 1];
		if(!Q_strcmp(entry->key, key)) return entry->value;
		++pos;
		if((s32)pos == localization.numindices) pos = 0;
	} while(pos != end);
	return NULL;
}


const s8* LOC_GetString(const s8 *key)
{ // Returns localized string if available, or input string otherwise
	const s8* value = LOC_GetRawString(key);
	return value ? value : key;
}

// Returns argument index(>= 0) and advances the string if it starts with a
// placeholder({} or {N}), otherwise returns a negative value and leaves the
// pointer unchanged
static s32 LOC_ParseArg(const s8 **pstr)
{
	const s8 *str = *pstr;
	if(*str != '{') return -1; // opening brace
	++str;
	s32 arg = 0; // optional index, defaulting to 0
	while(q_isdigit(*str)) arg = arg * 10 + *str++ - '0';
	if(*str != '}') return -1; // closing brace
	*pstr = ++str;
	return arg;
}

bool LOC_HasPlaceholders(const s8 *str)
{
	if(!localization.numindices) return 0;
	while(*str) {
		if(LOC_ParseArg(&str) >= 0) return 1;
		str++;
	}
	return 0;
}

// Replaces placeholders(of the form {} or {N}) with the corresponding arguments
// Returns number of written chars, excluding the NUL terminator
// If len > 0, output is always NUL-terminated
size_t LOC_Format(const s8 *format, const s8* (*getarg_fn)
		(s32 idx, void* userdata), void* userdata, s8* out, size_t len)
{
	size_t written = 0;
	s32 numargs = 0;
	if(!len) {
		Con_DPrintf("LOC_Format: no output space\n");
		return 0;
	}
	--len; // reserve space for the terminator
	while(*format && written < len) {
		const s8* insert;
		size_t space_left;
		size_t insert_len;
		s32 argindex = LOC_ParseArg(&format);
		if(argindex < 0) {
			out[written++] = *format++;
			continue;
		}
		insert = getarg_fn(argindex, userdata);
		space_left = len - written;
		insert_len = Q_strlen(insert);
		if(insert_len > space_left) {
		  Con_DPrintf("LOC_Format: overflow at argument #%d\n",numargs);
			insert_len = space_left;
		}
		Q_memcpy(out + written, insert, insert_len);
		written += insert_len;
	}
	if(*format) Con_DPrintf("LOC_Format: overflow\n");
	out[written] = 0;
	return written;
}
