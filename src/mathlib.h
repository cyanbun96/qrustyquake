// CyanBun96: compile with -flto to force the compiler to optimize this
// atsb: windows and linux compatible macro definitions
#if defined(_WIN32)
#define SIL static __forceinline // doesn't work with Makefile.w64? TODO add more ifdefs
#define RESTRICT
#define _USE_MATH_DEFINES
#include <math.h>
#else
#define SIL static inline __attribute__((always_inline))
#define RESTRICT __restrict
#endif
SIL f32 DotProduct(const f32 *restrict a, const f32 *restrict b)
	{f32 sum=0.0f;for(s32 i=0;i<3;++i)sum+=a[i]*b[i];return sum;}
SIL f64 DotProductF64(const f32 *a, const f32 *b)
	{f64 sum=0.0f;for(s32 i=0;i<3;++i)sum+=a[i]*b[i];return sum;}
SIL f32 DotProductU8(const u8 *a, const f32 *b)
	{f32 sum=0.0f;for(s32 i=0;i<3;++i)sum+=a[i]*b[i];return sum;}
SIL void VectorSubtract(const f32 *a, const f32 *b, f32 *c)
	{c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2];}
SIL void VectorAdd(const f32 *a, const f32 *b, f32 *c)
	{c[0] = a[0] + b[0]; c[1] = a[1] + b[1]; c[2] = a[2] + b[2];}
SIL void VectorCopy(const f32 *a, f32 *b)
	{b[0] = a[0]; b[1] = a[1]; b[2] = a[2];}
SIL f32 anglemod(f32 a)
	{a = (360.0 / 65536) * ((s32)(a * (65536 / 360.0)) & 65535); return a;}
SIL void VectorInverse(vec3_t v)
	{v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2];}
SIL void VectorScale(vec3_t in, vec_t scale, vec3_t out)
	{out[0] = in[0]*scale; out[1] = in[1]*scale; out[2] = in[2]*scale;}
SIL vec_t VectorLength(const vec3_t v)
	{return sqrt(DotProduct(v,v)); }
SIL s32 GreatestCommonDivisor(s32 a, s32 b)
	{while (b != 0) { s32 t = b; b = a % b; a = t; } return a;}

SIL void TransformVector(vec3_t in, vec3_t out)
{
	out[0] = DotProduct(in, vright);
	out[1] = DotProduct(in, vup);
	out[2] = DotProduct(in, vpn);
}


SIL f32 VectorNormalize(vec3_t v)
{
	f32 length = VectorLength(v);
	if (length) {
		f32 ilength = 1 / length;
		v[0] *= ilength;
		v[1] *= ilength;
		v[2] *= ilength;
	}
	return length;
}

SIL void ProjectPointOnPlane(vec3_t dst, const vec3_t p, const vec3_t normal)
{
	vec3_t n;
	f32 inv_denom = 1.0F / DotProduct(normal, normal);
	f32 d = DotProduct(normal, p) * inv_denom;
	n[0] = normal[0] * inv_denom;
	n[1] = normal[1] * inv_denom;
	n[2] = normal[2] * inv_denom;
	dst[0] = p[0] - d * n[0];
	dst[1] = p[1] - d * n[1];
	dst[2] = p[2] - d * n[2];
}

SIL void VectorMA(vec3_t veca, f32 scale, vec3_t vecb, vec3_t vecc)
{
	vecc[0] = veca[0] + scale * vecb[0];
	vecc[1] = veca[1] + scale * vecb[1];
	vecc[2] = veca[2] + scale * vecb[2];
}

SIL void CrossProduct(vec3_t v1, vec3_t v2, vec3_t cross)
{
	cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
	cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
	cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

SIL void R_ConcatRotations(f32 in1[3][3], f32 in2[3][3], f32 out[3][3])
{
	out[0][0]=in1[0][0]*in2[0][0]+in1[0][1]*in2[1][0]+in1[0][2]*in2[2][0];
	out[0][1]=in1[0][0]*in2[0][1]+in1[0][1]*in2[1][1]+in1[0][2]*in2[2][1];
	out[0][2]=in1[0][0]*in2[0][2]+in1[0][1]*in2[1][2]+in1[0][2]*in2[2][2];
	out[1][0]=in1[1][0]*in2[0][0]+in1[1][1]*in2[1][0]+in1[1][2]*in2[2][0];
	out[1][1]=in1[1][0]*in2[0][1]+in1[1][1]*in2[1][1]+in1[1][2]*in2[2][1];
	out[1][2]=in1[1][0]*in2[0][2]+in1[1][1]*in2[1][2]+in1[1][2]*in2[2][2];
	out[2][0]=in1[2][0]*in2[0][0]+in1[2][1]*in2[1][0]+in1[2][2]*in2[2][0];
	out[2][1]=in1[2][0]*in2[0][1]+in1[2][1]*in2[1][1]+in1[2][2]*in2[2][1];
	out[2][2]=in1[2][0]*in2[0][2]+in1[2][1]*in2[1][2]+in1[2][2]*in2[2][2];
}

SIL void R_ConcatTransforms(f32 in1[3][4], f32 in2[3][4], f32 out[3][4])
{
out[0][0]=in1[0][0]*in2[0][0]+in1[0][1]*in2[1][0]+in1[0][2]*in2[2][0];
out[0][1]=in1[0][0]*in2[0][1]+in1[0][1]*in2[1][1]+in1[0][2]*in2[2][1];
out[0][2]=in1[0][0]*in2[0][2]+in1[0][1]*in2[1][2]+in1[0][2]*in2[2][2];
out[0][3]=in1[0][0]*in2[0][3]+in1[0][1]*in2[1][3]+in1[0][2]*in2[2][3]+in1[0][3];
out[1][0]=in1[1][0]*in2[0][0]+in1[1][1]*in2[1][0]+in1[1][2]*in2[2][0];
out[1][1]=in1[1][0]*in2[0][1]+in1[1][1]*in2[1][1]+in1[1][2]*in2[2][1];
out[1][2]=in1[1][0]*in2[0][2]+in1[1][1]*in2[1][2]+in1[1][2]*in2[2][2];
out[1][3]=in1[1][0]*in2[0][3]+in1[1][1]*in2[1][3]+in1[1][2]*in2[2][3]+in1[1][3];
out[2][0]=in1[2][0]*in2[0][0]+in1[2][1]*in2[1][0]+in1[2][2]*in2[2][0];
out[2][1]=in1[2][0]*in2[0][1]+in1[2][1]*in2[1][1]+in1[2][2]*in2[2][1];
out[2][2]=in1[2][0]*in2[0][2]+in1[2][1]*in2[1][2]+in1[2][2]*in2[2][2];
out[2][3]=in1[2][0]*in2[0][3]+in1[2][1]*in2[1][3]+in1[2][2]*in2[2][3]+in1[2][3];
}

SIL s32 BoxOnPlaneSide(vec3_t emins, vec3_t emaxs, mplane_t *p)
{ // Returns 1, 2, or 1 + 2
        if(p->type < 3){ //lordhavoc
                if (p->dist <= emins[p->type]) return 1;
		else if(p->dist >= emaxs[p->type]) return 2;
		return 3;
	}
	f32 dist1 = 0, dist2 = 0;
	switch (p->signbits) { // general case
	case 0:
	dist1=p->normal[0]*emaxs[0]+p->normal[1]*emaxs[1]+p->normal[2]*emaxs[2];
	dist2=p->normal[0]*emins[0]+p->normal[1]*emins[1]+p->normal[2]*emins[2];
	break;
	case 1:
	dist1=p->normal[0]*emins[0]+p->normal[1]*emaxs[1]+p->normal[2]*emaxs[2];
	dist2=p->normal[0]*emaxs[0]+p->normal[1]*emins[1]+p->normal[2]*emins[2];
	break;
	case 2:
	dist1=p->normal[0]*emaxs[0]+p->normal[1]*emins[1]+p->normal[2]*emaxs[2];
	dist2=p->normal[0]*emins[0]+p->normal[1]*emaxs[1]+p->normal[2]*emins[2];
	break;
	case 3:
	dist1=p->normal[0]*emins[0]+p->normal[1]*emins[1]+p->normal[2]*emaxs[2];
	dist2=p->normal[0]*emaxs[0]+p->normal[1]*emaxs[1]+p->normal[2]*emins[2];
	break;
	case 4:
	dist1=p->normal[0]*emaxs[0]+p->normal[1]*emaxs[1]+p->normal[2]*emins[2];
	dist2=p->normal[0]*emins[0]+p->normal[1]*emins[1]+p->normal[2]*emaxs[2];
	break;
	case 5:
	dist1=p->normal[0]*emins[0]+p->normal[1]*emaxs[1]+p->normal[2]*emins[2];
	dist2=p->normal[0]*emaxs[0]+p->normal[1]*emins[1]+p->normal[2]*emaxs[2];
	break;
	case 6:
	dist1=p->normal[0]*emaxs[0]+p->normal[1]*emins[1]+p->normal[2]*emins[2];
	dist2=p->normal[0]*emins[0]+p->normal[1]*emaxs[1]+p->normal[2]*emaxs[2];
	break;
	case 7:
	dist1=p->normal[0]*emins[0]+p->normal[1]*emins[1]+p->normal[2]*emins[2];
	dist2=p->normal[0]*emaxs[0]+p->normal[1]*emaxs[1]+p->normal[2]*emaxs[2];
	break;
	default:
		Sys_Error("BoxOnPlaneSide:  Bad signbits");
		break;
	}
	s32 sides = 0;
	if (dist1 >= p->dist) sides = 1;
	if (dist2 < p->dist) sides |= 2;
	return sides;
}

SIL void AngleVectors(vec3_t angles, vec3_t forward, vec3_t right, vec3_t up)
{
	f32 angle = angles[YAW] * (M_PI * 2 / 360);
	f32 sy = sin(angle);
	f32 cy = cos(angle);
	angle = angles[PITCH] * (M_PI * 2 / 360);
	f32 sp = sin(angle);
	f32 cp = cos(angle);
	angle = angles[ROLL] * (M_PI * 2 / 360);
	f32 sr = sin(angle);
	f32 cr = cos(angle);
	forward[0] = cp * cy;
	forward[1] = cp * sy;
	forward[2] = -sp;
	right[0] = (-1 * sr * sp * cy + -1 * cr * -sy);
	right[1] = (-1 * sr * sp * sy + -1 * cr * cy);
	right[2] = -1 * sr * cp;
	up[0] = (cr * sp * cy + -sr * -sy);
	up[1] = (cr * sp * sy + -sr * cy);
	up[2] = cr * cp;
}

// Returns mathematically correct (floor-based) quotient and remainder for
// numer and denom, both of which should contain no fractional part. The
// quotient must fit in 32 bits.
SIL void FloorDivMod(f64 numer, f64 denom, s32 *quotient, s32 *rem)
{
	s32 q, r;
	if (numer >= 0.0) {
		f64 x = floor(numer / denom);
		q = (s32)x;
		r = (s32)floor(numer - (x * denom));
	} else { // perform operations with positive values, and fix mod
		f64 x = floor(-numer / denom); // to make floor-based
		q = -(s32)x;
		r = (s32)floor(-numer - (x * denom));
		if (r != 0) {
			q--;
			r = (s32)denom - r;
		}
	}
	*quotient = q;
	*rem = r;
}
#undef SIL
