#ifndef _TINY_UTIL_H
#define _TINY_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

// Type for float point operation.
// Can be float or double
#define FTYPE double

typedef unsigned char byte;

// some simple math function
int imin(int a, int b);

int imax(int a, int b);

FTYPE Fabs(FTYPE r);

FTYPE Fsqrt(FTYPE b);

FTYPE Fpow(FTYPE x, unsigned int a);

#ifdef __cplusplus
}
#endif

#endif
