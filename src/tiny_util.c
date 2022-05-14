#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include "tiny_util.h"

// some simple math function
int imin(int a, int b) { return ((a) < (b) ? (a) : (b)); }

int imax(int a, int b) { return ((a) >= (b) ? (a) : (b)); }

FTYPE Fabs(FTYPE r) { return (r > 0 ? r : (-r)); }

FTYPE Fsqrt(FTYPE b) {
  FTYPE err = 1e-6;
  FTYPE x = b / 2, nx = 0.5 * (x + b / x);
  unsigned char cnt = 0;
  while (Fabs(nx - x) > err) {
    if (cnt < 255) {
      x = nx;
      nx = 0.5 * (x + b / x);
    } else {
      break;
    }
    cnt++;
  }
  return nx;
}

FTYPE Fpow(FTYPE x, unsigned int a) {
  FTYPE res = 1;
  for (int i = 0; i < a; i++) {
    res = res * x;
  }
  return res;
}
