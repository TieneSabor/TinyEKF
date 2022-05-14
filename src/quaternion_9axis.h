#ifndef _QUATERNION_ESTIMATOR_H
#define _QUATERNION_ESTIMATOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "tiny_util.h"

void set_reference(FTYPE eax, FTYPE eay, FTYPE eaz, FTYPE emx, FTYPE emy,
                   FTYPE emz);

void set_measurement(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE bmx, FTYPE bmy,
                     FTYPE bmz);

void quaternion_init(void);

void quaternion_update(void);

void get_qe9_qk(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void get_qe9_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

#ifdef __cplusplus
}
#endif

#endif
