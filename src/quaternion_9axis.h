#ifndef _QUATERNION_ESTIMATOR_H
#define _QUATERNION_ESTIMATOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "tiny_util.h"

void qe9_set_reference(FTYPE eax, FTYPE eay, FTYPE eaz, FTYPE emx, FTYPE emy,
                       FTYPE emz);

void qe9_set_measurement(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE bmx, FTYPE bmy,
                         FTYPE bmz);

void qe9_quaternion_init(void);

void qe9_quaternion_update(void);

void qe9_get_qk(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void qe9_get_py(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

#ifdef __cplusplus
}
#endif

#endif
