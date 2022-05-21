#ifndef _QUATERNION_ESTIMATOR_6AXIS_H
#define _QUATERNION_ESTIMATOR_6AXIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "tiny_util.h"

void qe6_set_reference(FTYPE eax, FTYPE eay, FTYPE eaz);

void qe6_set_measurement(FTYPE bax, FTYPE bay, FTYPE baz);

void qe6_quaternion_init(void);

void qe6_quaternion_update(void);

void qe6_get_qk(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void qe6_get_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

#ifdef __cplusplus
}
#endif

#endif