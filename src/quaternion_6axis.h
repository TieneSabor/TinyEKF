#ifndef _QUATERNION_ESTIMATOR_6AXIS_H
#define _QUATERNION_ESTIMATOR_6AXIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "tiny_util.h"

void set_reference_QE6(FTYPE eax, FTYPE eay, FTYPE eaz);

void set_measurement_QE6(FTYPE bax, FTYPE bay, FTYPE baz);

void quaternion_init_QE6(void);

void quaternion_update_QE6(void);

void get_qk_QE6(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void get_rpy_QE6(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

#ifdef __cplusplus
}
#endif

#endif