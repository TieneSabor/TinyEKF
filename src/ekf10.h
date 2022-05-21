#ifndef _EKF10_H
#define _EKF10_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "stdio.h"
#include "tiny_util.h"

// Function return value.
#define EKF10_OK 0 // Done Successfully.
#define EKF10_ER 1 // Error

void ekf10_init(FTYPE wx, FTYPE wy, FTYPE wz);

void ekf10_set_Propagation_Noise(FTYPE Wwx, FTYPE Wwy, FTYPE Wwz, FTYPE Wbx,
                                 FTYPE Wby, FTYPE Wbz, FTYPE Wq1, FTYPE Wq2,
                                 FTYPE Wq3, FTYPE Wq4);

void ekf10_set_Observation_Noise(FTYPE Vwx, FTYPE Vwy, FTYPE Vwz, FTYPE Vq1,
                                 FTYPE Vq2, FTYPE Vq3, FTYPE Vq4);

/* x_{k+1|k} = f(x_k)
 * P_{k+1|k} = F*P_k*F^T + W
 */
void ekf10_predict(FTYPE dt);

/* K = P*H^T*(H*P*H^T + V)^(-1)
 * x_{k+1} = x_{x+1|x} + K(y - H*x_{K+1|k})
 * P_{k+1} = (1 - K*H)*P*(1 - K*H)^T + K*V*K^T
 */
void ekf10_update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx, FTYPE qy,
                  FTYPE qz);

/*
 */
void ekf10_covariance_predict(FTYPE dt);

/*
 */
void ekf10_innovation_update(void);

/*
 */
void ekf10_innovation_invert(void);

/*  */
void ekf10_KH_update(void);

/*  */
void ekf10_covariance_update(void);

/*  */
void ekf10_state_update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx,
                        FTYPE qy, FTYPE qz);

/*  */
void ekf10_get_angular_vel(FTYPE *wx, FTYPE *wy, FTYPE *wz);

/*  */
void ekf10_get_gyroscope_bias(FTYPE *bx, FTYPE *by, FTYPE *bz);

/*  */
void ekf10_get_quaternion(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void ekf10_get_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

// void printState(void);

#ifdef __cplusplus
}
#endif

#endif
