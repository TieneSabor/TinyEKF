#ifndef _EKF_H
#define _EKF_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "stdio.h"
#include "tiny_util.h"

// Function return value.
#define EKF_OK 0 // Done Successfully.
#define EKF_ER 1 // Error

// Some variable for tuning.
// Gyroscope latency in sec.
#define TAUX 1
#define TAUY 1
#define TAUZ 1

void ekf_init(FTYPE wx, FTYPE wy, FTYPE wz);

void set_Propagation_Noise(FTYPE Wwx, FTYPE Wwy, FTYPE Wwz, FTYPE Wbx,
                           FTYPE Wby, FTYPE Wbz, FTYPE Wq1, FTYPE Wq2,
                           FTYPE Wq3, FTYPE Wq4);

void set_Observation_Noise(FTYPE Vwx, FTYPE Vwy, FTYPE Vwz, FTYPE Vq1,
                           FTYPE Vq2, FTYPE Vq3, FTYPE Vq4);

/* x_{k+1|k} = f(x_k)
 * P_{k+1|k} = F*P_k*F^T + W
 */
void predict(FTYPE dt);

/* K = P*H^T*(H*P*H^T + V)^(-1)
 * x_{k+1} = x_{x+1|x} + K(y - H*x_{K+1|k})
 * P_{k+1} = (1 - K*H)*P*(1 - K*H)^T + K*V*K^T
 */
void update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx, FTYPE qy,
            FTYPE qz);

/*
 */
void covariance_predict(FTYPE dt);

/*
 */
void innovation_update(void);

/*
 */
void innovation_invert(void);

/*  */
void KH_update(void);

/*  */
void covariance_update(void);

/*  */
void state_update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx, FTYPE qy,
                  FTYPE qz);

/*  */
void get_angular_vel(FTYPE *wx, FTYPE *wy, FTYPE *wz);

/*  */
void get_gyroscope_bias(FTYPE *bx, FTYPE *by, FTYPE *bz);

/*  */
void get_quaternion(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void get_ekf_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

// void printState(void);

#ifdef __cplusplus
}
#endif

#endif
