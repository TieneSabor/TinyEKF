#ifndef _EKF7_H
#define _EKF7_H

#ifdef __cplusplus
extern "C" {
#endif

#include "math.h"
#include "stdio.h"
#include "tiny_util.h"

// Function return value.
#define EKF7_OK 0 // Done Successfully.
#define EKF7_ER 1 // Error

#define EKF7_COV_UPDATE 1
#define EKF7_COV_NOT_UPDATE 0

typedef struct matrix {
  byte dim; // dimision of the matrix, dim = m<<4 + n
  byte sps; // start position in the_space
} M;

void ekf7_set_state(FTYPE bx, FTYPE by, FTYPE bz, FTYPE qw, FTYPE qx, FTYPE qy,
                    FTYPE qz);

void ekf7_init(void);

void ekf7_set_Propagation_Noise(FTYPE Wbx, FTYPE Wby, FTYPE Wbz, FTYPE Wq1,
                                FTYPE Wq2, FTYPE Wq3, FTYPE Wq4);

void ekf7_set_Observation_Noise(FTYPE Vq1, FTYPE Vq2, FTYPE Vq3, FTYPE Vq4);

/* x_{k+1|k} = f(x_k)
 * P_{k+1|k} = F*P_k*F^T + W
 */
void ekf7_predict(FTYPE dt, FTYPE wx, FTYPE wy, FTYPE wz,
                  char covar_update_flag);

/* K = P*H^T*(H*P*H^T + V)^(-1)
 * x_{k+1} = x_{x+1|x} + K(y - H*x_{K+1|k})
 * P_{k+1} = (1 - K*H)*P*(1 - K*H)^T + K*V*K^T
 */
void ekf7_update(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz);

/*
 */
void ekf7_covariance_predict(FTYPE dt);

/*
 */
void ekf7_innovation_update(void);

/*
 */
void ekf7_innovation_invert(void);

/*  */
void ekf7_KH_update(void);

/*  */
void ekf7_covariance_update(void);

/*  */
void ekf7_state_update(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz);

/*  */
// void ekf7_get_angular_vel(FTYPE *wx, FTYPE *wy, FTYPE *wz);

/*  */
void ekf7_get_gyroscope_bias(FTYPE *bx, FTYPE *by, FTYPE *bz);

/*  */
void ekf7_get_quaternion(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);

void ekf7_get_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw);

// void printState(void);

#ifdef DBG_CODE
void printMat(M x);

void printSyM(M x);

void printF(void);

void printState(void);
#endif

#ifdef __cplusplus
}
#endif

#endif
