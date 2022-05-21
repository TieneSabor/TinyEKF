#ifndef _TINYEKF_HPP
#define _TINYEKF_HPP

#include "math.h"

#include <ekf10.h>
#include <quaternion_6axis.h>
#include <quaternion_9axis.h>

// Settings
#define QTN_ONLY 1
#define AX6_ONLY 2
#define AX9_ONLY 3

class TinyEKF {
public:
  TinyEKF(byte Setting, FTYPE wx = 0, FTYPE wy = 0, FTYPE wz = 0);
  // Set up parameters
  void Set_Propagation_Noise(FTYPE Wwx, FTYPE Wwy, FTYPE Wwz, FTYPE Wbx,
                             FTYPE Wby, FTYPE Wbz, FTYPE Wq1, FTYPE Wq2,
                             FTYPE Wq3, FTYPE Wq4);
  void Set_Observation_Noise(FTYPE Vwx, FTYPE Vwy, FTYPE Vwz, FTYPE Vq1,
                             FTYPE Vq2, FTYPE Vq3, FTYPE Vq4);
  void Set_References(FTYPE eax, FTYPE eay, FTYPE eaz, FTYPE emx = 0,
                      FTYPE emy = 1, FTYPE emz = 0);
  // Main operations
  void EKF_predict(FTYPE dt);
  void EKF_update(void);
  // Set up Inputs
  void Set_IMU_Measurements(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE wx, FTYPE wy,
                            FTYPE wz, FTYPE bmx = 0, FTYPE bmy = 1,
                            FTYPE bmz = 0);
  void Set_QTN_Measurements(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz = 0);
  // Get Outputs
  void Get_Angular_Velocity(FTYPE *wx, FTYPE *wy, FTYPE *wz);
  void Get_Gyroscope_Bias(FTYPE *bx, FTYPE *by, FTYPE *bz);
  void Get_Quaternions(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz);
  void Get_RPY(FTYPE *yaw, FTYPE *pitch, FTYPE *roll);

private:
  byte set_ = 0;
  FTYPE bax_, bay_, baz_, bmx_, bmy_, bmz_;
  FTYPE wx_, wy_, wz_;
  FTYPE qw_, qx_, qy_, qz_;
};

#endif