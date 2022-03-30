#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include <TinyEKF.hpp>

TinyEKF::TinyEKF(byte Setting, FTYPE wx, FTYPE wy, FTYPE wz) {
  set_ = Setting;
  ekf_init(wx, wy, wz);
  switch (set_) {
  case QTN_ONLY:
    break;
  case AX6_ONLY:
    quaternion_init_QE6();
    break;
  case AX9_ONLY:
    quaternion_init();
    break;
  }
}

void TinyEKF::Set_Propagation_Noise(FTYPE Wwx, FTYPE Wwy, FTYPE Wwz, FTYPE Wbx,
                                    FTYPE Wby, FTYPE Wbz, FTYPE Wq1, FTYPE Wq2,
                                    FTYPE Wq3, FTYPE Wq4) {
  set_Propagation_Noise(Wwx, Wwy, Wwz, Wbx, Wby, Wbz, Wq1, Wq2, Wq3, Wq4);
}

void TinyEKF::Set_Observation_Noise(FTYPE Vwx, FTYPE Vwy, FTYPE Vwz, FTYPE Vq1,
                                    FTYPE Vq2, FTYPE Vq3, FTYPE Vq4) {
  set_Observation_Noise(Vwx, Vwy, Vwz, Vq1, Vq2, Vq3, Vq4);
}

void TinyEKF::Set_References(FTYPE eax, FTYPE eay, FTYPE eaz, FTYPE emx = 0,
                             FTYPE emy = 0, FTYPE emz = 0) {
  switch (set_) {
  case QTN_ONLY:
    break;
  case AX9_ONLY:
    set_reference(eax, eay, eaz, emx, emy, emz);
    break;
  case AX6_ONLY:
    set_reference_QE6(eax, eay, eaz);
    break;
  }
}

void TinyEKF::EKF_predict(FTYPE dt) { predict(dt); }

void TinyEKF::EKF_update(void) {
  switch (set_) {
  case QTN_ONLY:
    update(wx_, wy_, wz_, qw_, qx_, qy_, qz_);
    break;
  case AX9_ONLY:
    set_measurement(bax_, bay_, baz_, bmx_, bmy_, bmz_);
    quaternion_update();
    get_qe_qk(&qw_, &qx_, &qy_, &qz_);
    update(wx_, wy_, wz_, qw_, qx_, qy_, qz_);
    break;
  case AX6_ONLY:
    set_measurement_QE6(bax_, bay_, baz_);
    quaternion_update_QE6();
    get_qk_QE6(&qw_, &qx_, &qy_, &qz_);
    update(wx_, wy_, wz_, qw_, qx_, qy_, qz_);
    break;
  }
}

void TinyEKF::Set_IMU_Measurements(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE wx,
                                   FTYPE wy, FTYPE wz, FTYPE bmx, FTYPE bmy,
                                   FTYPE bmz) {
  bax_ = bax;
  bay_ = bay;
  baz_ = baz;
  bmx_ = bmx;
  bmy_ = bmy;
  bmz_ = bmz;
  wx_ = wx;
  wy_ = wy;
  wz_ = wz;
}

void TinyEKF::Set_QTN_Measurements(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz = 0) {
  qw_ = qw;
  qx_ = qx;
  qy_ = qy;
  qz_ = qz;
}

void TinyEKF::Get_Angular_Velocity(FTYPE *wx, FTYPE *wy, FTYPE *wz) {
  get_angular_vel(wx, wy, wz);
}

void TinyEKF::Get_Gyroscope_Bias(FTYPE *bx, FTYPE *by, FTYPE *bz) {
  get_gyroscope_bias(bx, by, bz);
}

void TinyEKF::Get_Quaternions(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  get_quaternion(qw, qx, qy, qz);
}

void TinyEKF::Get_RPY(FTYPE *yaw, FTYPE *pitch, FTYPE *roll) {
  get_ekf_rpy(yaw, pitch, roll);
}