#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include "quaternion_6axis.h"

// Important types used here.
typedef unsigned char byte;
typedef struct matrix {
  byte dim; // dimision of the matrix, dim = m<<4 + n
  byte sps; // start position in the_space
} M;

// Definition of important constant.
// Should NOT be modified.
#define SPACE_SIZE 45 // maximum matrix element number
#define MAX_DIM 15    // maximum matrix operation
#define LAST4 0x0f    // 00001111
#define NUM_RF 3 // reference vector elements number (1 3D vector make it 3)
#define NUM_QT 3 // quaternion element number. We make qz always 0.

// Byte operation. For "dim" in matrix struct, it is composed of M(first 4 bit)
// and N(last 4 bit)
#define GET_M(x) (x >> 4)
#define GET_N(x) (x & 0x0f)
#define SET_B(m, n) ((m << 4) + n)

// The (x, y)th element for NOT symmetric matrix of m rows * n columns, starting
// with (0, 0). With O2 Opt., the add/mul. part in the index should be replaced
// by const. (Hopefully?)
#define MAT(mat, x, y) the_space_[mat.sps + y + x * (mat.dim & 0x0f)]
// Get the (x, y)th element of a symmetric matrix.
#define SYM(mat, x, y)                                                         \
  the_space_[mat.sps + imax(x, y) + imin(x, y) * (mat.dim & 0x0f) -            \
             (imin(x, y) * imin(x, y) + imin(x, y)) / 2]
// Get how much a NOT symmetric matrix take to store data
#define OCU(mat) (mat.dim >> 4) * (mat.dim & 0x0f)
// Get how much a symmetric matrix take to store data
#define OCS(mat)                                                               \
  0.5 * ((mat.dim >> 4) * ((mat.dim & 0x0f) - 1)) + (mat.dim >> 4)

// Every element of matrix in the quaternion estimator will be stored here.
FTYPE the_space_[SPACE_SIZE];
byte spscnt_;

// Matrix used for Quaternion estimation
M J_;   // Jacobian, 6*4
M PD_;  // Positive Definite Matrix of J^T*J, 4*4
M e_;   // Error between earth frame references and measured directions, 6*1
M tmp_; // Temp space for matrix multiplication, 6*1

// Reference vector value (earth frame) and measured vector value (body frame)
FTYPE eax_, eay_, eaz_;
FTYPE bax_, bay_, baz_;
// Quaternion Estimation
FTYPE qw_, qx_, qy_, qz_;

void set_reference_QE6(FTYPE eax, FTYPE eay, FTYPE eaz) {
  eax_ = eax;
  eay_ = eay;
  eaz_ = eaz;
}

void set_measurement_QE6(FTYPE bax, FTYPE bay, FTYPE baz) {
  bax_ = bax;
  bay_ = bay;
  baz_ = baz;
}

void quaternion_init_QE6(void) {
  // Fill all elements with 0
  for (int i = 0; i < SPACE_SIZE; i++) {
    the_space_[i] = 0;
  }
  // Initialize all structures.
  // Jacobian, non sysmmetric 6*4
  spscnt_ = 0;
  J_.dim = SET_B(NUM_RF, NUM_QT);
  J_.sps = spscnt_;
  spscnt_ += OCU(J_);
  // Positive Definite, symmetric 4*4
  PD_.dim = SET_B(NUM_QT, NUM_QT);
  PD_.sps = spscnt_;
  spscnt_ += OCS(PD_);
  // Error, vector 6*1
  e_.dim = SET_B(NUM_RF, 1);
  e_.sps = spscnt_;
  spscnt_ += OCU(e_);
  // Tmp space, vector 6*1
  tmp_.dim = SET_B(NUM_RF, 1);
  tmp_.sps = spscnt_;
  spscnt_ += OCU(tmp_);
  // Initialize q:
  qw_ = 1;
  qx_ = 0;
  qy_ = 0;
  qz_ = 0;
}

void set_J_QE6(void) {
  // qz_ must be 0
  qz_ = 0;
  // Set J by reference vectors and quaternions
  MAT(J_, 0, 0) = -2 * (qw_ * bax_ - qz_ * bay_ + qy_ * baz_);
  MAT(J_, 1, 0) = -2 * (qz_ * bax_ + qw_ * bay_ - qx_ * baz_);
  MAT(J_, 2, 0) = -2 * (-qy_ * bax_ + qx_ * bay_ + qw_ * baz_);

  MAT(J_, 0, 1) = -2 * (qx_ * bax_ + qy_ * bay_ + qz_ * baz_);
  MAT(J_, 1, 1) = -2 * (qy_ * bax_ - qx_ * bay_ - qw_ * baz_);
  MAT(J_, 2, 1) = -2 * (qz_ * bax_ + qw_ * bay_ - qx_ * baz_);

  MAT(J_, 0, 2) = -2 * (-qy_ * bax_ + qx_ * bay_ + qw_ * baz_);
  MAT(J_, 1, 2) = -2 * (qx_ * bax_ + qy_ * bay_ + qz_ * baz_);
  MAT(J_, 2, 2) = -2 * (-qw_ * bax_ + qz_ * bay_ - qy_ * baz_);
  // PD = J^T*J
  for (int i = 0; i < NUM_QT; i++) {
    for (int j = 0; j < NUM_QT; j++) {
      SYM(PD_, i, j) = 0;
      for (int k = 0; k < NUM_RF; k++) {
        SYM(PD_, i, j) += MAT(J_, k, i) * MAT(J_, k, j);
      }
    }
  }
}

void PD_invert_QE6(void) {
  SYM(PD_, 0, 0) = 1.0f / SYM(PD_, 0, 0);
  for (int i = 0; i < NUM_QT - 1; i++) {
    // w^T = -inv(A)*u^T
    // Here, the submatrix up to I(i, i) (that is, A) was already inverted.
    for (int j = 0; j <= i; j++) {
      MAT(tmp_, j, 0) = 0;
      for (int k = 0; k <= i; k++) {
        MAT(tmp_, j, 0) -= SYM(PD_, j, k) * SYM(PD_, k, i + 1);
      }
    }
    // y = inv(u*w^T+x)
    // Firstly u*w^T:
    FTYPE y = 0;
    for (int j = 0; j <= i; j++) {
      y += SYM(PD_, j, i + 1) * MAT(tmp_, j, 0);
    }
    SYM(PD_, i + 1, i + 1) = 1.0f / (y + SYM(PD_, i + 1, i + 1));
    // v^T = w^T * y
    for (int j = 0; j <= i; j++) {
      SYM(PD_, j, i + 1) = MAT(tmp_, j, 0) * SYM(PD_, i + 1, i + 1);
    }
    // B = inv(A) + w^T*v
    for (int j = 0; j <= i; j++) {
      for (int k = j; k <= i; k++) {
        SYM(PD_, j, k) = SYM(PD_, j, k) + MAT(tmp_, j, 0) * SYM(PD_, k, i + 1);
      }
    }
  }
  return;
}

void get_error_QE6(void) {
  FTYPE R11 = pow(qw_, 2) + pow(qx_, 2) - pow(qy_, 2) - pow(qz_, 2);
  FTYPE R12 = 2 * (qx_ * qy_ - qw_ * qz_);
  FTYPE R13 = 2 * (qx_ * qz_ + qw_ * qy_);
  FTYPE R21 = 2 * (qx_ * qy_ + qw_ * qz_);
  FTYPE R22 = pow(qw_, 2) + pow(qy_, 2) - pow(qx_, 2) - pow(qz_, 2);
  FTYPE R23 = 2 * (qy_ * qz_ - qw_ * qx_);
  FTYPE R31 = 2 * (qx_ * qz_ - qw_ * qy_);
  FTYPE R32 = 2 * (qy_ * qz_ + qw_ * qx_);
  FTYPE R33 = pow(qw_, 2) + pow(qz_, 2) - pow(qy_, 2) - pow(qx_, 2);

  MAT(e_, 0, 0) = eax_ - R11 * bax_ - R12 * bay_ - R13 * baz_;
  MAT(e_, 1, 0) = eay_ - R21 * bax_ - R22 * bay_ - R23 * baz_;
  MAT(e_, 2, 0) = eaz_ - R31 * bax_ - R32 * bay_ - R33 * baz_;
}

void quaternion_update_QE6(void) {
  set_J_QE6();
  PD_invert_QE6();
  get_error_QE6();
  // get dq = PD*J^T*e
  for (int i = 0; i < NUM_QT; i++) {
    // dq_i
    MAT(tmp_, i, 0) = 0;
    for (int j = 0; j < NUM_QT; j++) {
      for (int k = 0; k < NUM_RF; k++) {
        MAT(tmp_, i, 0) += SYM(PD_, i, j) * MAT(J_, k, j) * MAT(e_, k, 0);
      }
    }
  }
  // q_k+1 = q_k + dq
  qw_ -= MAT(tmp_, 0, 0);
  qx_ -= MAT(tmp_, 1, 0);
  qy_ -= MAT(tmp_, 2, 0);
  qz_ = 0;
}

void get_qk_QE6(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  *qw = qw_;
  *qx = qx_;
  *qy = qy_;
  *qz = qz_;
}

void get_rpy_QE6(FTYPE *yaw, FTYPE *pitch, FTYPE *roll) {
  *yaw = atan2(2 * (qw_ * qx_ + qy_ * qz_), 1 - 2 * (qx_ * qx_ + qy_ * qy_));
  *pitch = asin(2 * (qw_ * qy_ - qz_ * qx_));
  *roll = atan2(2 * (qw_ * qz_ + qx_ * qy_), 1 - 2 * (qy_ * qy_ + qz_ * qz_));
}
