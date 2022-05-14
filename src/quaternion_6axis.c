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
  byte sps; // start position in thqe6_e_space
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
#define MAT(mat, x, y) qe6_spacqe6_e_[mat.sps + y + x * (mat.dim & 0x0f)]
// Get the (x, y)th element of a symmetric matrix.
#define SYM(mat, x, y)                                                         \
  qe6_spacqe6_e_[mat.sps + imax(x, y) + imin(x, y) * (mat.dim & 0x0f) -        \
                 (imin(x, y) * imin(x, y) + imin(x, y)) / 2]
// Get how much a NOT symmetric matrix take to store data
#define OCU(mat) (mat.dim >> 4) * (mat.dim & 0x0f)
// Get how much a symmetric matrix take to store data
#define OCS(mat)                                                               \
  0.5 * ((mat.dim >> 4) * ((mat.dim & 0x0f) - 1)) + (mat.dim >> 4)

// Every element of matrix in the quaternion estimator will be stored here.
FTYPE qe6_spacqe6_e_[SPACE_SIZE];
byte qe6_spscnt_;

// Matrix used for Quaternion estimation
M qe6_J_;   // Jacobian, 6*4
M qe6_PD_;  // Positive Definite Matrix of J^T*J, 4*4
M qe6_e_;   // Error between earth frame references and measured directions, 6*1
M qe6_tmp_; // Temp space for matrix multiplication, 6*1

// Reference vector value (earth frame) and measured vector value (body frame)
FTYPE qe6_eax_, qe6_eay_, qe6_eaz_;
FTYPE qe6_bax_, qe6_bay_, qe6_baz_;
// Quaternion Estimation
FTYPE qe6_qw_, qe6_qx_, qe6_qy_, qe6_qz_;

void set_reference_QE6(FTYPE eax, FTYPE eay, FTYPE eaz) {
  qe6_eax_ = eax;
  qe6_eay_ = eay;
  qe6_eaz_ = eaz;
}

void set_measurement_QE6(FTYPE bax, FTYPE bay, FTYPE baz) {
  qe6_bax_ = bax;
  qe6_bay_ = bay;
  qe6_baz_ = baz;
}

void quaternion_init_QE6(void) {
  // Fill all elements with 0
  for (int i = 0; i < SPACE_SIZE; i++) {
    qe6_spacqe6_e_[i] = 0;
  }
  // Initialize all structures.
  // Jacobian, non sysmmetric 6*4
  qe6_spscnt_ = 0;
  qe6_J_.dim = SET_B(NUM_RF, NUM_QT);
  qe6_J_.sps = qe6_spscnt_;
  qe6_spscnt_ += OCU(qe6_J_);
  // Positive Definite, symmetric 4*4
  qe6_PD_.dim = SET_B(NUM_QT, NUM_QT);
  qe6_PD_.sps = qe6_spscnt_;
  qe6_spscnt_ += OCS(qe6_PD_);
  // Error, vector 6*1
  qe6_e_.dim = SET_B(NUM_RF, 1);
  qe6_e_.sps = qe6_spscnt_;
  qe6_spscnt_ += OCU(qe6_e_);
  // Tmp space, vector 6*1
  qe6_tmp_.dim = SET_B(NUM_RF, 1);
  qe6_tmp_.sps = qe6_spscnt_;
  qe6_spscnt_ += OCU(qe6_tmp_);
  // Initialize q:
  qe6_qw_ = 1;
  qe6_qx_ = 0;
  qe6_qy_ = 0;
  qe6_qz_ = 0;
}

void set_qe6_J_QE6(void) {
  // qe6_qz_ must be 0
  qe6_qz_ = 0;
  // Set J by reference vectors and quaternions
  MAT(qe6_J_, 0, 0) =
      -2 * (qe6_qw_ * qe6_bax_ - qe6_qz_ * qe6_bay_ + qe6_qy_ * qe6_baz_);
  MAT(qe6_J_, 1, 0) =
      -2 * (qe6_qz_ * qe6_bax_ + qe6_qw_ * qe6_bay_ - qe6_qx_ * qe6_baz_);
  MAT(qe6_J_, 2, 0) =
      -2 * (-qe6_qy_ * qe6_bax_ + qe6_qx_ * qe6_bay_ + qe6_qw_ * qe6_baz_);

  MAT(qe6_J_, 0, 1) =
      -2 * (qe6_qx_ * qe6_bax_ + qe6_qy_ * qe6_bay_ + qe6_qz_ * qe6_baz_);
  MAT(qe6_J_, 1, 1) =
      -2 * (qe6_qy_ * qe6_bax_ - qe6_qx_ * qe6_bay_ - qe6_qw_ * qe6_baz_);
  MAT(qe6_J_, 2, 1) =
      -2 * (qe6_qz_ * qe6_bax_ + qe6_qw_ * qe6_bay_ - qe6_qx_ * qe6_baz_);

  MAT(qe6_J_, 0, 2) =
      -2 * (-qe6_qy_ * qe6_bax_ + qe6_qx_ * qe6_bay_ + qe6_qw_ * qe6_baz_);
  MAT(qe6_J_, 1, 2) =
      -2 * (qe6_qx_ * qe6_bax_ + qe6_qy_ * qe6_bay_ + qe6_qz_ * qe6_baz_);
  MAT(qe6_J_, 2, 2) =
      -2 * (-qe6_qw_ * qe6_bax_ + qe6_qz_ * qe6_bay_ - qe6_qy_ * qe6_baz_);
  // PD = J^T*J
  for (int i = 0; i < NUM_QT; i++) {
    for (int j = 0; j < NUM_QT; j++) {
      SYM(qe6_PD_, i, j) = 0;
      for (int k = 0; k < NUM_RF; k++) {
        SYM(qe6_PD_, i, j) += MAT(qe6_J_, k, i) * MAT(qe6_J_, k, j);
      }
    }
  }
}

void qe6_PD_invert_QE6(void) {
  SYM(qe6_PD_, 0, 0) = 1.0f / SYM(qe6_PD_, 0, 0);
  for (int i = 0; i < NUM_QT - 1; i++) {
    // w^T = -inv(A)*u^T
    // Here, the submatrix up to I(i, i) (that is, A) was already inverted.
    for (int j = 0; j <= i; j++) {
      MAT(qe6_tmp_, j, 0) = 0;
      for (int k = 0; k <= i; k++) {
        MAT(qe6_tmp_, j, 0) -= SYM(qe6_PD_, j, k) * SYM(qe6_PD_, k, i + 1);
      }
    }
    // y = inv(u*w^T+x)
    // Firstly u*w^T:
    FTYPE y = 0;
    for (int j = 0; j <= i; j++) {
      y += SYM(qe6_PD_, j, i + 1) * MAT(qe6_tmp_, j, 0);
    }
    SYM(qe6_PD_, i + 1, i + 1) = 1.0f / (y + SYM(qe6_PD_, i + 1, i + 1));
    // v^T = w^T * y
    for (int j = 0; j <= i; j++) {
      SYM(qe6_PD_, j, i + 1) = MAT(qe6_tmp_, j, 0) * SYM(qe6_PD_, i + 1, i + 1);
    }
    // B = inv(A) + w^T*v
    for (int j = 0; j <= i; j++) {
      for (int k = j; k <= i; k++) {
        SYM(qe6_PD_, j, k) =
            SYM(qe6_PD_, j, k) + MAT(qe6_tmp_, j, 0) * SYM(qe6_PD_, k, i + 1);
      }
    }
  }
  return;
}

void get_error_QE6(void) {
  FTYPE R11 =
      pow(qe6_qw_, 2) + pow(qe6_qx_, 2) - pow(qe6_qy_, 2) - pow(qe6_qz_, 2);
  FTYPE R12 = 2 * (qe6_qx_ * qe6_qy_ - qe6_qw_ * qe6_qz_);
  FTYPE R13 = 2 * (qe6_qx_ * qe6_qz_ + qe6_qw_ * qe6_qy_);
  FTYPE R21 = 2 * (qe6_qx_ * qe6_qy_ + qe6_qw_ * qe6_qz_);
  FTYPE R22 =
      pow(qe6_qw_, 2) + pow(qe6_qy_, 2) - pow(qe6_qx_, 2) - pow(qe6_qz_, 2);
  FTYPE R23 = 2 * (qe6_qy_ * qe6_qz_ - qe6_qw_ * qe6_qx_);
  FTYPE R31 = 2 * (qe6_qx_ * qe6_qz_ - qe6_qw_ * qe6_qy_);
  FTYPE R32 = 2 * (qe6_qy_ * qe6_qz_ + qe6_qw_ * qe6_qx_);
  FTYPE R33 =
      pow(qe6_qw_, 2) + pow(qe6_qz_, 2) - pow(qe6_qy_, 2) - pow(qe6_qx_, 2);

  MAT(qe6_e_, 0, 0) =
      qe6_eax_ - R11 * qe6_bax_ - R12 * qe6_bay_ - R13 * qe6_baz_;
  MAT(qe6_e_, 1, 0) =
      qe6_eay_ - R21 * qe6_bax_ - R22 * qe6_bay_ - R23 * qe6_baz_;
  MAT(qe6_e_, 2, 0) =
      qe6_eaz_ - R31 * qe6_bax_ - R32 * qe6_bay_ - R33 * qe6_baz_;
}

void quaternion_update_QE6(void) {
  set_qe6_J_QE6();
  qe6_PD_invert_QE6();
  get_error_QE6();
  // get dq = PD*J^T*e
  for (int i = 0; i < NUM_QT; i++) {
    // dq_i
    MAT(qe6_tmp_, i, 0) = 0;
    for (int j = 0; j < NUM_QT; j++) {
      for (int k = 0; k < NUM_RF; k++) {
        MAT(qe6_tmp_, i, 0) +=
            SYM(qe6_PD_, i, j) * MAT(qe6_J_, k, j) * MAT(qe6_e_, k, 0);
      }
    }
  }
  // q_k+1 = q_k + dq
  qe6_qw_ -= MAT(qe6_tmp_, 0, 0);
  qe6_qx_ -= MAT(qe6_tmp_, 1, 0);
  qe6_qy_ -= MAT(qe6_tmp_, 2, 0);
  qe6_qz_ = 0;
}

void get_qk_QE6(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  *qw = qe6_qw_;
  *qx = qe6_qx_;
  *qy = qe6_qy_;
  *qz = qe6_qz_;
}

void get_rpy_QE6(FTYPE *roll, FTYPE *pitch, FTYPE *yaw) {
  *roll = atan2(2 * (qe6_qw_ * qe6_qx_ + qe6_qy_ * qe6_qz_),
               1 - 2 * (qe6_qx_ * qe6_qx_ + qe6_qy_ * qe6_qy_));
  *pitch = asin(2 * (qe6_qw_ * qe6_qy_ - qe6_qz_ * qe6_qx_));
  *yaw = atan2(2 * (qe6_qw_ * qe6_qz_ + qe6_qx_ * qe6_qy_),
                1 - 2 * (qe6_qy_ * qe6_qy_ + qe6_qz_ * qe6_qz_));
}
