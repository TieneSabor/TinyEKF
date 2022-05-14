#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include "quaternion_9axis.h"

// Important types used here.
typedef struct matrix {
  byte dim; // dimision of the matrix, dim = m<<4 + n
  byte sps; // start position in thqe9_e_space
} M;

// Definition of important constant.
// Should NOT be modified.
#define SPACE_SIZE 55 // maximum matrix element number
#define MAX_DIM 15    // maximum matrix operation
#define LAST4 0x0f    // 00001111
#define NUM_RF 6 // reference vector elements number (2 3D vector make it 6)
#define NUM_QT 4 // quaternion element number

// Byte operation. For "dim" in matrix struct, it is composed of M(first 4 bit)
// and N(last 4 bit)
#define GET_M(x) (x >> 4)
#define GET_N(x) (x & 0x0f)
#define SET_B(m, n) ((m << 4) + n)

// The (x, y)th element for NOT symmetric matrix of m rows * n columns, starting
// with (0, 0). With O2 Opt., the add/mul. part in the index should be replaced
// by const. (Hopefully?)
#define MAT(mat, x, y) qe9_spacqe9_e_[mat.sps + y + x * (mat.dim & 0x0f)]
// Get the (x, y)th element of a symmetric matrix.
#define SYM(mat, x, y)                                                         \
  qe9_spacqe9_e_[mat.sps + imax(x, y) + imin(x, y) * (mat.dim & 0x0f) -        \
                 (imin(x, y) * imin(x, y) + imin(x, y)) / 2]
// Get how much a NOT symmetric matrix take to store data
#define OCU(mat) (mat.dim >> 4) * (mat.dim & 0x0f)
// Get how much a symmetric matrix take to store data
#define OCS(mat)                                                               \
  0.5 * ((mat.dim >> 4) * ((mat.dim & 0x0f) - 1)) + (mat.dim >> 4)

// Every element of matrix in the quaternion estimator will be stored here.
FTYPE qe9_spacqe9_e_[SPACE_SIZE];
byte qe9_spscnt_;

// Matrix used for Quaternion estimation
M qe9_J_;   // Jacobian, 6*4
M qe9_PD_;  // Positive Definite Matrix of J^T*J, 4*4
M qe9_e_;   // Error between earth frame references and measured directions, 6*1
M qe9_tmp_; // Temp space for matrix multiplication, 6*1

// Reference vector value (earth frame) and measured vector value (body frame)
FTYPE qe9_eax_, qe9_eay_, qe9_eaz_, qe9_emx_, qe9_emy_, qe9_emz_;
FTYPE qe9_bax_, qe9_bay_, qe9_baz_, qe9_bmx_, qe9_bmy_, qe9_bmz_;
// Quaternion Estimation
FTYPE qe9_qw_, qe9_qx_, qe9_qy_, qe9_qz_;

void set_reference(FTYPE eax, FTYPE eay, FTYPE eaz, FTYPE emx, FTYPE emy,
                   FTYPE emz) {
  qe9_eax_ = eax;
  qe9_eay_ = eay;
  qe9_eaz_ = eaz;
  qe9_emx_ = emx;
  qe9_emy_ = emy;
  qe9_emz_ = emz;
}

void set_measurement(FTYPE bax, FTYPE bay, FTYPE baz, FTYPE bmx, FTYPE bmy,
                     FTYPE bmz) {
  qe9_bax_ = bax;
  qe9_bay_ = bay;
  qe9_baz_ = baz;
  qe9_bmx_ = bmx;
  qe9_bmy_ = bmy;
  qe9_bmz_ = bmz;
}

void quaternion_init(void) {
  // Fill all elements with 0
  for (int i = 0; i < SPACE_SIZE; i++) {
    qe9_spacqe9_e_[i] = 0;
  }
  // Initialize all structures.
  // Jacobian, non sysmmetric 6*4
  qe9_spscnt_ = 0;
  qe9_J_.dim = SET_B(NUM_RF, NUM_QT);
  qe9_J_.sps = qe9_spscnt_;
  qe9_spscnt_ += OCU(qe9_J_);
  // Positive Definite, symmetric 4*4
  qe9_PD_.dim = SET_B(NUM_QT, NUM_QT);
  qe9_PD_.sps = qe9_spscnt_;
  qe9_spscnt_ += OCS(qe9_PD_);
  // Error, vector 6*1
  qe9_e_.dim = SET_B(NUM_RF, 1);
  qe9_e_.sps = qe9_spscnt_;
  qe9_spscnt_ += OCU(qe9_e_);
  // Tmp space, vector 6*1
  qe9_tmp_.dim = SET_B(NUM_RF, 1);
  qe9_tmp_.sps = qe9_spscnt_;
  qe9_spscnt_ += OCU(qe9_tmp_);
  // Initialize q:
  qe9_qw_ = 1;
  qe9_qx_ = 0;
  qe9_qy_ = 0;
  qe9_qz_ = 0;
}

void set_J(void) {
  // Set J by reference vectors and quaternions
  MAT(qe9_J_, 0, 0) =
      -2 * (qe9_qw_ * qe9_bmx_ - qe9_qz_ * qe9_bmy_ + qe9_qy_ * qe9_bmz_);
  MAT(qe9_J_, 1, 0) =
      -2 * (qe9_qz_ * qe9_bmx_ + qe9_qw_ * qe9_bmy_ - qe9_qx_ * qe9_bmz_);
  MAT(qe9_J_, 2, 0) =
      -2 * (-qe9_qy_ * qe9_bmx_ + qe9_qx_ * qe9_bmy_ + qe9_qw_ * qe9_bmz_);
  MAT(qe9_J_, 3, 0) =
      -2 * (qe9_qw_ * qe9_bax_ - qe9_qz_ * qe9_bay_ + qe9_qy_ * qe9_baz_);
  MAT(qe9_J_, 4, 0) =
      -2 * (qe9_qz_ * qe9_bax_ + qe9_qw_ * qe9_bay_ - qe9_qx_ * qe9_baz_);
  MAT(qe9_J_, 5, 0) =
      -2 * (-qe9_qy_ * qe9_bax_ + qe9_qx_ * qe9_bay_ + qe9_qw_ * qe9_baz_);

  MAT(qe9_J_, 0, 1) =
      -2 * (qe9_qx_ * qe9_bmx_ + qe9_qy_ * qe9_bmy_ + qe9_qz_ * qe9_bmz_);
  MAT(qe9_J_, 1, 1) =
      -2 * (qe9_qy_ * qe9_bmx_ - qe9_qx_ * qe9_bmy_ - qe9_qw_ * qe9_bmz_);
  MAT(qe9_J_, 2, 1) =
      -2 * (qe9_qz_ * qe9_bmx_ + qe9_qw_ * qe9_bmy_ - qe9_qx_ * qe9_bmz_);
  MAT(qe9_J_, 3, 1) =
      -2 * (qe9_qx_ * qe9_bax_ + qe9_qy_ * qe9_bay_ + qe9_qz_ * qe9_baz_);
  MAT(qe9_J_, 4, 1) =
      -2 * (qe9_qy_ * qe9_bax_ - qe9_qx_ * qe9_bay_ - qe9_qw_ * qe9_baz_);
  MAT(qe9_J_, 5, 1) =
      -2 * (qe9_qz_ * qe9_bax_ + qe9_qw_ * qe9_bay_ - qe9_qx_ * qe9_baz_);

  MAT(qe9_J_, 0, 2) =
      -2 * (-qe9_qy_ * qe9_bmx_ + qe9_qx_ * qe9_bmy_ + qe9_qw_ * qe9_bmz_);
  MAT(qe9_J_, 1, 2) =
      -2 * (qe9_qx_ * qe9_bmx_ + qe9_qy_ * qe9_bmy_ + qe9_qz_ * qe9_bmz_);
  MAT(qe9_J_, 2, 2) =
      -2 * (-qe9_qw_ * qe9_bmx_ + qe9_qz_ * qe9_bmy_ - qe9_qy_ * qe9_bmz_);
  MAT(qe9_J_, 3, 2) =
      -2 * (-qe9_qy_ * qe9_bax_ + qe9_qx_ * qe9_bay_ + qe9_qw_ * qe9_baz_);
  MAT(qe9_J_, 4, 2) =
      -2 * (qe9_qx_ * qe9_bax_ + qe9_qy_ * qe9_bay_ + qe9_qz_ * qe9_baz_);
  MAT(qe9_J_, 5, 2) =
      -2 * (-qe9_qw_ * qe9_bax_ + qe9_qz_ * qe9_bay_ - qe9_qy_ * qe9_baz_);

  MAT(qe9_J_, 0, 3) =
      -2 * (-qe9_qz_ * qe9_bmx_ - qe9_qw_ * qe9_bmy_ + qe9_qx_ * qe9_bmz_);
  MAT(qe9_J_, 1, 3) =
      -2 * (qe9_qw_ * qe9_bmx_ - qe9_qz_ * qe9_bmy_ + qe9_qy_ * qe9_bmz_);
  MAT(qe9_J_, 2, 3) =
      -2 * (qe9_qx_ * qe9_bmx_ + qe9_qy_ * qe9_bmy_ + qe9_qz_ * qe9_bmz_);
  MAT(qe9_J_, 3, 3) =
      -2 * (-qe9_qz_ * qe9_bax_ - qe9_qw_ * qe9_bay_ + qe9_qx_ * qe9_baz_);
  MAT(qe9_J_, 4, 3) =
      -2 * (qe9_qw_ * qe9_bax_ - qe9_qz_ * qe9_bay_ + qe9_qy_ * qe9_baz_);
  MAT(qe9_J_, 5, 3) =
      -2 * (qe9_qx_ * qe9_bax_ + qe9_qy_ * qe9_bay_ + qe9_qz_ * qe9_baz_);
  // PD = J^T*J
  for (int i = 0; i < NUM_QT; i++) {
    for (int j = 0; j < NUM_QT; j++) {
      SYM(qe9_PD_, i, j) = 0;
      for (int k = 0; k < NUM_RF; k++) {
        SYM(qe9_PD_, i, j) += MAT(qe9_J_, k, i) * MAT(qe9_J_, k, j);
      }
    }
  }
}

void qe9_PD_invert(void) {
  SYM(qe9_PD_, 0, 0) = 1.0f / SYM(qe9_PD_, 0, 0);
  for (int i = 0; i < NUM_QT - 1; i++) {
    // w^T = -inv(A)*u^T
    // Here, the submatrix up to I(i, i) (that is, A) was already inverted.
    for (int j = 0; j <= i; j++) {
      MAT(qe9_tmp_, j, 0) = 0;
      for (int k = 0; k <= i; k++) {
        MAT(qe9_tmp_, j, 0) -= SYM(qe9_PD_, j, k) * SYM(qe9_PD_, k, i + 1);
      }
    }
    // y = inv(u*w^T+x)
    // Firstly u*w^T:
    FTYPE y = 0;
    for (int j = 0; j <= i; j++) {
      y += SYM(qe9_PD_, j, i + 1) * MAT(qe9_tmp_, j, 0);
    }
    SYM(qe9_PD_, i + 1, i + 1) = 1.0f / (y + SYM(qe9_PD_, i + 1, i + 1));
    // v^T = w^T * y
    for (int j = 0; j <= i; j++) {
      SYM(qe9_PD_, j, i + 1) = MAT(qe9_tmp_, j, 0) * SYM(qe9_PD_, i + 1, i + 1);
    }
    // B = inv(A) + w^T*v
    for (int j = 0; j <= i; j++) {
      for (int k = j; k <= i; k++) {
        SYM(qe9_PD_, j, k) =
            SYM(qe9_PD_, j, k) + MAT(qe9_tmp_, j, 0) * SYM(qe9_PD_, k, i + 1);
      }
    }
  }
  return;
}

void get_error(void) {
  FTYPE R11 =
      pow(qe9_qw_, 2) + pow(qe9_qx_, 2) - pow(qe9_qy_, 2) - pow(qe9_qz_, 2);
  FTYPE R12 = 2 * (qe9_qx_ * qe9_qy_ - qe9_qw_ * qe9_qz_);
  FTYPE R13 = 2 * (qe9_qx_ * qe9_qz_ + qe9_qw_ * qe9_qy_);
  FTYPE R21 = 2 * (qe9_qx_ * qe9_qy_ + qe9_qw_ * qe9_qz_);
  FTYPE R22 =
      pow(qe9_qw_, 2) + pow(qe9_qy_, 2) - pow(qe9_qx_, 2) - pow(qe9_qz_, 2);
  FTYPE R23 = 2 * (qe9_qy_ * qe9_qz_ - qe9_qw_ * qe9_qx_);
  FTYPE R31 = 2 * (qe9_qx_ * qe9_qz_ - qe9_qw_ * qe9_qy_);
  FTYPE R32 = 2 * (qe9_qy_ * qe9_qz_ + qe9_qw_ * qe9_qx_);
  FTYPE R33 =
      pow(qe9_qw_, 2) + pow(qe9_qz_, 2) - pow(qe9_qy_, 2) - pow(qe9_qx_, 2);

  MAT(qe9_e_, 0, 0) =
      qe9_emx_ - R11 * qe9_bmx_ - R12 * qe9_bmy_ - R13 * qe9_bmz_;
  MAT(qe9_e_, 1, 0) =
      qe9_emy_ - R21 * qe9_bmx_ - R22 * qe9_bmy_ - R23 * qe9_bmz_;
  MAT(qe9_e_, 2, 0) =
      qe9_emz_ - R31 * qe9_bmx_ - R32 * qe9_bmy_ - R33 * qe9_bmz_;

  MAT(qe9_e_, 3, 0) =
      qe9_eax_ - R11 * qe9_bax_ - R12 * qe9_bay_ - R13 * qe9_baz_;
  MAT(qe9_e_, 4, 0) =
      qe9_eay_ - R21 * qe9_bax_ - R22 * qe9_bay_ - R23 * qe9_baz_;
  MAT(qe9_e_, 5, 0) =
      qe9_eaz_ - R31 * qe9_bax_ - R32 * qe9_bay_ - R33 * qe9_baz_;
}

void quaternion_update(void) {
  set_J();
  qe9_PD_invert();
  get_error();
  // get dq = PD*J^T*e
  for (int i = 0; i < NUM_QT; i++) {
    // dq_i
    MAT(qe9_tmp_, i, 0) = 0;
    for (int j = 0; j < NUM_QT; j++) {
      for (int k = 0; k < NUM_RF; k++) {
        MAT(qe9_tmp_, i, 0) +=
            SYM(qe9_PD_, i, j) * MAT(qe9_J_, k, j) * MAT(qe9_e_, k, 0);
      }
    }
  }
  // q_k+1 = q_k + dq
  qe9_qw_ -= MAT(qe9_tmp_, 0, 0);
  qe9_qx_ -= MAT(qe9_tmp_, 1, 0);
  qe9_qy_ -= MAT(qe9_tmp_, 2, 0);
  qe9_qz_ -= MAT(qe9_tmp_, 3, 0);
}

void get_qe9_qk(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  *qw = qe9_qw_;
  *qx = qe9_qx_;
  *qy = qe9_qy_;
  *qz = qe9_qz_;
}

void get_qe9_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw) {
  *roll = atan2(2 * (qe9_qw_ * qe9_qx_ + qe9_qy_ * qe9_qz_),
                1 - 2 * (qe9_qx_ * qe9_qx_ + qe9_qy_ * qe9_qy_));
  *pitch = asin(2 * (qe9_qw_ * qe9_qy_ - qe9_qz_ * qe9_qx_));
  *yaw = atan2(2 * (qe9_qw_ * qe9_qz_ + qe9_qx_ * qe9_qy_),
               1 - 2 * (qe9_qy_ * qe9_qy_ + qe9_qz_ * qe9_qz_));
}
