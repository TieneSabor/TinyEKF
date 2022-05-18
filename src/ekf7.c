#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include "ekf7.h"

#ifdef __cplusplus
extern "C" {
#endif

// #define DBG_CODE

#ifdef DBG_CODE
#define PRT_DBG 1
#define EKF_DBG(format, ...)                                                   \
  if (PRT_DBG) {                                                               \
    fprintf(stderr, "(%s:%d) " format "\n", __FILE__, __LINE__,                \
            ##__VA_ARGS__);                                                    \
  }
#endif

// Definition of important constant.
// Should NOT be modified.
#define SPACE_SIZE 130 // maximum matrix element number
#define MAX_DIM 15     // maximum matrix operation
#define LAST4 0x0f     // 00001111
#define NUM_ST 7       // state number
#define NUM_OB 4       // observation number

// Byte operation. For "dim" in matrix struct, it is composed of M(first 4 bit)
// and N(last 4 bit)
#define GET_M(x) (x >> 4)
#define GET_N(x) (x & 0x0f)
#define SET_B(m, n) ((m << 4) + n)

// The (x, y)th element for NOT symmetric matrix of m rows * n columns, starting
// with (0, 0). With O2 Opt., the add/mul. part in the index should be replaced
// by const. (Hopefully?)
#define MAT(mat, x, y) ekf7_space_[mat.sps + y + x * (mat.dim & 0x0f)]
// Get the (x, y)th element of a symmetric matrix.
#define SYM(mat, x, y)                                                         \
  ekf7_space_[mat.sps + imax(x, y) + imin(x, y) * (mat.dim & 0x0f) -           \
              (imin(x, y) * imin(x, y) + imin(x, y)) / 2]
// Get how much a NOT symmetric matrix take to store data
#define OCU(mat) (mat.dim >> 4) * (mat.dim & 0x0f)
// Get how much a symmetric matrix take to store data
#define OCS(mat)                                                               \
  0.5 * ((mat.dim >> 4) * ((mat.dim & 0x0f) - 1)) + (mat.dim >> 4)

// Some special variable stored in space
#define BX ekf7_space_[ekf7_state_.sps + 0]
#define BY ekf7_space_[ekf7_state_.sps + 1]
#define BZ ekf7_space_[ekf7_state_.sps + 2]
#define QW ekf7_space_[ekf7_state_.sps + 3]
#define QX ekf7_space_[ekf7_state_.sps + 4]
#define QY ekf7_space_[ekf7_state_.sps + 5]
#define QZ ekf7_space_[ekf7_state_.sps + 6]

// Every element of matrix in the EKF will be stored here.
FTYPE ekf7_space_[SPACE_SIZE];

// Matrix used for EKF estimation
M ekf7_state_;
M ekf7_error_;
M ekf7_covar_;
M ekf7_W_;
M ekf7_V_;
M ekf7_I_;
M ekf7_KH_;
M ekf7_tmp_;

// Some flags for checking initialization
byte ekf7_init_called_ = 0;
byte ekf7_PN_set_ = 0;
byte ekf7_ON_set_ = 0;

// Important parameters
#define JACSIZ 10 // size used for jacobian matrix
#define FQW ekf7_space_[0]
#define FQX ekf7_space_[1]
#define FQY ekf7_space_[2]
#define FQZ ekf7_space_[3]
#define FBX ekf7_space_[4]
#define FBY ekf7_space_[5]
#define FBZ ekf7_space_[6]
#define FWX ekf7_space_[7]
#define FWY ekf7_space_[8]
#define FWZ ekf7_space_[9]

byte ekf7_spscnt_;

void ekf7_set_state(FTYPE bx, FTYPE by, FTYPE bz, FTYPE qw, FTYPE qx, FTYPE qy,
                    FTYPE qz) {
  BX = bx;
  BY = by;
  BZ = bz;
  QW = qw;
  QX = qx;
  QY = qy;
  QZ = qz;
}

// Set linearlize transition matrix F.
void ekf7_set_F(FTYPE WX, FTYPE WY, FTYPE WZ) {
  if (!ekf7_init_called_) {
    return;
  }
  // Set up quaternion in F
  FQW = QW;
  FQX = QX;
  FQY = QY;
  FQZ = QZ;
  // Set up gyroscope in F
  FBX = BX;
  FBY = BY;
  FBZ = BZ;
  // Set up angular velocity in F
  FWX = WX;
  FWY = WY;
  FWZ = WZ;
}

FTYPE ekf7_get_F(FTYPE dt, byte i, byte j) {
  FTYPE res = 0;
  // All element that i < 3 is 0
  if (i > 3) {
    // -qw
    if ((i == (j + 4)) && (i > 3)) {
      res = -(FQW * 0.5);
    }
    // -qx
    else if (((i == 6) && (j == 1))) {
      res = -(FQX * 0.5);
    }
    // qx
    else if (((j == 3) && (i == 0)) || ((i == 5) && (j == 2))) {
      res = (FQX * 0.5);
    }
    // -qy
    else if ((i == 4) && (j == 2)) {
      res = -(FQY * 0.5);
    }
    // qy
    else if (((i == 3) && (j == 1)) || ((i == 6) && (j == 0))) {
      res = (FQY * 0.5);
    }
    // -qz
    else if ((i == 5) && (j == 0)) {
      res = -(FQZ * 0.5);
    }
    // qz
    else if (((i == 3) && (j == 2)) || ((i == 4) && (j == 1))) {
      res = (FQZ * 0.5);
    }
    // wx-bx
    else if (((i == 4) && (j == 3)) || ((i == 5) && (j == 6))) {
      res = ((FWX - FBX) * 0.5);
    }
    // bx-wx
    else if (((i == 3) && (j == 4)) || ((i == 6) && (j == 5))) {
      res = ((FBX - FWX) * 0.5);
    }
    // wy-by
    else if (((i == 6) && (j == 4)) || ((i == 5) && (j == 3))) {
      res = ((FWY - FBY) * 0.5);
    }
    // by-wy
    else if (((i == 4) && (j == 6)) || ((i == 3) && (j == 5))) {
      res = ((FBY - FWY) * 0.5);
    }
    // wz-bz
    else if (((i == 4) && (j == 5)) || ((i == 6) && (j == 3))) {
      res = ((FWZ - FBZ) * 0.5);
    }
    // bz-wz
    else if (((i == 5) && (j == 4)) || ((i == 3) && (j == 6))) {
      res = ((FBZ - FWZ) * 0.5);
    }
    // others are 0
    else {
      res = 0;
    }
  }
  /*
  // i == 3
  if (i == 3) {
    switch (j) {
    case 0:
      res = (FQX * 0.5);
      break;
    case 1:
      res = (FQY * 0.5);
      break;
    case 2:
      res = (FQZ * 0.5);
      break;
    case 4:
      res = ((FBX - FWX) * 0.5);
      break;
    case 5:
      res = ((FBY - FWY) * 0.5);
      break;
    case 6:
      res = ((FBZ - FWZ) * 0.5);
      break;
    default:
      res = 0;
      break;
    }
  }
  // i == 4
  else if (i == 4) {
    switch (j) {
    case 1:
      res = -(FQZ * 0.5);
      break;
    case 2:
      res = (FQY * 0.5);
      break;
    case 3:
      res = ((FWX - FBX) * 0.5);
      break;
    case 5:
      res = ((FBZ - FWZ) * 0.5);
      break;
    case 6:
      res = ((FWY - FBY) * 0.5);
      break;
    default:
      res = 0;
      break;
    }
  }
  // i == 5
  else if (i == 5) {
    switch (j) {
    case 0:
      res = (FQZ * 0.5);
      break;
    case 2:
      res = -(FQX * 0.5);
      break;
    case 3:
      res = ((FWY - FBY) * 0.5);
      break;
    case 4:
      res = ((FWZ - FBZ) * 0.5);
      break;
    case 6:
      res = ((FBX - FWX) * 0.5);
      break;
    default:
      res = 0;
      break;
    }
  }
  // i == 6
  else if (i == 6) {
    switch (j) {
    case 0:
      res = -(FQY * 0.5);
      break;
    case 1:
      res = (FQX * 0.5);
      break;
    case 3:
      res = ((FWZ - FBZ) * 0.5);
      break;
    case 4:
      res = ((FBY - FWY) * 0.5);
      break;
    case 5:
      res = ((FWX - FBX) * 0.5);
      break;
    default:
      res = 0;
      break;
    }
  } else {
    res = 0;
  }
  */

  // return I + Fdt
  res *= dt;
  return res;
}

FTYPE ekf7_get_F_plus_I(FTYPE dt, byte i, byte j) {
  if (i == j) {
    return 1.0f + ekf7_get_F(dt, i, j);
  } else {
    return ekf7_get_F(dt, i, j);
  }
}

FTYPE ekf7_get_H(byte i, byte j) {
  if ((j >= 3) && ((i + 3) == j)) {
    return 1.0f;
  } else {
    return 0.0f;
  }
}

void ekf7_init(void) {
  ekf7_init_called_ = 1;
  // Fill all elements with 0
  for (int i = 0; i < SPACE_SIZE; i++) {
    ekf7_space_[i] = 0;
  }
  // Initialize all structures.
  // State, a vector
  // int spscnt = 0 + JACSIZ;
  ekf7_spscnt_ = 0 + JACSIZ;
  ekf7_state_.dim = SET_B(NUM_ST, 1);
  ekf7_state_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_state_);
  // Covariance, a symmetric matrix
  ekf7_covar_.dim = SET_B(NUM_ST, NUM_ST);
  ekf7_covar_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCS(ekf7_covar_);
  // Propogation niose, a diagnal matrix
  ekf7_W_.dim = SET_B(NUM_ST, 1);
  ekf7_W_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_W_);
  // Observation noise, a diagnal matrix
  ekf7_V_.dim = SET_B(NUM_OB, 1);
  ekf7_V_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_V_);
  // Innovation matrix, a symmetric matrix
  ekf7_I_.dim = SET_B(NUM_OB, NUM_OB);
  ekf7_I_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCS(ekf7_I_);
  // KH = PH^T(I^-1)H, not a symmetric matrix
  ekf7_KH_.dim = SET_B(NUM_ST, NUM_ST);
  ekf7_KH_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_KH_);
  // error = (H^+)y - x
  ekf7_error_.dim = SET_B(NUM_ST, 1);
  ekf7_error_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_error_);
  // for matrix multiplation
  ekf7_tmp_.dim = SET_B(NUM_ST, 1);
  ekf7_tmp_.sps = ekf7_spscnt_;
  ekf7_spscnt_ += OCU(ekf7_tmp_);
  printf("spscnt: %d\r\n", ekf7_spscnt_);
  // Initialize state
  MAT(ekf7_state_, 3, 0) = 1; // qw = 1;
  // Initialize covariance P.
  // for (int i = 0; i < NUM_ST; i++) {
  //  SYM(ekf7_covar_, i, i) = 1.0;
  //}
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = 0; j < NUM_ST; j++) {
      SYM(ekf7_covar_, i, j) = 1e-2;
    }
  }
  // Initialize Jacobian F.
  ekf7_set_F(0, 0, 0);
}

void ekf7_set_Propagation_Noise(FTYPE Wbx, FTYPE Wby, FTYPE Wbz, FTYPE Wq1,
                                FTYPE Wq2, FTYPE Wq3, FTYPE Wq4) {
  if (!ekf7_init_called_) {
    return;
  }
  // Fill the value
  MAT(ekf7_W_, 0, 0) = Wbx;
  MAT(ekf7_W_, 1, 0) = Wby;
  MAT(ekf7_W_, 2, 0) = Wbz;
  MAT(ekf7_W_, 3, 0) = Wq1;
  MAT(ekf7_W_, 4, 0) = Wq2;
  MAT(ekf7_W_, 5, 0) = Wq3;
  MAT(ekf7_W_, 6, 0) = Wq4;
  // Set flag
  ekf7_PN_set_ = 1;
}

void ekf7_set_Observation_Noise(FTYPE Vq1, FTYPE Vq2, FTYPE Vq3, FTYPE Vq4) {
  if (!ekf7_init_called_) {
    return;
  }
  // Fill the value
  MAT(ekf7_V_, 0, 0) = Vq1;
  MAT(ekf7_V_, 1, 0) = Vq2;
  MAT(ekf7_V_, 2, 0) = Vq3;
  MAT(ekf7_V_, 3, 0) = Vq4;
  // Set flag
  ekf7_ON_set_ = 1;
}

void ekf7_predict(FTYPE dt, FTYPE wx, FTYPE wy, FTYPE wz,
                  char covar_update_flag) {
  // Update the Jacobian.
  ekf7_set_F(wx, wy, wz);
  // printF();
  // Update state.
  FTYPE nqw =
      QW + 0.5 * dt * (-(wx - BX) * QX - (wy - BY) * QY - (wz - BZ) * QZ);
  FTYPE nqx =
      QX + 0.5 * dt * ((wx - BX) * QW + (wz - BZ) * QY - (wy - BY) * QZ);
  FTYPE nqy =
      QY + 0.5 * dt * ((wy - BY) * QW - (wz - BZ) * QX + (wx - BX) * QZ);
  FTYPE nqz =
      QZ + 0.5 * dt * ((wz - BZ) * QW + (wy - BY) * QX - (wx - BX) * QY);
  FTYPE norm = sqrt(pow(nqw, 2) + pow(nqx, 2) + pow(nqy, 2) + pow(nqz, 2));
  // FTYPE norm = 1;
  QW = nqw / norm;
  QX = nqx / norm;
  QY = nqy / norm;
  QZ = nqz / norm;
  // Update the covariance
  // P = FPF^T + W
  if (covar_update_flag == EKF7_COV_UPDATE) {
    ekf7_covariance_predict(dt);
  }
}

void ekf7_covariance_predict(FTYPE dt) {
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = i; j < NUM_ST; j++) {
      MAT(ekf7_tmp_, j, 0) = 0;
      if (i == j) {
        MAT(ekf7_tmp_, j, 0) = MAT(ekf7_W_, i, 0);
      }
      for (int k = 0; k < NUM_ST; k++) {
        for (int l = 0; l < NUM_ST; l++) {
          MAT(ekf7_tmp_, j, 0) += ekf7_get_F_plus_I(dt, i, k) *
                                  SYM(ekf7_covar_, k, l) *
                                  ekf7_get_F_plus_I(dt, j, l);
        }
      }
    }
    for (int j = i; j < NUM_ST; j++) {
      SYM(ekf7_covar_, i, j) = MAT(ekf7_tmp_, j, 0);
    }
  }
}

void ekf7_innovation_update(void) {
  for (int i = 0; i < NUM_OB; i++) {
    for (int j = i; j < NUM_OB; j++) {
      MAT(ekf7_tmp_, j, 0) = 0;
      if (i == j) {
        MAT(ekf7_tmp_, j, 0) = MAT(ekf7_V_, i, 0);
      }
      for (int k = 0; k < NUM_ST; k++) {
        for (int l = 0; l < NUM_ST; l++) {
          MAT(ekf7_tmp_, j, 0) +=
              ekf7_get_H(i, k) * SYM(ekf7_covar_, k, l) * ekf7_get_H(j, l);
        }
      }
    }
    for (int j = i; j < NUM_ST; j++) {
      SYM(ekf7_I_, i, j) = MAT(ekf7_tmp_, j, 0);
    }
  }
}

// Matrix Inversion Routine
void ekf7_innovation_invert(void) {
  SYM(ekf7_I_, 0, 0) = 1.0f / SYM(ekf7_I_, 0, 0);
  for (int i = 0; i < NUM_OB - 1; i++) {
    // w^T = -inv(A)*u^T
    // Here, the submatrix up to I(i, i) (that is, A) was already inverted.
    for (int j = 0; j <= i; j++) {
      MAT(ekf7_tmp_, j, 0) = 0;
      for (int k = 0; k <= i; k++) {
        MAT(ekf7_tmp_, j, 0) -= SYM(ekf7_I_, j, k) * SYM(ekf7_I_, k, i + 1);
      }
    }
    // y = inv(u*w^T+x)
    // Firstly u*w^T:
    FTYPE y = 0;
    for (int j = 0; j <= i; j++) {
      y += SYM(ekf7_I_, j, i + 1) * MAT(ekf7_tmp_, j, 0);
    }
    SYM(ekf7_I_, i + 1, i + 1) = 1.0f / (y + SYM(ekf7_I_, i + 1, i + 1));
    // v^T = w^T * y
    for (int j = 0; j <= i; j++) {
      SYM(ekf7_I_, j, i + 1) =
          MAT(ekf7_tmp_, j, 0) * SYM(ekf7_I_, i + 1, i + 1);
    }
    // B = inv(A) + w^T*v
    for (int j = 0; j <= i; j++) {
      for (int k = j; k <= i; k++) {
        SYM(ekf7_I_, j, k) =
            SYM(ekf7_I_, j, k) + MAT(ekf7_tmp_, j, 0) * SYM(ekf7_I_, k, i + 1);
      }
    }
  }
  return;
}

void ekf7_KH_update(void) {
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = 0; j < NUM_ST; j++) {
      // printf("\r\n-------- KH update (%d, %d) --------\r\n", i, j);
      // printSyM(ekf7_covar_);
      MAT(ekf7_KH_, i, j) = 0;
      for (int k = 0; k < NUM_OB; k++) {
        // printf("\r\n-------- KH update K: %d --------\r\n", k);
        // printSyM(ekf7_covar_);
        if ((j >= 3) && (j < NUM_ST)) {
          MAT(ekf7_KH_, i, j) +=
              SYM(ekf7_covar_, i, k + 3) * SYM(ekf7_I_, k, j - 3);
        }
      }
    }
  }
}

// P = (I-KH)P
void ekf7_covariance_update(void) {
  // EKF_DBG("\r\n-------- Covariance --------\r\n")
  // printSyM(ekf7_covar_);
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = i; j < NUM_ST; j++) {
      FTYPE Pij = SYM(ekf7_covar_, i, j);
      SYM(ekf7_covar_, i, j) = 0;
      for (int k = 0; k < NUM_ST; k++) {
        if (i == k) {
          SYM(ekf7_covar_, i, j) += (1 - MAT(ekf7_KH_, i, k)) * Pij;
        } else {
          SYM(ekf7_covar_, i, j) +=
              (-MAT(ekf7_KH_, i, k)) * SYM(ekf7_covar_, k, j);
        }
      }
    }
  }
}

// x = x + K(y-Hx) = x + KH(x_err)
void ekf7_state_update(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz) {
  // error
  MAT(ekf7_error_, 0, 0) = 0 - BX;
  MAT(ekf7_error_, 1, 0) = 0 - BY;
  MAT(ekf7_error_, 2, 0) = 0 - BZ;
  MAT(ekf7_error_, 3, 0) = qw - QW;
  MAT(ekf7_error_, 4, 0) = qx - QX;
  MAT(ekf7_error_, 5, 0) = qy - QY;
  MAT(ekf7_error_, 6, 0) = qz - QZ;
  // update
  for (int i = 0; i < NUM_ST; i++) {
    // MAT(ekf7_state_, i, 0) = 0;
    for (int j = 0; j < NUM_ST; j++) {
      MAT(ekf7_state_, i, 0) += MAT(ekf7_KH_, i, j) * MAT(ekf7_error_, j, 0);
    }
  }
}

void ekf7_update(FTYPE qw, FTYPE qx, FTYPE qy, FTYPE qz) {
  // get KH:
  // EKF_DBG("\r\n-------- before update inno. --------\r\n")
  // printSyM(ekf7_I_);
  ekf7_innovation_update();
  float bfi[7][7];
  for (int i = 0; i < NUM_OB; i++) {
    for (int j = 0; j < NUM_OB; j++) {
      bfi[i][j] = SYM(ekf7_I_, i, j);
    }
  }
  ekf7_innovation_invert();
  double sbe[7][7];
  ekf7_KH_update();
  ekf7_state_update(qw, qx, qy, qz);
  ekf7_covariance_update();
}

/*  */
void ekf7_get_gyroscope_bias(FTYPE *bx, FTYPE *by, FTYPE *bz) {
  *bx = BX;
  *by = BY;
  *bz = BZ;
}

/*  */
void ekf7_get_quaternion(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  *qw = QW;
  *qx = QX;
  *qy = QY;
  *qz = QZ;
}

void ekf7_get_rpy(FTYPE *roll, FTYPE *pitch, FTYPE *yaw) {
  *roll = atan2(2 * (QW * QX + QY * QZ), 1 - 2 * (QX * QX + QY * QY));
  *pitch = asin(2 * (QW * QY - QZ * QX));
  *yaw = atan2(2 * (QW * QZ + QX * QY), 1 - 2 * (QY * QY + QZ * QZ));
}

#ifdef DBG_CODE
void printMat(M x) {
  if (!PRT_DBG)
    return;
  for (int i = 0; i < GET_M(x.dim); i++) {
    for (int j = 0; j < GET_N(x.dim); j++) {
      printf("|%f", MAT(x, i, j));
    }
    printf("|\r\n");
  }
}

void printSyM(M x) {
  if (!PRT_DBG)
    return;
  printf("[");
  for (int i = 0; i < GET_M(x.dim); i++) {
    printf("[");
    for (int j = 0; j < GET_N(x.dim); j++) {
      printf("%f,", SYM(x, i, j));
    }
    printf("],\r\n");
  }
  printf("]\r\n");
}

void printF(void) {
  if (!PRT_DBG)
    return;
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = 0; j < NUM_ST; j++) {
      printf("|%f", ekf7_get_F(0.01, i, j));
    }
    printf("|\r\n");
  }
}

void printState(void) {
  for (int i = 0; i < NUM_ST; i++) {
    printf("%f,", MAT(ekf7_state_, i, 0));
  }
  printf("\r\n");
}
#endif

#ifdef __cplusplus
}
#endif
