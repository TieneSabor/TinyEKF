#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#ifndef __linux__
#include "WProgram.h"
#endif
#endif

#include "ekf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct matrix {
  byte dim; // dimision of the matrix, dim = m<<4 + n
  byte sps; // start position in the_space
} M;

// Definition of important constant.
// Should NOT be modified.
#define SPACE_SIZE 245 // maximum matrix element number
#define MAX_DIM 15     // maximum matrix operation
#define LAST4 0x0f     // 00001111
#define NUM_ST 10      // state number
#define NUM_OB 7       // observation number

// Byte operation. For "dim" in matrix struct, it is composed of M(first 4 bit)
// and N(last 4 bit)
#define GET_M(x) (x >> 4)
#define GET_N(x) (x & 0x0f)
#define SET_B(m, n) ((m << 4) + n)

// The (x, y)th element for NOT symmetric matrix of m rows * n columns, starting
// with (0, 0). With O2 Opt., the add/mul. part in the index should be replaced
// by const. (Hopefully?)
#define MAT(mat, x, y) ekf_space_[mat.sps + y + x * (mat.dim & 0x0f)]
// Get the (x, y)th element of a symmetric matrix.
#define SYM(mat, x, y)                                                         \
  ekf_space_[mat.sps + imax(x, y) + imin(x, y) * (mat.dim & 0x0f) -            \
             (imin(x, y) * imin(x, y) + imin(x, y)) / 2]
// Get how much a NOT symmetric matrix take to store data
#define OCU(mat) (mat.dim >> 4) * (mat.dim & 0x0f)
// Get how much a symmetric matrix take to store data
#define OCS(mat)                                                               \
  0.5 * ((mat.dim >> 4) * ((mat.dim & 0x0f) - 1)) + (mat.dim >> 4)

// Some special variable stored in space
#define WX ekf_space_[state_.sps + 0]
#define WY ekf_space_[state_.sps + 1]
#define WZ ekf_space_[state_.sps + 2]
#define Q1 ekf_space_[state_.sps + 6]
#define Q2 ekf_space_[state_.sps + 7]
#define Q3 ekf_space_[state_.sps + 8]
#define Q4 ekf_space_[state_.sps + 9]

// Every element of matrix in the EKF will be stored here.
FTYPE ekf_space_[SPACE_SIZE];

// Matrix used for EKF estimation
M state_;
M error_;
M covar_;
M W_;
M V_;
M I_;
M KH_;
M tmp_;

// Some flags for checking initialization
byte init_called_ = 0;
byte PN_set_ = 0;
byte ON_set_ = 0;

// Important parameters
#define JACSIZ 10 // size used for jacobian matrix
#define INVTAUX ekf_space_[0]
#define INVTAUY ekf_space_[1]
#define INVTAUZ ekf_space_[2]
#define FQ1 ekf_space_[3]
#define FQ2 ekf_space_[4]
#define FQ3 ekf_space_[5]
#define FQ4 ekf_space_[6]
#define FWX ekf_space_[7]
#define FWY ekf_space_[8]
#define FWZ ekf_space_[9]

byte spscnt_;

// Set linearlize transition matrix F.
void set_F(int init) {
  if (!init_called_) {
    return;
  }
  if (init) {
    // Gyroscope latency.
    INVTAUX = -1 / TAUX;
    INVTAUY = -1 / TAUY;
    INVTAUZ = -1 / TAUZ;
  }
  // Set up quaternion in F
  FQ1 = Q1;
  FQ2 = Q2;
  FQ3 = Q3;
  FQ4 = Q4;
  // Set up angular velocity in F
  FWX = WX;
  FWY = WY;
  FWZ = WZ;
}

FTYPE get_F(FTYPE dt, byte i, byte j) {
  FTYPE res = 0;
  // q1
  if ((i == (j + 7)) && (i > 6)) {
    res = (ekf_space_[3] * 0.5);
  }
  // q2
  else if ((j == 9) && (i == 1)) {
    res = (ekf_space_[4] * 0.5);
  }
  // -q2
  else if (((i == 6) && (j == 0)) || ((i == 8) && (j == 2))) {
    res = -(ekf_space_[4] * 0.5);
  }
  // q3
  else if ((i == 7) && (j == 2)) {
    res = (ekf_space_[5] * 0.5);
  }
  // -q3
  else if (((i == 6) && (j == 1)) || ((i == 9) && (j == 0))) {
    res = -(ekf_space_[5] * 0.5);
  }
  // q4
  else if ((i == 8) && (j == 0)) {
    res = (ekf_space_[6] * 0.5);
  }
  // -q4
  else if (((i == 6) && (j == 2)) || ((i == 7) && (j == 1))) {
    res = -(ekf_space_[6] * 0.5);
  }
  // wx
  else if (((i == 7) && (j == 6)) || ((i == 8) && (j == 9))) {
    res = (ekf_space_[7] * 0.5);
  }
  // -wx
  else if (((i == 6) && (j == 7)) || ((i == 9) && (j == 8))) {
    res = -(ekf_space_[7] * 0.5);
  }
  // wy
  else if (((i == 9) && (j == 7)) || ((i == 8) && (j == 6))) {
    res = (ekf_space_[8] * 0.5);
  }
  // -wy
  else if (((i == 7) && (j == 9)) || ((i == 6) && (j == 8))) {
    res = -(ekf_space_[8] * 0.5);
  }
  // wz
  else if (((i == 7) && (j == 8)) || ((i == 9) && (j == 6))) {
    res = (ekf_space_[9] * 0.5);
  }
  // -wz
  else if (((i == 8) && (j == 7)) || ((i == 6) && (j == 9))) {
    res = -(ekf_space_[9] * 0.5);
  }
  // others are 0
  else {
    res = 0;
  }
  // return I + Fdt
  res *= dt;
  return res;
}

FTYPE get_F_plus_I(FTYPE dt, byte i, byte j) {
  if (i == j) {
    return 1 + get_F(dt, i, j);
  } else {
    return get_F(dt, i, j);
  }
}

FTYPE get_H(byte i, byte j) {
  if ((j < 3) && (i == j)) {
    return 1.0f;
  } else if ((j >= 3) && ((i + 3) == j)) {
    return 1.0f;
  } else {
    return 0.0f;
  }
}

void ekf_init(FTYPE wx, FTYPE wy, FTYPE wz) {
  init_called_ = 1;
  // Fill all elements with 0
  for (int i = 0; i < SPACE_SIZE; i++) {
    ekf_space_[i] = 0;
  }
  // Initialize all structures.
  // State, a vector
  // int spscnt = 0 + JACSIZ;
  spscnt_ = 0 + JACSIZ;
  state_.dim = SET_B(NUM_ST, 1);
  state_.sps = spscnt_;
  spscnt_ += OCU(state_);
  // Covariance, a symmetric matrix
  covar_.dim = SET_B(NUM_ST, NUM_ST);
  covar_.sps = spscnt_;
  spscnt_ += OCS(covar_);
  // Propogation niose, a diagnal matrix
  W_.dim = SET_B(NUM_ST, 1);
  W_.sps = spscnt_;
  spscnt_ += OCU(W_);
  // Observation noise, a diagnal matrix
  V_.dim = SET_B(NUM_OB, 1);
  V_.sps = spscnt_;
  spscnt_ += OCU(V_);
  // Innovation matrix, a symmetric matrix
  I_.dim = SET_B(NUM_OB, NUM_OB);
  I_.sps = spscnt_;
  spscnt_ += OCS(I_);
  // KH = PH^T(I^-1)H, not a symmetric matrix
  KH_.dim = SET_B(NUM_ST, NUM_ST);
  KH_.sps = spscnt_;
  spscnt_ += OCU(KH_);
  // error = (H^+)y - x
  error_.dim = SET_B(NUM_ST, 1);
  error_.sps = spscnt_;
  spscnt_ += OCU(error_);
  // for matrix multiplation
  tmp_.dim = SET_B(NUM_ST, 1);
  tmp_.sps = spscnt_;
  spscnt_ += OCU(tmp_);
  // Initialize state with first angular speed measurement
  MAT(state_, 0, 0) = wx;
  MAT(state_, 1, 0) = wy;
  MAT(state_, 2, 0) = wz;
  MAT(state_, 6, 0) = 1; // qw = 1;
  // Initialize covariance P.
  // for (int i = 0; i < NUM_ST; i++) {
  //  SYM(covar_, i, i) = 1.0;
  //}
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = 0; j < NUM_ST; j++) {
      SYM(covar_, i, j) = 1e-2;
    }
  }
  // Initialize Jacobian F.
  set_F(1);
}

void set_Propagation_Noise(FTYPE Wwx, FTYPE Wwy, FTYPE Wwz, FTYPE Wbx,
                           FTYPE Wby, FTYPE Wbz, FTYPE Wq1, FTYPE Wq2,
                           FTYPE Wq3, FTYPE Wq4) {
  if (!init_called_) {
    return;
  }
  // Fill the value
  MAT(W_, 0, 0) = Wwz;
  MAT(W_, 1, 0) = Wwy;
  MAT(W_, 2, 0) = Wwz;
  MAT(W_, 3, 0) = Wbx;
  MAT(W_, 4, 0) = Wby;
  MAT(W_, 5, 0) = Wbz;
  MAT(W_, 6, 0) = Wq1;
  MAT(W_, 7, 0) = Wq2;
  MAT(W_, 8, 0) = Wq3;
  MAT(W_, 9, 0) = Wq4;
  // Set flag
  PN_set_ = 1;
}

void set_Observation_Noise(FTYPE Vwx, FTYPE Vwy, FTYPE Vwz, FTYPE Vq1,
                           FTYPE Vq2, FTYPE Vq3, FTYPE Vq4) {
  if (!init_called_) {
    return;
  }
  // Fill the value
  MAT(V_, 0, 0) = Vwz;
  MAT(V_, 1, 0) = Vwy;
  MAT(V_, 2, 0) = Vwz;
  MAT(V_, 3, 0) = Vq1;
  MAT(V_, 4, 0) = Vq2;
  MAT(V_, 5, 0) = Vq3;
  MAT(V_, 6, 0) = Vq4;
  // Set flag
  ON_set_ = 1;
}

void predict(FTYPE dt) {
  // Update the Jacobian.
  set_F(0);
  // Update state.
  FTYPE nqw = MAT(state_, 6, 0) + 0.5 * dt * (-WX * Q2 - WY * Q3 - WZ * Q4);
  FTYPE nqx = MAT(state_, 7, 0) + 0.5 * dt * (WX * Q1 + WZ * Q3 - WY * Q4);
  FTYPE nqy = MAT(state_, 8, 0) + 0.5 * dt * (WY * Q1 - WZ * Q2 + WX * Q4);
  FTYPE nqz = MAT(state_, 9, 0) + 0.5 * dt * (WZ * Q1 + WY * Q2 - WX * Q3);
  FTYPE norm = sqrt(pow(nqw, 2) + pow(nqx, 2) + pow(nqy, 2) + pow(nqz, 2));
  //FTYPE norm = 1;
  MAT(state_, 6, 0) = nqw / norm;
  MAT(state_, 7, 0) = nqx / norm;
  MAT(state_, 8, 0) = nqy / norm;
  MAT(state_, 9, 0) = nqz / norm;
  // Update the covariance
  // P = FPF^T + W
  covariance_predict(dt);
}

void covariance_predict(FTYPE dt) {
  for (int i = 0; i < NUM_OB; i++) {
    for (int j = i; j < NUM_OB; j++) {
      MAT(tmp_, j, 0) = 0;
      if (i == j) {
        MAT(tmp_, j, 0) = MAT(W_, i, 0);
      }
      for (int k = 0; k < NUM_ST; k++) {
        for (int l = 0; l < NUM_ST; l++) {
          MAT(tmp_, j, 0) += get_F_plus_I(dt, i, k) * SYM(covar_, k, l) *
                             get_F_plus_I(dt, j, l);
        }
      }
    }
    for (int j = i; j < NUM_ST; j++) {
      SYM(covar_, i, j) = MAT(tmp_, j, 0);
    }
  }
}

void innovation_update(void) {
  for (int i = 0; i < NUM_OB; i++) {
    for (int j = i; j < NUM_OB; j++) {
      MAT(tmp_, j, 0) = 0;
      if (i == j) {
        MAT(tmp_, j, 0) = MAT(V_, i, 0);
      }
      for (int k = 0; k < NUM_ST; k++) {
        for (int l = 0; l < NUM_ST; l++) {
          MAT(tmp_, j, 0) += get_H(i, k) * SYM(covar_, k, l) * get_H(j, l);
        }
      }
    }
    for (int j = i; j < NUM_ST; j++) {
      SYM(I_, i, j) = MAT(tmp_, j, 0);
    }
  }
}

// Matrix Inversion Routine
void innovation_invert(void) {
  SYM(I_, 0, 0) = 1.0f / SYM(I_, 0, 0);
  for (int i = 0; i < NUM_OB - 1; i++) {
    // w^T = -inv(A)*u^T
    // Here, the submatrix up to I(i, i) (that is, A) was already inverted.
    for (int j = 0; j <= i; j++) {
      MAT(tmp_, j, 0) = 0;
      for (int k = 0; k <= i; k++) {
        MAT(tmp_, j, 0) -= SYM(I_, j, k) * SYM(I_, k, i + 1);
      }
    }
    // y = inv(u*w^T+x)
    // Firstly u*w^T:
    FTYPE y = 0;
    for (int j = 0; j <= i; j++) {
      y += SYM(I_, j, i + 1) * MAT(tmp_, j, 0);
    }
    SYM(I_, i + 1, i + 1) = 1.0f / (y + SYM(I_, i + 1, i + 1));
    // v^T = w^T * y
    for (int j = 0; j <= i; j++) {
      SYM(I_, j, i + 1) = MAT(tmp_, j, 0) * SYM(I_, i + 1, i + 1);
    }
    // B = inv(A) + w^T*v
    for (int j = 0; j <= i; j++) {
      for (int k = j; k <= i; k++) {
        SYM(I_, j, k) = SYM(I_, j, k) + MAT(tmp_, j, 0) * SYM(I_, k, i + 1);
      }
    }
  }
  return;
}

void KH_update(void) {
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = 0; j < NUM_ST; j++) {
      // printf("\r\n-------- KH update (%d, %d) --------\r\n", i, j);
      // printSyM(covar_);
      MAT(KH_, i, j) = 0;
      for (int k = 0; k < 7; k++) {
        // printf("\r\n-------- KH update K: %d --------\r\n", k);
        // printSyM(covar_);
        if ((j >= 0) && (j < 3)) {
          if ((k >= 0) && (k < 3)) {
            MAT(KH_, i, j) +=
                (SYM(covar_, i, k) + SYM(covar_, i, k + 3)) * SYM(I_, k, j);
          } else if ((k >= 3) && (k < 7)) {
            MAT(KH_, i, j) += SYM(covar_, i, k + 3) * SYM(I_, k, j);
          }
        } else if ((j >= 3) && (j < 10)) {
          if ((k >= 0) && (k < 3)) {
            MAT(KH_, i, j) +=
                (SYM(covar_, i, k) + SYM(covar_, i, k + 3)) * SYM(I_, k, j - 3);
          } else if ((k >= 3) && (k < 7)) {
            MAT(KH_, i, j) += SYM(covar_, i, k + 3) * SYM(I_, k, j - 3);
          }
        }
      }
    }
  }
}

// P = (I-KH)P
void covariance_update(void) {
  // EKF_DBG("\r\n-------- Covariance --------\r\n")
  // printSyM(covar_);
  for (int i = 0; i < NUM_ST; i++) {
    for (int j = i; j < NUM_ST; j++) {
      FTYPE Pij = SYM(covar_, i, j);
      SYM(covar_, i, j) = 0;
      for (int k = 0; k < NUM_ST; k++) {
        if (i == k) {
          SYM(covar_, i, j) += (1 - MAT(KH_, i, k)) * Pij;
        } else {
          SYM(covar_, i, j) += (-MAT(KH_, i, k)) * SYM(covar_, k, j);
        }
      }
    }
  }
}

// x = x + K(y-Hx) = x + KH(x_err)
void state_update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx, FTYPE qy,
                  FTYPE qz) {
  // error
  MAT(error_, 0, 0) = wx - MAT(state_, 0, 0);
  MAT(error_, 1, 0) = wy - MAT(state_, 1, 0);
  MAT(error_, 2, 0) = wz - MAT(state_, 2, 0);
  MAT(error_, 3, 0) = 0 - MAT(state_, 3, 0);
  MAT(error_, 4, 0) = 0 - MAT(state_, 4, 0);
  MAT(error_, 5, 0) = 0 - MAT(state_, 5, 0);
  MAT(error_, 6, 0) = qw - MAT(state_, 6, 0);
  MAT(error_, 7, 0) = qx - MAT(state_, 7, 0);
  MAT(error_, 8, 0) = qy - MAT(state_, 8, 0);
  MAT(error_, 9, 0) = qz - MAT(state_, 9, 0);
  // update
  for (int i = 0; i < NUM_ST; i++) {
    // MAT(state_, i, 0) = 0;
    for (int j = 0; j < NUM_ST; j++) {
      MAT(state_, i, 0) += MAT(KH_, i, j) * MAT(error_, j, 0);
    }
  }
}

void update(FTYPE wx, FTYPE wy, FTYPE wz, FTYPE qw, FTYPE qx, FTYPE qy,
            FTYPE qz) {
  // get KH:
  // EKF_DBG("\r\n-------- before update inno. --------\r\n")
  // printSyM(I_);
  innovation_update();
  float bfi[7][7];
  for (int i = 0; i < NUM_OB; i++) {
    for (int j = 0; j < NUM_OB; j++) {
      bfi[i][j] = SYM(I_, i, j);
    }
  }
  innovation_invert();
  double sbe[7][7];
  KH_update();
  state_update(wx, wy, wz, qw, qx, qy, qz);
  covariance_update();
}

/*  */
void get_angular_vel(FTYPE *wx, FTYPE *wy, FTYPE *wz) {
  *wx = WX;
  *wy = WY;
  *wz = WZ;
}

/*  */
void get_gyroscope_bias(FTYPE *bx, FTYPE *by, FTYPE *bz) {
  *bx = MAT(state_, 3, 0);
  *by = MAT(state_, 4, 0);
  *bz = MAT(state_, 5, 0);
}

/*  */
void get_quaternion(FTYPE *qw, FTYPE *qx, FTYPE *qy, FTYPE *qz) {
  *qw = Q1;
  *qx = Q2;
  *qy = Q3;
  *qz = Q4;
}

void get_ekf_rpy(FTYPE *yaw, FTYPE *pitch, FTYPE *roll) {
  *yaw = atan2(2 * (Q1 * Q2 + Q3 * Q4), 1 - 2 * (Q2 * Q2 + Q3 * Q3));
  *pitch = asin(2 * (Q1 * Q3 - Q4 * Q2));
  *roll = atan2(2 * (Q1 * Q4 + Q2 * Q3), 1 - 2 * (Q3 * Q3 + Q4 * Q4));
}

#ifdef __cplusplus
}
#endif
