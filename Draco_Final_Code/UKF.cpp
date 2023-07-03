#include <BasicLinearAlgebra.h>
#include <ElementStorage.h>
#include <math.h>

using namespace BLA;

//UKF
int kappa = 0;  //All values are multiplied by 10^6
float alfa = 0.005;
int beta = 6;
double lambda = (6 + kappa) * sq(alfa) - 6;
double W0m = lambda / (6 + lambda);
double W0c = lambda / (6 + lambda) + (1 - sq(alfa) + beta);
double W = 1 / (2 * (6 + lambda));


extern ArrayMatrix<6, 1, double> x_aposteriori;
extern ArrayMatrix<3, 1, double> z;
extern double timeInt;

ArrayMatrix<6, 1, double> x_apriori = {5, 0, 0, 0, 0, 0 };

ArrayMatrix<3, 1, double> z_apriori = 0;
ArrayMatrix<3, 1, double> v = 0;

ArrayMatrix<7, 13, double> x_sigma = 0;
ArrayMatrix<7, 13, double> y_sigma = 0;
ArrayMatrix<3, 13, double> z_sigma = 0;

ArrayMatrix<6, 6, double> P_apriori = 0;
ArrayMatrix<6, 6, double> P_aposteriori;
ArrayMatrix<6, 6, double> sP_aposteriori = 0;
ArrayMatrix<6, 12, double> matrix_W = 0;
ArrayMatrix<6, 13, double> Wi_prime = 0;
ArrayMatrix<3, 3, double> P_vv = 0;
ArrayMatrix<3, 3, double> P_zz = 0;
ArrayMatrix<6, 3, double> P_xz = 0;
ArrayMatrix<6, 3, double> K = 0;

ArrayMatrix<6, 6, double> Q;
ArrayMatrix<3, 3, double> R;

double rate_sensornoise = 0.033;
void UKF_setup() {
  Q.Fill(0);
  R.Fill(0);
  P_aposteriori.Fill(0);
  for (int i = 0; i < 6; i++) {
    Q(i, i) = 0.03;
  }
  for (int i = 0; i < 4; i++) {
    R(i, i) = rate_sensornoise;
  }
}

BLA::ArrayMatrix<4, 1, double> quatinv(BLA::ArrayMatrix<4, 1, double> q1) {
  double quatmagnsq = sq(q1(0)) + sq(q1(1)) + sq(q1(2)) + sq(q1(3));

  BLA::ArrayMatrix<4, 1, double> q2 = { q1(0) / quatmagnsq, -q1(1) / quatmagnsq, -q1(2) / quatmagnsq, -q1(3) / quatmagnsq };
  return q2;
}

BLA::ArrayMatrix<4, 1, double> quatnorm(BLA::ArrayMatrix<4, 1, double> q1) {
  double quatmagn = sqrt(sq(q1(0)) + sq(q1(1)) + sq(q1(2)) + sq(q1(3)));
  BLA::ArrayMatrix<4, 1, double> q2 = { q1(0) / quatmagn, q1(1) / quatmagn, q1(2) / quatmagn, q1(3) / quatmagn };
  return q2;
}

BLA::ArrayMatrix<4, 1, double> quatmult(BLA::ArrayMatrix<4, 1, double> q1, BLA::ArrayMatrix<4, 1, double> q2) {
  BLA::ArrayMatrix<4, 1, double> q3 = { q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3),
                                        q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
                                        q1(0) * q2(2) - q1(1) * q2(3) + q1(2) * q2(0) + q1(3) * q2(1),
                                        q1(0) * q2(3) + q1(1) * q2(2) - q1(2) * q2(1) + q1(3) * q2(0) };
  return q3;
}

BLA::ArrayMatrix<4, 1, double> toQuat(BLA::ArrayMatrix<3, 1, double> AxisAngle) {
  double angle = sqrt(sq(AxisAngle(0)) + sq(AxisAngle(1)) + sq(AxisAngle(2)));
  BLA::ArrayMatrix<4, 1, double> Quat = { cos(radians(angle / 2)), AxisAngle(0) * sin(radians(angle / 2)), AxisAngle(1) * sin(radians(angle / 2)), AxisAngle(2) * sin(radians(angle / 2)) };
  if (angle != 0) {
    Quat.Submatrix<3, 1>(1, 0) = Quat.Submatrix<3, 1>(1, 0) / angle;
  }
  return Quat;
}


BLA::ArrayMatrix<3, 1, double> toAxisAngle(BLA::ArrayMatrix<4, 1, double> Quat) {
  double angle = 2 * acos(Quat(0));
  BLA::ArrayMatrix<3, 1> AxisAngle = { Quat(1) * degrees(angle), Quat(2) * degrees(angle), Quat(3) * degrees(angle) };

  if (angle != 0) {
    AxisAngle = AxisAngle / sin(angle / 2);
  }
  return AxisAngle;
}


using namespace BLA;
void setProcessCovariance() {
  double attitude_noise = (0.5 * timeInt);
  for (int i = 0; i < 3; i++) {
    Q(i, i) = 0.1;  //x_aposteriori(i,0) * attitude_noise; //don't know what happened here
  }
}

void meas_meancov() {
  z_apriori = z_sigma.Submatrix<3, 1>(0, 0) * W0m;
  for (int i = 1; i < 13; i++) {
    z_apriori += z_sigma.Submatrix<3, 1>(0, i) * W;
  }
  P_zz = ((z_sigma.Submatrix<3, 1>(0, 0) - z_apriori) * ~(z_sigma.Submatrix<3, 1>(0, 0) - z_apriori)) * W0c;
  for (int i = 1; i < 13; i++) {
    P_zz += ((z_sigma.Submatrix<3, 1>(0, i) - z_apriori) * ~(z_sigma.Submatrix<3, 1>(0, i) - z_apriori)) * W;
  }
  P_vv = P_zz + R;
}



void ap_meancov() {
  ArrayMatrix<3, 13, double> rot_ei;  //big matrix to store
  ArrayMatrix<3, 1, double> rot_e;

  ArrayMatrix<4, 1, double> qt = toQuat(x_aposteriori.Submatrix<3, 1>(0, 0));

  do {
    ArrayMatrix<4, 1, double> qtinv = quatinv(qt);
    rot_e.Fill(0);

    for (int i = 0; i < 13; i++) {
      ArrayMatrix<4, 1, double> ei = quatmult(y_sigma.Submatrix<4, 1>(0, i), qtinv);
      ei = quatnorm(ei);
      rot_ei.Submatrix<3, 1>(0, i) = toAxisAngle(ei);

      rot_e += rot_ei.Submatrix<3, 1>(0, i) / 13;
    }

    ArrayMatrix<4, 1, double> qe = toQuat(rot_e);

    qt = quatmult(qe, qt);
    qt = quatnorm(qt);

  } while (abs(rot_e(0)) > 0.01);

  ArrayMatrix<4, 1, double> qtinv = quatinv(qt);

  ArrayMatrix<3, 1, double> rot_mean = toAxisAngle(qt);

  ArrayMatrix<3, 1, double> angmean = y_sigma.Submatrix<3, 1>(4, 0) * W0m;
  for (int i = 1; i < 13; i++) {
    angmean += y_sigma.Submatrix<3, 1>(4, i) * W;
  }

  ArrayMatrix<3, 13, double> angdiff;
  for (int i = 0; i < 13; i++) {
    angdiff.Submatrix<3, 1>(0, i) = y_sigma.Submatrix<3, 1>(4, i) - angmean;
  }

  ArrayMatrix<3, 13, double> rotdiff;
  for (int i = 0; i < 13; i++) {
    rotdiff.Submatrix<3, 1>(0, i) = toAxisAngle(quatmult(y_sigma.Submatrix<4, 1>(0, i), qtinv));
  }

  x_apriori = rot_mean && angmean;  //vertical concatenate

  Wi_prime = rot_ei && angdiff;

  P_apriori = (Wi_prime.Submatrix<6, 1>(0, 0) * ~Wi_prime.Submatrix<6, 1>(0, 0)) * W0c;
  for (int i = 1; i < 13; i++) {
    P_apriori += (Wi_prime.Submatrix<6, 1>(0, i) * ~Wi_prime.Submatrix<6, 1>(0, i)) * W;
  }
}

void cross_correlation() {
  P_xz = (Wi_prime.Submatrix<6, 1>(0, 0) * ~(z_sigma.Submatrix<3, 1>(0, 0) - z_apriori)) * W0c;
  for (int i = 1; i < 13; i++) {
    P_xz += (Wi_prime.Submatrix<6, 1>(0, i) * ~(z_sigma.Submatrix<3, 1>(0, i) - z_apriori)) * W;
  }
}

void process_model() {
  ArrayMatrix<4, 1, double> qk = toQuat(x_aposteriori.Submatrix<3, 1>(0, 0));
  ArrayMatrix<4, 1, double> qw;
  x_sigma.Submatrix<7, 1>(0, 0) = qk && x_aposteriori.Submatrix<3, 1>(3, 0);

  for (int i = 1; i < 13; i++) {
    qw = toQuat(matrix_W.Submatrix<3, 1>(0, i - 1));
    qw = quatnorm(qw);
    x_sigma.Submatrix<7, 1>(0, i) = quatmult(qk, qw) && (x_aposteriori.Submatrix<3, 1>(3, 0) + matrix_W.Submatrix<3, 1>(3, i - 1));
  }
}

void timestep() {
  ArrayMatrix<4, 1> qgrad;
  for (int i = 0; i < 13; i++) {
    qgrad = toQuat(x_sigma.Submatrix<3, 1>(4, i) * timeInt);
    qgrad = quatnorm(qgrad);
    y_sigma.Submatrix<7, 1>(0, i) = quatmult(x_sigma.Submatrix<4, 1>(0, i), qgrad) && x_sigma.Submatrix<3, 1>(4, i);
  }
}

void UKF_full_step() {
  setProcessCovariance();

  BLA::ArrayMatrix<6, 6, double> x = (P_aposteriori + Q) * (6 + lambda);
  sP_aposteriori.Fill(0);
  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < j - 1; i++) {
      double sum = x(i, j);

      for (int k = 0; k < i - 1; k++) {
        sum -= sP_aposteriori(k, i) * sP_aposteriori(k, j);
      }
      sP_aposteriori(i, j) = sum / sP_aposteriori(i, i);
    }
    double sum = x(j, j);
    for (int k = 0; k < j - 1; k++) {
      sum -= sP_aposteriori(k, j) * sP_aposteriori(k, j);
    }
    sP_aposteriori(j, j) = sqrt(sum);
  }



  for (int i = 0; i < 6; i++) {
    matrix_W.Submatrix<6, 1>(0, i) = sP_aposteriori.Submatrix<6, 1>(0, i);
    matrix_W.Submatrix<6, 1>(0, i + 6) = sP_aposteriori.Submatrix<6, 1>(0, i) * (-1);
  }

  process_model();
  timestep();
  ap_meancov();

  z_sigma = y_sigma.Submatrix<3, 13>(4, 0);  //measurement model

  meas_meancov();

  v = z - z_apriori;

  cross_correlation();

  auto P_vv_inv = P_vv;
  Invert(P_vv_inv);
  K = P_xz * P_vv_inv;

  x_aposteriori = x_apriori + K * v;

  P_aposteriori = P_apriori - K * P_vv * ~K;
}


