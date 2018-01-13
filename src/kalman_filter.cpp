#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif
KalmanFilter::KalmanFilter() {}
KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double px(x_(0)),py(x_(1)),vx(x_(2)),vy(x_(3));
  double c1= px*px+py*py;
  double c2 = sqrt(c1);
  VectorXd z_pred(3);
  z_pred(0)=c2;
  z_pred(1)=atan2(py,px);
  z_pred(2)=(px*vx+py*vy)/c2;
  VectorXd y = z - z_pred;
  if(y(1)>0) {
    while(y(1)>M_PI)
      y(1) = y(1) - 2*M_PI;
  } else {
    while(y(1)<-M_PI)
      y(1) = y(1) + 2*M_PI;
  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
}
