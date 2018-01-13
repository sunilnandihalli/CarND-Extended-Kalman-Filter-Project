#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:

  Eigen::VectorXd x_;  // state vector
  Eigen::MatrixXd P_;  // state covariance matrix
  Eigen::MatrixXd F_;  // state transition matrix
  Eigen::MatrixXd Q_;  // process covariance matrix
  Eigen::MatrixXd H_;  // measurement matrix
  Eigen::MatrixXd R_;  // measurement covariance matrix
  KalmanFilter();
  virtual ~KalmanFilter();
  void Predict();
  void Update(const Eigen::VectorXd &z /*measurement at k+1*/);
  void UpdateEKF(const Eigen::VectorXd &z);
  void UpdatePosterior(const Eigen::VectorXd &y);
};

#endif /* KALMAN_FILTER_H_ */
