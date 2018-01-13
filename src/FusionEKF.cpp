#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ <<1,0,0,0,0,1,0,0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (!is_initialized_) {
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4,4);
    ekf_.F_ = MatrixXd(4,4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double ro(measurement_pack.raw_measurements_(0)),
	theta(measurement_pack.raw_measurements_(1)),
	ro_dot(measurement_pack.raw_measurements_(2));
      double ct(cos(theta)),st(sin(theta));
      ekf_.x_ << ro*ct,ro*st,ro_dot/ct,0;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.P_ = Hj_.transpose()*R_radar_*Hj_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_<<measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),0.0,0.0;
      
      ekf_.P_ = H_laser_.transpose()*R_laser_*H_laser_;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }
  auto dt = 1e-6*(measurement_pack.timestamp_ - previous_timestamp_);

  ekf_.F_ <<
    1,0,dt,0,
    0,1,0,dt,
    0,0,1,0,
    0,0,0,1;
  MatrixXd G = MatrixXd(4,2);
  G<<
    0.5*dt*dt,0,
    0,0.5*dt*dt,
    dt,0,
    0,dt;
  MatrixXd Qv = MatrixXd(2,2);
  Qv <<9,0,0,9;
  ekf_.Q_ = G*Qv*G.transpose();
  ekf_.Predict();
  previous_timestamp_ = measurement_pack.timestamp_;  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_= R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

}
