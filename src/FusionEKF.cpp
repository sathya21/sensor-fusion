#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

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
  ekf_.x_ = VectorXd(4);

	//the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
			 0, 1, 0, 1,
			 0, 0, 1, 0,
			 0, 0, 0, 1;

  ekf_.P_ =  MatrixXd(4, 4);

  ekf_.P_ << 1, 0, 0, 0,
			 0, 1, 0, 0,
			 0, 0, 1000, 0,
			 0, 0, 0, 1000;

  //measurement matrix
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.H_ << 1, 0, 0, 0,
			  0, 1, 0, 0;


  ekf_.Q_ =  MatrixXd(4, 4);

  ekf_.R_ = MatrixXd(2, 2);
  ekf_.R_ << 0.0225, 0,
			 0, 0.0225;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    	float x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
    	float y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
    	ekf_.x_ << x, y, 0, 0;

      }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    }


    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/


    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  	previous_timestamp_ = measurement_pack.timestamp_;

  	float dt_2 = dt * dt;
  	float dt_3 = dt_2 * dt;
  	float dt_4 = dt_3 * dt;

  	//Modify the F matrix so that the time is integrated
  	ekf_.F_(0, 2) = dt;
  	ekf_.F_(1, 3) = dt;

  	float noise_ax = 9;
  	float noise_ay = 9;

  	ekf_.Q_ = MatrixXd(4, 4);
  	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
  			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
  			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/



  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  float x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
	  float y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
	  VectorXd tempx  = VectorXd(3);
	  tempx << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],measurement_pack.raw_measurements_[2];

	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(tempx);

  } else {
    // Laser updates
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
