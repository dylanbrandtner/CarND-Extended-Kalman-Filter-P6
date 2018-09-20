#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

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
    
    // measurement matrix - laser
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
                
    // measurement matrix - radar
    Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1; 

    // create a 4D state vector, we don't know yet the values of the x state
    VectorXd x = VectorXd(4);
    x << 1,1,1,1;

    //state covariance matrix P
    MatrixXd P = MatrixXd(4, 4);
    P << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;


    //the initial transition matrix F
    MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
          
    // initial process covariance matrix Q
    MatrixXd Q = MatrixXd(4, 4);
    Q << 1, 0, 1, 0,
          0, 1, 0, 1,
          1, 0, 1, 0,
          0, 1, 0, 1;

    
    ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
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
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        
        // store timestamp
        previous_timestamp_ = measurement_pack.timestamp_;

        // Initialize state with first measurements
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
          // Convert radar from polar to cartesian coordinates and initialize state.
          float range = measurement_pack.raw_measurements_(0);
          float bearing = measurement_pack.raw_measurements_(1);
          float velocity = measurement_pack.raw_measurements_(2);
          ekf_.x_ <<  range * cos(bearing), range * sin(bearing), velocity * cos(bearing), velocity * sin(bearing);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
          ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
  
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   
    //Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    // Noise for covariance matrix
    int noise_ax = 9;
    int noise_ay = 9;    
    
    // Setup covariance matrix
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << pow(dt,4)/4*noise_ax, 0, pow(dt,3)/2*noise_ax, 0,
              0, pow(dt,4)/4*noise_ay, 0, pow(dt,3)/2*noise_ay,
              pow(dt,3)/2*noise_ax, 0, pow(dt,2)*noise_ax, 0,
              0, pow(dt,3)/2*noise_ay, 0, pow(dt,2)*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Set H/R matrices appropriately and update
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "sensor_type_: " << measurement_pack.sensor_type_ << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
