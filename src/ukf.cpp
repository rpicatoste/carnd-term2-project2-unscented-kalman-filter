#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);
	x_.fill(0.0);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);
	P_.fill(0.0);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 3;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	// TODO: Complete the initialization. See ukf.h for other member properties.
	// Hint: one or more values initialized above might be wildly off...
	time_us_ = 0;
	n_x_ = 5;
	n_aug_ = 7;
	lambda_ = 0.0;
	is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack)
{
	/**
	TODO:

	Complete this function! Make sure you switch between lidar and radar
	measurements.
	*/
	if (!is_initialized_) {
		this->FirstMeasurement( measurement_pack );
		return;
	}



	// Update

	// Use the sensor type to perform the update step.
	// Update the state and covariance matrices.
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

		std::cout << "Radar update! = " << measurement_pack.raw_measurements_.transpose() << std::endl;

		this->x_ = this->ConvertRadarToCartesian( measurement_pack.raw_measurements_ );
	}
	else {

		std::cout << "Laser update! = " << measurement_pack.raw_measurements_.transpose() << std::endl;

		this->x_(0) = measurement_pack.raw_measurements_(0);
		this->x_(1) = measurement_pack.raw_measurements_(1);

	}
}


void UKF::FirstMeasurement(MeasurementPackage measurement_pack)
{
	/**
	 * Initialize the state ekf_.x_ with the first measurement.
	 * Initialize the covariance matrix.
	 * Remember: you'll need to convert radar from polar to cartesian coordinates.
	 */

	this->P_ << MatrixXd::Identity( this->P_.rows(), this->P_.rows() );

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Convert radar from polar to cartesian coordinates and initialize state.
		this->x_ = this->ConvertRadarToCartesian( measurement_pack.raw_measurements_ );

	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

		this->x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0, 0;

	}
	else {
		// Wrong sensor
		assert(0);
	}

	this->time_us_ = measurement_pack.timestamp_;

	this->is_initialized_ = true;

	std::cout << "UKF initialized, x = " << this->x_.transpose() << std::endl;
	std::cout << "ukf.P_:\n " << this->P_ << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/
}





VectorXd UKF::ConvertCartesianToRadar(const VectorXd &cartesian)
{
	VectorXd radar_measurement(3);

	float px   = cartesian(0);
	float py   = cartesian(1);
	float v    = cartesian(2);
	float yaw  = cartesian(3);
	float dyaw = cartesian(4);

	radar_measurement(0) = sqrt( px*px + py*py );
	radar_measurement(1) = atan2( py, px );

	if( fabs(radar_measurement(0)) < 1.0e-6 ){
		radar_measurement(2) = 0.0;
		std::cout << "ConvertCartesianToRadar: Near zero division (by " << radar_measurement(0) << ")! - x = " << cartesian.transpose() << std::endl;
	}
	else{
		radar_measurement(2) = (px*v*cos(yaw) + py*v*sin(yaw)) / radar_measurement(0);
	}

	return radar_measurement;
}

VectorXd UKF::ConvertRadarToCartesian(const VectorXd &radar_measurement)
{
	// This function will calculate the cartesian coordinates from the radar coordinate.
	// However, the speed information in radar coordinates is incomplete, and only
	// the position will be returned (the cartesian speeds will be 0).
	VectorXd cartesian( this->n_x_ );

	float rho = radar_measurement(0);
	float phi = radar_measurement(1);

	cartesian(0) = rho * cos(phi);
	cartesian(1) = rho * sin(phi);
	cartesian(2) = 0.0;
	cartesian(3) = 0.0;
	cartesian(4) = 0.0;

	return cartesian;
}

float UKF::LimitAngle(float angle)
{
	// Limit the angular error to +/- PI
	while( angle > M_PI ){
		angle -= 2*M_PI;
	}
	while( angle < -M_PI ){
		angle += 2*M_PI;
	}
	return angle;
}

