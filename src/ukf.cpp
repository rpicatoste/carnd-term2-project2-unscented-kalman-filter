#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define DEBUGGING (1)

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
	std_a_ = 1.2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.6;



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

	// Initialization.
	this->time_us_ = 0;
	this->n_x_ = 5;
	this->n_aug_ = this->n_x_ + 2;
	this->n_sigma_points_ = 2*this->n_aug_+1;
	this->lambda_ = 3.0 - this->n_x_;
	this->is_initialized_ = false;

	this->Xsig_pred_ = MatrixXd( this->n_x_, this->n_sigma_points_);

	//set weights
	this->weights_ = VectorXd( this->n_sigma_points_ );
	this->weights_(0) = this->lambda_ / (this->lambda_ + this->n_aug_);

	for(int ii = 1; ii < this->n_sigma_points_; ii++){
		this->weights_(ii) = 1 / (2*(this->lambda_ + this->n_aug_));
	}

	// Measurement covariance matrices
	R_laser_ = MatrixXd(2,2);
	R_laser_ << this->std_laspx_*this->std_laspx_,  0,
				0, 						    		this->std_laspy_*this->std_laspy_;

	R_radar_ = MatrixXd(3,3);
	R_radar_ << this->std_radr_*this->std_radr_, 0,									  0,
				0, 								 this->std_radphi_*this->std_radphi_, 0,
				0,								 0, 								  this->std_radrd_*this->std_radrd_;

	NIS_laser_ = 0.0;
	NIS_radar_ = 0.0;

}


UKF::~UKF() {}



/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack)
{

	/*****************************************************************************
	 *  Initialization and time calculation.
	 ****************************************************************************/
	this->DebugPrint("ProcessMeasurement: Start ----------------------------------------", measurement_pack.raw_measurements_.transpose() );
	if (!is_initialized_) {
		this->FirstMeasurement( measurement_pack );
		return;
	}

	if ( ( (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) && (this->use_radar_ == false) ) ||
		 ( (measurement_pack.sensor_type_ == MeasurementPackage::LASER) && (this->use_laser_ == false) ) )
	{
		this->DebugPrint("Measurement avoided.");
		// Skip if the sensor used is meant to be ignored.
		return;
	}

	// Compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - this->time_us_) / 1000000.0;	//dt - expressed in seconds

	// If the simulator is restarted, the delta time will be negative, and we need to consider the measurement
	// as the initialization one.
	if(dt < 0.0){
		std::cout << "Restart of the simulation detected: The Kalman filter is being initialized." << std::endl;
		this->FirstMeasurement(measurement_pack);
		return;
	}
	else if (dt > 10.0){
		std::cout << "WARNING: dt was above 10 seconds: " << dt << " sec. The Kalman filter is being initialized." << std::endl;
		this->FirstMeasurement(measurement_pack);
		return;

	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	this->Prediction( dt );

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/
	// Use the sensor type to perform the update step.
	// Update the state and covariance matrices.
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		this->UpdateRadar( measurement_pack );
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		this->UpdateLidar( measurement_pack );
	}
	else{
		assert(0);
	}

	this->time_us_ = measurement_pack.timestamp_;
	this->DebugPrint("ProcessMeasurement: Leaving\n");
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
	// Estimate the object's location. Modify the state vector, x_.
	// Predict sigma points, the state, and the state covariance matrix.

	// Generate sigma points.
	this->DebugPrint("Prediction: Generate sigma points. Initial state:", this->x_.transpose() );

	// Create augmented mean state
	VectorXd x_aug = VectorXd( this->n_aug_ );
	x_aug.head(this->n_x_) = this->x_;
	x_aug.segment(this->n_x_, 2) << 0, 0;

	// Create augmented covariance matrix
	MatrixXd P_aug = MatrixXd( this->n_aug_, this->n_aug_ );
	P_aug.fill(0.0);
	P_aug.topLeftCorner(this->n_x_, this->n_x_) = this->P_;
	P_aug.block(this->n_x_, this->n_x_, 2, 2) << 	this->std_a_*this->std_a_,	0,
													0,							this->std_yawdd_*this->std_yawdd_;

	// Create square root matrix: sqrt( (lambda + n_x) * P )
	MatrixXd A = P_aug.llt().matrixL();
	A = sqrt(this->lambda_ + this->n_aug_) * A;

	// Create augmented sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, this->n_sigma_points_ );
	MatrixXd sigma_points_half_1( this->n_aug_, this->n_aug_);
	MatrixXd sigma_points_half_2( this->n_aug_, this->n_aug_);

	sigma_points_half_1 = x_aug.rowwise().replicate( this->n_aug_ ) + A;
	sigma_points_half_2 = x_aug.rowwise().replicate( this->n_aug_ ) - A;

	Xsig_aug.col(0) = x_aug;
	Xsig_aug.block( 0, 1,         		 this->n_aug_, this->n_aug_) = sigma_points_half_1;
	Xsig_aug.block( 0, this->n_aug_ + 1, this->n_aug_, this->n_aug_) = sigma_points_half_2;

	// Predict sigma points
	this->DebugPrint("Prediction: Predict sigma points");
	for (int ii = 0; ii< this->n_sigma_points_; ii++)
	{
		// Extract values for better readability
		double p_x      = Xsig_aug(0,ii);
		double p_y      = Xsig_aug(1,ii);
		double v        = Xsig_aug(2,ii);
		double yaw      = Xsig_aug(3,ii);
		double yawd     = Xsig_aug(4,ii);
		double nu_a     = Xsig_aug(5,ii);
		double nu_yawdd = Xsig_aug(6,ii);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		// Add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		this->Xsig_pred_(0,ii) = px_p;
		this->Xsig_pred_(1,ii) = py_p;
		this->Xsig_pred_(2,ii) = v_p;
		this->Xsig_pred_(3,ii) = yaw_p;
		this->Xsig_pred_(4,ii) = yawd_p;
	}

	// Predict mean and covariance.
	this->DebugPrint("Prediction: Predict mean.");

	this->x_.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		this->x_ = this->x_ + this->weights_(ii) * this->Xsig_pred_.col(ii);
	}

	this->DebugPrint("Prediction: Predict covariance.");
	this->P_.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
	    // State difference
	    VectorXd x_diff = this->Xsig_pred_.col(ii) - this->x_;
	    x_diff(3) = this->LimitAngle( x_diff(3) );
	    this->P_ = this->P_ + this->weights_(ii) * x_diff * x_diff.transpose() ;
	}
	this->DebugPrint("Prediction: Leave. State predicted: ", this->x_.transpose() );

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack)
{
	// Use lidar data to update the belief about the object's
	// position. Modify the state vector, x_, and covariance, P_.

	this->DebugPrint("Laser update! = ",  measurement_pack.raw_measurements_.transpose() );

	// Predict measurement
	// Size of the measurement vector for lidar.
	int n_z = 2;

	// Create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, this->n_sigma_points_);

	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);

	// Transform sigma points into measurement space.
	this->DebugPrint("UpdateLidar: transform sigma points into measurement space.");
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		Zsig(0, ii) = this->Xsig_pred_(0, ii);
		Zsig(1, ii) = this->Xsig_pred_(1, ii);
	}

	// Calculate mean predicted measurement
	this->DebugPrint("UpdateLidar: Calculate mean predicted measurement.");
	z_pred.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		z_pred = z_pred + this->weights_(ii) * Zsig.col(ii);
	}

	// Calculate measurement covariance matrix S
	this->DebugPrint("UpdateLidar: Calculate measurement covariance matrix S.");

	S.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		// Measurement difference
		VectorXd z_diff = Zsig.col(ii) - z_pred;
		S = S + this->weights_(ii) * z_diff * z_diff.transpose() ;
	}
	S = S + R_laser_;

	// Update state
	this->DebugPrint("UpdateLidar: Update state.");
	VectorXd z = VectorXd(n_z);
	z << measurement_pack.raw_measurements_;

	// Create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd( this->n_x_, n_z );

	// Calculate cross correlation matrix
	this->DebugPrint("UpdateLidar: Calculate cross correlation matrix.");
	Tc.fill( 0.0 );
	for(int ii = 0; ii < this->n_sigma_points_; ii++){

		VectorXd z_diff = Zsig.col(ii) - z_pred;
		VectorXd x_diff = this->Xsig_pred_.col(ii).head(this->n_x_) - this->x_;
		Tc = Tc + this->weights_(ii) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc*S.inverse();

	// Update state mean and covariance matrix
	this->DebugPrint("UpdateLidar: Update state mean and covariance matrix.");
	// Residual
	VectorXd z_diff = z - z_pred;
	this->x_ = this->x_ + K * z_diff;
	this->P_ = this->P_ - K * S * K.transpose();

	this->NIS_laser_ =  z_diff.transpose() * S * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack)
{
//	Complete this function! Use radar data to update the belief about the object's
//	position. Modify the state vector, x_, and covariance, P_.

	this->DebugPrint("Radar update! = ",  measurement_pack.raw_measurements_.transpose() );

	// Predict measurement
	this->DebugPrint("UpdateRadar: Predict measurement.");
	// Size of the measurement vector for radar.
	int n_z = 3;

	// Create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, this->n_sigma_points_);

	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);

	// Transform sigma points into measurement space.
	this->DebugPrint("UpdateRadar: transform sigma points into measurement space.");
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		Zsig.col(ii) = this->ConvertCartesianToRadar( Xsig_pred_.col(ii).head( this->n_x_ ) );
	}

	// Calculate mean predicted measurement
	this->DebugPrint("UpdateRadar: Calculate mean predicted measurement.");
	z_pred.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		z_pred = z_pred + this->weights_(ii) * Zsig.col(ii);
	}

	// Calculate measurement covariance matrix S
	this->DebugPrint("UpdateRadar: Calculate measurement covariance matrix S.");

	S.fill(0.0);
	for(int ii = 0; ii < this->n_sigma_points_; ii++){
		// Measurement difference
		VectorXd z_diff = Zsig.col(ii) - z_pred;
		S = S + this->weights_(ii) * z_diff * z_diff.transpose() ;
	}
	S = S + R_radar_;

	// Update state
	this->DebugPrint("UpdateRadar: Update state.");
	VectorXd z = VectorXd(n_z);
	z << measurement_pack.raw_measurements_;

	// Create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd( this->n_x_, n_z );

	// Calculate cross correlation matrix
	this->DebugPrint("UpdateRadar: Calculate cross correlation matrix.");
	Tc.fill( 0.0 );
	for(int ii = 0; ii < this->n_sigma_points_; ii++){

		VectorXd z_diff = Zsig.col(ii) - z_pred;
		//angle normalization
		z_diff(1) = this->LimitAngle( z_diff(1) );

		VectorXd x_diff = this->Xsig_pred_.col(ii).head(this->n_x_) - this->x_;

		Tc = Tc + this->weights_(ii) * x_diff * z_diff.transpose();

	}

	//calculate Kalman gain K;
	MatrixXd K = Tc*S.inverse();

	// Update state mean and covariance matrix
	this->DebugPrint("UpdateRadar: Update state mean and covariance matrix.");
	// Residual
	VectorXd z_diff = z - z_pred;
	z_diff(1) = this->LimitAngle( z_diff(1) );
	this->x_ = this->x_ + K * z_diff;
	this->P_ = this->P_ - K * S * K.transpose();

	this->NIS_radar_ =  z_diff.transpose() * S * z_diff;

	this->DebugPrint("UpdateRadar: Leave.");
}




void UKF::FirstMeasurement(MeasurementPackage measurement_pack)
{
	// Initialize the state ekf_.x_ with the first measurement.
	// Initialize the covariance matrix.

	this->P_ <<   0.5, 0,   0,  0,  0,
				  0,   0.5, 0,  0,  0,
				  0,   0,  10,  0,  0,
				  0,   0,   0, 10,  0,
				  0,   0,   0,  0, 10;

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
	angle = fmod(angle + M_PI, 2 * M_PI);
	if (angle < 0)
		angle += 2 * M_PI;

	return angle - M_PI;
}


void UKF::DebugPrint(const char* message)
{
	if(DEBUGGING){
		std::cout << "DEBUG: " << message << std::endl;
	}
}

void UKF::DebugPrint(const char* message, const MatrixXd &matrix_to_print)
{
	if(DEBUGGING){
		std::cout << "DEBUG: " << message << std::endl;
		std::cout << matrix_to_print << std::endl;
	}
}
