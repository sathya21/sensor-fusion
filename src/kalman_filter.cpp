#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {
}

MatrixXd CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	float epsilon = 0.0001;
	if (px < epsilon && py < epsilon) {
		px = epsilon;
		py = epsilon;
	}
	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1 * c2);

	cout << c3;

	//check division by zero
	if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj.setZero();
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0, -(py / c1), (px / c1), 0, 0, py
			* (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py
			/ c2;

	return Hj;
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
		MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() {
	/**
	 TODO:
	 * predict the state
	 */

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
	cout << P_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	 TODO:
	 * update the state by using Kalman Filter equations
	 */

	cout << "H_";
	cout << H_;
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	cout << "z_pred";
	cout << z_pred;

	cout << "z";
	cout << z;
	cout << "yx";

	MatrixXd Ht = H_.transpose();

	cout << "Ht";
	cout << "\nHt\n";
	cout << Ht;

	cout << "\nP_\n";
	cout << P_;

	cout << "\nR_\n";
	cout << R_;

	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	cout << "PHt\n";

	MatrixXd K = PHt * Si; //new estimate

	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	 TODO:
	 * update the state by using Extended Kalman Filter equations
	 */
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	float epsilon = 0.0001;
	if (px < epsilon && py < epsilon) {
		px = epsilon;
		py = epsilon;
	}
	float rho = sqrt(pow(px, 2) + pow(py, 2));
	float phi = atan2(py, px);

	if (phi < -3.14) {
		phi = phi + 2 * 3.14;
		cout << "\nphi\n";
		cout << phi;
	} else if (phi > 3.14) {
		phi = phi - 2 * 3.14;
		cout << "\nphi\n";
		cout << phi;
	}

	if (rho == 0) {
		rho = epsilon;
	}

	float rho_dot = (px * vx + py * vy) / rho;

	VectorXd z_pred = VectorXd(3);
	z_pred << rho, phi, rho_dot;
	cout << "\nz\n";
	cout << z;

	MatrixXd Hj = CalculateJacobian(x_);
	cout << "\nhj\n";
	cout << Hj;
	cout << "\nx_\n";
	cout << x_;
	//z_pred = H_ * x_;
	cout << "\nz_pred\n";
	cout << z_pred;

	cout << "\nz\n";
	cout << z;
	VectorXd y = z - z_pred;

	cout << "\nzy\n";

	MatrixXd Ht = Hj.transpose();

	cout << "\nHt\n";
	cout << Ht;

	cout << "\nP_\n";
	cout << P_;

	cout << "\ntestR_\n";
	cout << R_;

	MatrixXd test = Hj * P_ * Ht;
	cout << "\ntest\n";
	cout << test;

	MatrixXd S = Hj * P_ * Ht + R_;

	cout << "\nS\n";
	cout << S;

	MatrixXd Si = S.inverse();
	cout << "\nsi\n";
	cout << Si;

	MatrixXd PHt = P_ * Ht;
	cout << "\npHt\n";
	cout << PHt;

	MatrixXd K = PHt * Si;

	cout << "\nK\n";
	cout << K;

	cout << "\ny\n";

	cout << y;
	x_ = x_ + (K * y);

	cout << "\nx_\n";

	cout << x_;

	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - (K * Hj)) * P_;

}

