#pragma once

#include "mesh.h"

class Shape_Function {
public:

	// function
	Shape_Function();
	~Shape_Function();

	void local_N(const double& xi, const double& eta);
	void matrix();

	// member
	double N1;
	double N2;
	double N3;
	double N4;
	Eigen::Matrix<double, 2, 8> N;
};

class Strain_Displacement {
public:

	// function
	Strain_Displacement();
	~Strain_Displacement();

	void local_dN(const double& xi, const double& eta);
	void Jacobian(double x[], double y[]);
	void matrix();

	// member
	double dN1xi;
	double dN1eta;
	double dN2xi;
	double dN2eta;
	double dN3xi;
	double dN3eta;
	double dN4xi;
	double dN4eta;

	double J11;
	double J12;
	double J21;
	double J22;
	double det_J;

	Eigen::Matrix<double, 3, 8> B;
};
