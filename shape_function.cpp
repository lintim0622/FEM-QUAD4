#pragma once

#include "shape_function.h"

Shape_Function::Shape_Function() : N1{}, N2{}, N3{}, N4{} {

	N.setZero();
	// std::cout << "Shape_Function\n";

};

Shape_Function::~Shape_Function() {};

void Shape_Function::local_N(const double& xi, const double& eta) {

	N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
	N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
	N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
	N4 = 0.25 * (1.0 - xi) * (1.0 + eta);

}

void Shape_Function::matrix() {

	N <<  N1, 0.0,  N2, 0.0,  N3, 0.0,  N4, 0.0,
		 0.0,  N1, 0.0,  N2, 0.0,  N3, 0.0,  N4;

}

Strain_Displacement::Strain_Displacement() {

	dN1xi  = 0.0;
	dN1eta = 0.0;
	dN2xi  = 0.0;
	dN2eta = 0.0;
	dN3xi  = 0.0;
	dN3eta = 0.0;
	dN4xi  = 0.0;
	dN4eta = 0.0;

	J11 = 0.0;
	J12 = 0.0;
	J21 = 0.0;
	J22 = 0.0;
	det_J = 0.0;

	B.setZero();

};

Strain_Displacement::~Strain_Displacement() {};

void Strain_Displacement::local_dN(const double& xi, const double& eta) {

	dN1xi  = 0.25 * (-1.0) * (1.0 - eta);
	dN1eta = 0.25 * (1.0 - xi) * (-1.0);

	dN2xi  = 0.25 * (1.0) * (1.0 - eta);
	dN2eta = 0.25 * (1.0 + xi) * (-1.0);

	dN3xi  = 0.25 * (1.0) * (1.0 + eta);
	dN3eta = 0.25 * (1.0 + xi) * (1.0);

	dN4xi  = 0.25 * (-1.0) * (1.0 + eta);
	dN4eta = 0.25 * (1.0 - xi) * (1.0);

}

void Strain_Displacement::Jacobian(double x[], double y[]) {

	double& x1 = x[0];
	double& x2 = x[1];
	double& x3 = x[2];
	double& x4 = x[3];
	double& y1 = y[0];
	double& y2 = y[1];
	double& y3 = y[2];
	double& y4 = y[3];

	J11 = x1 * dN1xi  + x2 * dN2xi  + x3 * dN3xi  + x4 * dN4xi;
	J12 = y1 * dN1xi  + y2 * dN2xi  + y3 * dN3xi  + y4 * dN4xi;
	J21 = x1 * dN1eta + x2 * dN2eta + x3 * dN3eta + x4 * dN4eta;
	J22 = y1 * dN1eta + y2 * dN2eta + y3 * dN3eta + y4 * dN4eta;

	det_J = J11 * J22 - J12 * J21;

}

void Strain_Displacement::matrix() {

	double det_J_inv = 1.0 / det_J;

	double dN1x = det_J_inv * ( J22 * dN1xi - J12 * dN1eta );
	double dN1y = det_J_inv * (-J21 * dN1xi + J11 * dN1eta );

	double dN2x = det_J_inv * ( J22 * dN2xi - J12 * dN2eta );
	double dN2y = det_J_inv * (-J21 * dN2xi + J11 * dN2eta );

	double dN3x = det_J_inv * ( J22 * dN3xi - J12 * dN3eta );
	double dN3y = det_J_inv * (-J21 * dN3xi + J11 * dN3eta );

	double dN4x = det_J_inv * ( J22 * dN4xi - J12 * dN4eta );
	double dN4y = det_J_inv * (-J21 * dN4xi + J11 * dN4eta );

	B << dN1x,  0.0, dN2x,  0.0, dN3x,  0.0, dN4x,  0.0,
		  0.0, dN1y,  0.0, dN2y,  0.0, dN3y,  0.0, dN4y,
		 dN1y, dN1x, dN2y, dN2x, dN3y, dN3x, dN4y, dN4x;

}	