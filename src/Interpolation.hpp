/*
 * Interpolation.h
 *
 *  Created on: Apr 15, 2014
 *      Author: lurker
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "Util.hpp"

class Interpolation {
public:
	Interpolation();
	virtual ~Interpolation();

	double Linear1D(Vector<double> data, Grid_1D grid, double x);
	double Linear2D(Matrix<double> data, Grid_2D grid, double x, double y);
	double Linear3D(Tensor<double> data, Grid_3D grid, double x, double y, double z);

	double CbSpline1D(Vector<double> data, Grid_1D grid, double x);
	double CbSpline2D(Matrix<double> data, Grid_2D grid, double x, double y);
	double CbSpline3D(Tensor<double> data, Grid_3D grid, double x, double y, double z);
};

#endif /* INTERPOLATION_H_ */
