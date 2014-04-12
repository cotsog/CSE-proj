/*
 * Interpolation.h
 *
 *  Created on: Apr 11, 2014
 *      Author: lurker
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#define SPLINE_TYPE int

#include "Util.h"

class Interpolation {
public:
	Interpolation();
	virtual ~Interpolation();
	// will use Eigen/unsupported/Spline for interpolation
	// divergence free interpolation is also required.
};

#endif /* INTERPOLATION_H_ */
