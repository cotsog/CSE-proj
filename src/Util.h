/*
 * Util.h
 *
 *  Created on: Apr 10, 2014
 *      Author: lurker
 */

#ifndef UTIL_H_
#define UTIL_H_
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <vector>
// for 3D
//#include "boost/multi_array.hpp"
#include <valarray>
// for 1D and 2D
#include "Eigen/Dense"

#define GO_IN std::cout<<"go in"<<std::endl;
#define GO_OUT std::cout<<"go out"<<std::endl;

#define ORDER3 3
#define ORDER5 5
#define ORDER7 7
#define ORDER9 9

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#define VECTOR Eigen::VectorXd
#define MATRIX Eigen::MatrixXd
// used for 3D, experimental use
// can also try boost::multi_array, vector, but not optimized.
#define TENSOR std::valarray<std::valarray<std::valarray<double>>>

class Grid_1D{
public:
	Grid_1D(double a,double b,int& s){
		start_x = a;
		final_x = b;
		size_x  = s;
	};
	virtual ~Grid_1D(){
	};
	double start_x;
	double final_x;
	size_t size_x;
};

class Grid_2D{
public:
	Grid_2D(double a, double b, int&s, double c, double d, int&t){
		start_x = a;start_y = c;
		final_x = b;final_y = d;
		size_x  = s;size_y  = t;

	}
	virtual ~Grid_2D(){
	}
	double start_x;
	double final_x;
	double start_y;
	double final_y;
	size_t size_x;
	size_t size_y;
};

class Grid_3D{
public:
	Grid_3D(double a, double b ,int& s, double c, double d, int& t,
			double e, double f, int& r){
		start_x = a;start_y = c;start_z = e;
		final_x = b;final_y = d;final_z = f;
		size_x  = s;size_y  = t;size_z = r;
	}
	virtual ~Grid_3D(){

	}
	double start_x;
	double final_x;
	double start_y;
	double final_y;
	double start_z;
	double final_z;
	size_t size_x;
	size_t size_y;
	size_t size_z;
};
// 1D constant
#define VEL_1D 1.0
// 2D constant
#define VEL_2D_X 1.0
#define VEL_2D_Y 1.0
// 3D constant
#define VEL_3D_X 1.0
#define VEL_3D_Y 1.0
#define VEL_3D_Z 1.0

#define PRINT(X) std::cout<< X << std::endl;


#endif /* UTIL_H_ */
