/*
 * Util.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: lurker
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_
/* MATH */
#include <cmath>
#include <cstddef>
#include <valarray>
#include <algorithm>
/* IO */
#include <iostream>
#include <fstream>
#include <iomanip>
/* DEBUG SYSTEM */
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
/* TYPES */
#include "Grid.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Tensor.hpp"
/* OPENMP */
#include <omp.h>

using namespace std;

#define CHK cout<<"CHECK"<<endl;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

typedef double (*FuncVel_1D)(double, double);
typedef double (*FuncVel_2D)(double, double ,double);
typedef double (*FuncVel_3D)(double, double ,double, double);

typedef const double Const;


#endif /* UTIL_HPP_ */
