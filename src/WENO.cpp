/*
 * WENO.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: lurker
 */

#include "WENO.h"

WENO::WENO() {
	// TODO Auto-generated constructor stub

}

WENO::~WENO() {
	// TODO Auto-generated destructor stub
}

void WENO::PWENO_1D(int flag, int ORDER, VECTOR init_state, Grid_1D grid, int time_step, double delta_t){
	// vector is row-majored.
	std::size_t length = init_state.rows();

	VECTOR Solution = VECTOR::Zero(length + ORDER);
	VECTOR AUX      = VECTOR::Zero(length + ORDER);

	if (flag == 0){
		// constant velocity
#ifndef VEL_1D
#define VEL_1D 1.0
#endif
		int LEFT,RIGHT, ROTATE;
		double xshift = VEL_1D*delta_t/(grid.final_x- grid.start_x)*grid.size_x;
		// DIRECTION OF PROPAGATION
		EXTEND(xshift, LEFT, RIGHT, ORDER, ROTATE);
		MATRIX C(ORDER,ORDER);
		VECTOR XI(ORDER);
		VECTOR SUB_1(ORDER);
		VECTOR SUB_2(ORDER);

		MAKEXI(xshift, XI, ORDER);

		// Assembly C according to x_shift
		// Now finished with x_shift in (-.5,.5)
		ASSEMBLE(C,ORDER,xshift);
		// set solution domain, the other points are zero
		Solution.segment(LEFT,length) = init_state;
		// updates time_step times
		for (int i = 0; i < time_step; i++){
			// periodical boundary condition, difference in index as (length - 1)
			// complete Solution as full solution
			Solution.segment(0,LEFT) = Solution.segment(length - 1,LEFT);
			Solution.segment(length + LEFT,RIGHT) = Solution.segment(LEFT + 1,RIGHT);
			// end of boundary condition
			// updating according to time propagating
			for (std::size_t j = LEFT; j < length+LEFT; j ++){
				// extract information from full Solution
				SUB_1= Solution.segment(j-LEFT + 1,ORDER);
				SUB_2= Solution.segment(j-LEFT,ORDER);
				// x_shift > 0, and x_shift < 0
				// up-wind-similar scheme
				// update AUX in solution domain
				AUX(j) = Solution(j) - (xshift)*(SUB_1-SUB_2).transpose()*C*XI;
			}
			Solution.segment(LEFT,length) = AUX.segment(LEFT,length);
			// During this time step, if ROTATE is non-zero, then a rotation of solution is needed.
			// Direction of rotation
			if (ROTATE > 0){
				// MOVE TO RIGHT SIDE
				AUX.segment(LEFT+ROTATE, length - ROTATE) = Solution.segment(LEFT, length - ROTATE);
				// periodical
				AUX.segment(LEFT,ROTATE) = Solution.segment(LEFT + length - ROTATE - 1, ROTATE);
				Solution.segment(LEFT,length) = AUX.segment(LEFT,length);
			}
			else if (ROTATE < 0){
				// MOVE TO LEFT SIDE, ROTATE IS NEGATIVE
				AUX.segment(LEFT,length + ROTATE) = Solution.segment(LEFT -  ROTATE, length + ROTATE);
				AUX.segment(length + LEFT + ROTATE, -ROTATE) = Solution.segment(LEFT+1, -ROTATE);
				Solution.segment(LEFT, length) = AUX.segment(LEFT,length);
			}
		}
		// write into file or not
		// Here is test output. L2 norm
		std::cout << sqrt((Solution.segment(LEFT,length) - init_state).squaredNorm()/init_state.rows());
	}
	else if (flag == 1){
		//variable velocity
		std::cout << "not implemented.\n" << std::endl;
	}
}

void WENO::PWENO_2D(int flag, int ORDER, MATRIX init_state, Grid_2D grid, int time_step, double delta_t){
	// use two component splitting method on each time step
	// consider x-direction first then y_direction. First order splitting would be enough.
	std::size_t len_row = init_state.rows();
	std::size_t len_col = init_state.cols();
	MATRIX Solution = MATRIX::Zero(len_row + ORDER, len_col + ORDER);
	MATRIX AUX      = MATRIX::Zero(len_row + ORDER, len_col + ORDER);
	if (flag == 0){
#ifndef VEL_2D_X
#define VEL_2D_X 1.0
#endif
#ifndef VEL_2D_Y
#define VEL_2D_Y 1.0
#endif
		int LEFT_X,RIGHT_X,LEFT_Y,RIGHT_Y, ROTATE_X, ROTATE_Y;
		double xshift = VEL_2D_X*delta_t/(grid.final_x- grid.start_x)*grid.size_x;
		double yshift = VEL_2D_Y*delta_t/(grid.final_y - grid.start_y)*grid.size_y;
		// SOLUTION DOMAIN TO VIRTUAL SOLUTION DOMAIN
		EXTEND(xshift, LEFT_X, RIGHT_X,ORDER, ROTATE_X);
		EXTEND(yshift, LEFT_Y, RIGHT_Y, ORDER, ROTATE_Y);

		// EXTENDED DOMAIN
		MATRIX C_X(ORDER,ORDER);
		MATRIX C_Y(ORDER,ORDER);
		VECTOR XI_X(ORDER);
		VECTOR XI_Y(ORDER);
		VECTOR SUB_1(ORDER);
		VECTOR SUB_2(ORDER);

		MAKEXI(xshift, XI_X, ORDER);
		MAKEXI(yshift, XI_Y, ORDER);

		ASSEMBLE(C_X,ORDER,xshift);
		ASSEMBLE(C_Y,ORDER,yshift);
		// Initialize solution
		Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = init_state;
		// update Solution
		for (int i = 0; i < time_step; i++){
			// BLOCK-WISE EXTENSION, ROW_BASE
			Solution.block(0, LEFT_Y, LEFT_X, len_col) = Solution.block(len_row - 1, LEFT_Y, LEFT_X, len_col);
			Solution.block(len_row + LEFT_X, LEFT_Y, RIGHT_X, len_col) = Solution.block(LEFT_X + 1, LEFT_Y, RIGHT_X, len_col);
//			Solution.block(LEFT_X, 0, len_row, LEFT_Y) = Solution.block(LEFT_X, len_col - 1, len_row, LEFT_Y);
//			Solution.block(LEFT_X, len_col + LEFT_Y, len_row, RIGHT_Y) = Solution.block(LEFT_X, LEFT_Y + 1, len_row, RIGHT_Y);
			// TWO STEP SPLITTING METHOD.
			// X-DIRECTION
			// TODO: USE block－wise operation
			for (std::size_t j_x = LEFT_X; j_x < len_row+LEFT_X; j_x++){
				for (std::size_t j_y = LEFT_Y; j_y < len_col + LEFT_Y; j_y++){
					SUB_1= Solution.block(j_x-LEFT_X + 1,j_y, ORDER,1);
					SUB_2= Solution.block(j_x-LEFT_X    ,j_y, ORDER,1);
					AUX(j_x,j_y) = Solution(j_x,j_y) - (xshift)*(SUB_1-SUB_2).transpose()*C_X*XI_X;
				}
			}
			Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);
			// X-DIRECTION COMPLETED

			// EXTEND SOLUTION BLOCK-WISE, COL_BASE
//			Solution.block(0, LEFT_Y, LEFT_X, len_col) = Solution.block(len_row - 1, LEFT_Y, LEFT_X, len_col);
//			Solution.block(len_row + LEFT_X, LEFT_Y, RIGHT_X, len_col) = Solution.block(LEFT_X + 1, LEFT_Y, RIGHT_X, len_col);
			Solution.block(LEFT_X, 0, len_row, LEFT_Y) = Solution.block(LEFT_X, len_col - 1, len_row, LEFT_Y);
			Solution.block(LEFT_X, len_col + LEFT_Y, len_row, RIGHT_Y) = Solution.block(LEFT_X, LEFT_Y + 1, len_row, RIGHT_Y);
			// Y-DIRECTION
			// TODO: USE block-wise operation
			for (std::size_t j_x = LEFT_X; j_x < len_row+LEFT_X; j_x++){
				for (std::size_t j_y = LEFT_Y; j_y < len_col + LEFT_Y; j_y++){
					SUB_1= Solution.block(j_x, j_y-LEFT_Y + 1,1, ORDER).transpose();
					SUB_2= Solution.block(j_x, j_y-LEFT_Y,1, ORDER).transpose();
					AUX(j_x,j_y) = Solution(j_x,j_y) - (yshift)*(SUB_1-SUB_2).transpose()*C_Y*XI_Y;
				}
			}
			// update step. should be accelerated.
			Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);

			if (ROTATE_X > 0){
				// MOVE TO RIGHT SIDE
				AUX.block(LEFT_X+ROTATE_X,LEFT_Y, len_row - ROTATE_X, len_col) = Solution.block(LEFT_X, LEFT_Y, len_row - ROTATE_X, len_col);
				// periodical
				AUX.block(LEFT_X,LEFT_Y,ROTATE_X,len_col) = Solution.block(LEFT_X + len_row - ROTATE_X - 1, LEFT_Y,ROTATE_X, len_col);
				Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);
			}
			else if (ROTATE_X < 0){
				// MOVE TO LEFT SIDE, ROTATE IS NEGATIVE
				AUX.block(LEFT_X,LEFT_Y,len_row + ROTATE_X, len_col) = Solution.block(LEFT_X -  ROTATE_X, LEFT_Y, len_row + ROTATE_X, len_col);
				AUX.block(len_row + LEFT_X + ROTATE_X,LEFT_Y, -ROTATE_X, len_col) = Solution.block(LEFT_X+1,LEFT_Y, -ROTATE_X, len_col);
				Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);
			}

			if (ROTATE_Y > 0){
				// MOVE TO RIGHT SIDE
				AUX.block(LEFT_X, LEFT_Y+ROTATE_Y,len_row, len_col - ROTATE_Y) = Solution.block(LEFT_X, LEFT_Y, len_row, len_col - ROTATE_Y);
				// periodical
				AUX.block(LEFT_X, LEFT_Y,len_row, ROTATE_Y) = Solution.block(LEFT_X, LEFT_Y + len_col - ROTATE_Y - 1, len_row, ROTATE_Y);
				Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);
			}
			else if (ROTATE_Y < 0){
				// MOVE TO LEFT SIDE, ROTATE IS NEGATIVE
				AUX.block(LEFT_X, LEFT_Y,len_row, len_col + ROTATE_Y) = Solution.block(LEFT_X, LEFT_Y -  ROTATE_Y, len_row, len_col + ROTATE_Y);
				AUX.block(LEFT_X, len_col + LEFT_Y + ROTATE_Y, len_row,  -ROTATE_Y) = Solution.block(LEFT_X, LEFT_Y+1,len_row, -ROTATE_Y);
				Solution.block(LEFT_X,LEFT_Y,len_row, len_col) = AUX.block(LEFT_X,LEFT_Y,len_row, len_col);
			}
		}
		std::cout << sqrt((Solution.block(LEFT_X,LEFT_Y,len_row, len_col)  - init_state).squaredNorm()/init_state.rows()/init_state.cols());
	}
	else if (flag == 1){
		std::cout << "not implemented yet.\n" << std::endl;
	}
}
//void WENO::PWENO_3D(int flag, int ORDER, TENSOR init_state, Grid_3D grid, int time_step, double delta_t){
//
//}
//void WENO::NPWENO_1D(int flag, int ORDER, VECTOR init_state, Grid_1D grid, int time_step, double delta_t){
//	// PML(ABC) needed.
//}
//void WENO::NPWENO_2D(int flag,int ORDER,  MATRIX init_state, Grid_2D grid, int time_step, double delta_t){
//
//}
//void WENO::NPWENO_3D(int flag, int ORDER,  TENSOR init_state, Grid_3D grid, int time_step, double delta_t){
//
//}

void WENO::ASSEMBLE(MATRIX &C, int ORDER, double xshift){
	if (ORDER == 3){
		if (ph(xshift)){
			//CL
			C << -1./6,     0, 1./6,
				   5./6, 1./2, -1./3,
				   1./3, -1./2, 1./6;
		}
		else if(nh(xshift)){
			//CR
			C << 1./3, -1./2, 1./6,
				  5./6, 1./2, -1./3,
				 -1./6,     0, 1./6;
		}

	}
	if (ORDER == 5){
		if (ph(xshift)){
		//CL
		C << 1./30,     0., -1./24,      0.,  1./120,
			-13./60, -1./24,   1./4,   1./24,  -1./30,
			47./60,   5./8,  -1./3,   -1./8,   1./20,
			9./20,  -5./8,  1./12,    1./8,  -1./30,
			-1./20,  1./24,  1./24,  -1./24,  1./120;
		}
		else if (nh(xshift)){
			C << -1./20,   1./24,  1./24,  -1./24,  1./120,
				 9./20,  -5./8,  1./12,    1./8,  -1./30,
				47./60,   5./8,  -1./3,   -1./8,   1./20,
			   -13./60, -1./24,   1./4,   1./24,  -1./30,
				 1./30,     0., -1./24,      0.,  1./120;
		}
	}
	if (ORDER == 7){
		if (ph(xshift)){
		//CL
		C << -1./140, 0., 7./720, 0, -1./360, 0., 1./5040,
				5./84, 1./180, -19./240, -1./144, 1./48, 1./720, -1./840,
				-101./420, -5./72, 7./24, 11./144, -13./240, -1./144, 1./336,
				319./420, 49./72, -23./72, -7./36, 23./360, 1./72, -1./252,
				107./210, -49./72, 1./48, 7./36, -1./30, -1./72, 1./336,
				-19./210, 5./72, 7./80, -11./144, 1./240, 1./144, -1./840,
				1./105, -1./180, -1./90, 1./144, 1./720, -1./720, 1./5040;
		}
		else if (nh(xshift)){
			C << 1./105, -1./180, -1./90, 1./144, 1./720, -1./720, 1./5040,
					-19./210, 5./72, 7./80, -11./144, 1./240, 1./144, -1./840,
					107./210, -49./72, 1./48, 7./36, -1./30, -1./72, 1./336,
					319./420, 49./72, -23./72, -7./36, 23./360, 1./72, -1./252,
					-101./420, -5./72, 7./24, 11./144, -13./240, -1./144, 1./336,
					5./84, 1./180, -19./240, -1./144, 1./48, 1./720, -1./840,
					-1./140, 0., 7./720, 0, -1./360, 0., 1./5040;
		}
	}

	if (ORDER == 9){
		if (ph(xshift)){
			//CL
			C <<        1./630,         0.,  -41./18144,        0., 13./17280,        0., -1./12096,        0.,  1./362880,
					 -41./2520,   -1./1120, 2081./90720,   7./5760,   -1./135,  -1./2880, 23./30240,  1./40320,  -1./45360,
					 199./2520,   17./1440,  -281./2592, -89./5760, 139./4320,  11./2880, -17./6048,  -1./5760,   1./12960,
					-641./2520, -127./1440, 4097./12960, 587./5760,  -29./432, -41./2880,167./30240,   1./1920,   -1./6480,
					1879./2520,   205./288,  -797./2592,  -91./384, 587./8640,    5./192, -19./3024,  -1./1152,    1./5184,
					  275./504,  -205./288,   -59./2592,   91./384, -29./1080,   -5./192,  25./6048,   1./1152,   -1./6480,
					  -61./504,  127./1440, 1637./12960,-587./5760, -17./4320,  41./2880,-43./30240,  -1./1920,   1./12960,
					   11./504,  -17./1440, -491./18144,  89./5760,  11./2160, -11./2880,   1./6048,   1./5760,  -1./45360,
					   -1./504,    1./1120,   59./22680,  -7./5760,-11./17280,   1./2880,  1./60480, -1./40320,  1./362880;
		}
		else if (nh(xshift)){
			C << -1./504,    1./1120,   59./22680,  -7./5760,-11./17280,   1./2880,  1./60480, -1./40320,  1./362880,
					 11./504,  -17./1440, -491./18144,  89./5760,  11./2160, -11./2880,   1./6048,   1./5760,  -1./45360,
					 -61./504,  127./1440, 1637./12960,-587./5760, -17./4320,  41./2880,-43./30240,  -1./1920,   1./12960,
					 275./504,  -205./288,   -59./2592,   91./384, -29./1080,   -5./192,  25./6048,   1./1152,   -1./6480,
					 1879./2520,   205./288,  -797./2592,  -91./384, 587./8640,    5./192, -19./3024,  -1./1152,    1./5184,
					 -641./2520, -127./1440, 4097./12960, 587./5760,  -29./432, -41./2880,167./30240,   1./1920,   -1./6480,
					 199./2520,   17./1440,  -281./2592, -89./5760, 139./4320,  11./2880, -17./6048,  -1./5760,   1./12960,
					 -41./2520,   -1./1120, 2081./90720,   7./5760,   -1./135,  -1./2880, 23./30240,  1./40320,  -1./45360,
					    1./630,         0.,  -41./18144,        0., 13./17280,        0., -1./12096,        0.,  1./362880;
		}
	}
}

bool WENO::ph(double xshift){
	if ((xshift >= 0) && (xshift <= 0.5))
		return true;
	else
		return false;
}

bool WENO::nh(double xshift){
	if ((xshift <= 0) && (xshift >= -.5))
		return true;
	else
		return false;
}

void WENO::EXTEND(double& _shift, int& _LEFT, int& _RIGHT, int _ORDER, int& _ROTATE){
	// _ROTATE will give the solution a rotation
	_ROTATE = floor(_shift + 0.5);
	_shift = _shift - _ROTATE;
	if (ph(_shift)){
				_LEFT = (_ORDER + 1)/2;
				_RIGHT = (_ORDER - 1)/2;
			}
	else  if(nh(_shift)){
		_LEFT = (_ORDER - 1)/2;
		_RIGHT = (_ORDER + 1)/2;
	}
	else{
		std::cout << "exception threw out." << std::endl;
		exit(1);
	}
}

void WENO::MAKEXI(double _shift, VECTOR& _XI, int ORDER){
	_XI(0) = 1.0;
	for (int i = 1; i < ORDER; i++){
		_XI(i) = fabs(_shift)*_XI(i-1);
	}
}
