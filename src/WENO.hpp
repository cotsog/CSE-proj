/*
 * WENO.h
 *
 *  Created on: Apr 15, 2014
 *      Author: lurker
 */

#include "Util.hpp"

// only work for periodical boundary condition

#ifndef WENO_H_
#define WENO_H_

class WENO {
public:
	WENO();
	virtual ~WENO();
	// base methods
	void WENO1D(int flag, Vector<double> init_state, Grid_1D grid, int Order, int time_step, double delta_t);
	void WENO2D(int flag, Matrix<double> init_state, Grid_2D grid, int Order, int time_step, double delta_t);
	void WENO3D(int flag, Tensor<double> init_state, Grid_3D grid, int Order, int time_step, double delta_t);

	// variable velocity
	void VWENO1D(Vector<double> init_state, Grid_1D grid,  int Order,int time_step, double delta_t, FuncVel_1D Vel_X);
	void VWENO2D(Matrix<double> init_state, Grid_2D grid,  int Order,int time_step, double delta_t, FuncVel_2D Vel_X, FuncVel_2D Vel_Y);
	void VWENO3D(Tensor<double> init_state, Grid_3D grid,  int Order,int time_step, double delta_t, FuncVel_3D Vel_X, FuncVel_3D Vel_Y, FuncVel_3D Vel_Z);

	// constant velocity
	void CWENO1D(Vector<double> init_state, Grid_1D grid,  int Order,int time_step, double delta_t, Const Vel_X);
	void CWENO2D(Matrix<double> init_state, Grid_2D grid,  int Order,int time_step ,double delta_t, Const Vel_X, Const Vel_Y);
	void CWENO3D(Tensor<double> init_state, Grid_3D grid,  int Order,int time_step, double delta_t, Const Vel_X, Const Vel_Y, Const Vel_Z);

	// WENO matrix
	void Assign(Matrix<double>&,Matrix<double>&, int Order);
	void MakeXi(double xi, Vector<double>& Xi, int Order);

	// RK
	double RK_2D(FuncVel_2D Vel_Unknown, int axis, double x, double y ,double time_start, double time_end);

	// Splitting
	void Split_2D(FuncVel_2D Vel_Unknown, int axis, Matrix<double>& Unknown_shift, Matrix<int>& Unknown_rotate,
			Matrix<double>& Unknown_xi, Grid_2D grid, double delta_t, double time_start, double time_end);

	void Update_2D(int axis, Matrix<double>& Solution, Matrix<double>& Aux, Matrix<double> xi, Matrix<int> rotate, Vector<double>& Xi, int Order,Matrix<double> CL, Matrix<double> CR);
};

#endif /* WENO_H_ */
