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

	// RK 4th order, will produce 3rd order accuracy in locating position.
	// the highest accuracy will be 3 at most.
	double RK_2D(FuncVel_2D Vel_Unknown, int axis, double x, double y ,double time_start, double time_end);
	double RK_3D(FuncVel_3D Vel_Unknown, int axis, double x, double y, double z, double time_start, double time_end);

	// Splitting
	void Split_2D(FuncVel_2D Vel_Unknown, int axis, Matrix<double>& Unknown_shift, Matrix<int>& Unknown_rotate,
			Matrix<double>& Unknown_xi, Grid_2D grid, double delta_t, double time_start, double time_end);
	// update
	void Update_2D(int axis, Matrix<double>& Solution, Matrix<double>& Aux, Matrix<double>& xi, Matrix<int>& rotate, Vector<double>& Xi,
			int Order,Matrix<double>& CL, Matrix<double>& CR);
	// Split scheme, can be updated. Already have first order, Strang, Yoshida
	Matrix<double> SplitScheme(int order, int dimenstion);
	
	// Splitting 
	void Split_3D(FuncVel_3D _Vel, int _axis, Tensor<double>& _shift, Tensor<int>& _rotate, 
			Tensor<double>& _xi, Grid_3D grid, double delta_t, double time_start, double time_end);
	// update	
	void Update_3D(int _axis, Tensor<double>& Solution, Tensor<double>& Aux, Tensor<double> _xi, Tensor<int> _rotate, Vector<double>& Xi,
			int Order, Matrix<double> CL, Matrix<double> CR);
	long timer()
	{
	    struct timeval tv;
	    gettimeofday(&tv, NULL);
	    return (tv.tv_sec * 1000000) + tv.tv_usec;
	};

};

#endif /* WENO_H_ */
