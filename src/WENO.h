/*
 * WENO.h
 *
 *  Created on: Apr 10, 2014
 *      Author: lurker
 */

#ifndef WENO_H_
#define WENO_H_

#include "Util.h"

class WENO {
public:
	WENO();
	virtual ~WENO();
	// WENO SCHEME
	// PERIODICAL BOUNDARY CONDITION
	// CONST VECTOR FIELD
	void PWENO_1D(int flag, int ORDER, const VECTOR init_state, Grid_1D grid, int time_step, double delta_t);
	void PWENO_2D(int flag, int ORDER, const MATRIX init_state, Grid_2D grid, int time_step, double delta_t);
	void PWENO_3D(int flag, int ORDER, const TENSOR init_state, Grid_3D grid, int time_step, double delta_t);
	// VARIABLE VECTOR FIELD
	void VPWENO_1D(int flag, int ORDER, const VECTOR init_state, Grid_1D grid, int time_step, double delta_t, VECTOR& Solution, VECTOR& AUX);
	void VPWENO_2D(int flag, int ORDER, const MATRIX init_state, Grid_2D grid, int time_step, double delta_t, MATRIX& Solution, MATRIX& AUX);
	void VPWENO_3D(int flag, int ORDER, const TENSOR init_state, Grid_3D grid, int time_step, double delta_t, TENSOR& Solution, TENSOR& AUX);

	//ABC BOUNDARY CONDITION, PML needed.
	void NPWENO_1D(int flag,int ORDER, const VECTOR init_state, Grid_1D grid, int time_step, double delta_t);
	void NPWENO_2D(int flag, int ORDER, const MATRIX init_state, Grid_2D grid, int time_step, double delta_t);
	void NPWENO_3D(int flag, int ORDER, const TENSOR init_state, Grid_3D grid, int time_step, double delta_t);
	// UTILITIES for VECTOR/MATRIX STUFF
	void ASSEMBLE(MATRIX &C, int ORDER, double xshift);
	void EXTEND(double& , int&, int& ,int, int&);
	void MAKEXI(double, VECTOR&, int);

	// UTILITIES for quick implementation
	bool ph(double);
	bool nh(double);

};

#endif /* WENO_H_ */
