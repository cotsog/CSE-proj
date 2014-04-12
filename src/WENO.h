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
	void PWENO_1D(int flag, int ORDER, const VECTOR init_state, Grid_1D grid, int time_step, double delta_t);
	void PWENO_2D(int flag, int ORDER, const MATRIX init_state, Grid_2D grid, int time_step, double delta_t);
	void PWENO_3D(int flag, int ORDER, const TENSOR init_state, Grid_3D grid, int time_step, double delta_t);
	void NPWENO_1D(int flag,int ORDER, const VECTOR init_state, Grid_1D grid, int time_step, double delta_t);
	void NPWENO_2D(int flag, int ORDER, const MATRIX init_state, Grid_2D grid, int time_step, double delta_t);
	void NPWENO_3D(int flag, int ORDER, const TENSOR init_state, Grid_3D grid, int time_step, double delta_t);
	// UTILITIES for MATRIX STUFF
	void ASSEMBLE(MATRIX &C, int ORDER, double xshift);

	// UTILITIES for quick implementation
	bool ph(double);
	bool nh(double);

};

#endif /* WENO_H_ */
