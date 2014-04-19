//============================================================================
// Name        : WENO.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Util.hpp"
#include "WENO.h"
using namespace std;

int test1D() {
	std::cout<< std::setw(8) << "N" << std::setw(10) << "steps" << std::setw(15)
		<< "delta_t" << std::setw(20) << "ORDER3" << std::setw(12) << "ORDER5" <<
		std::setw(15) << "ORDER7" << std::setw(15)<< "ORDER9" << '\n';
		int N = 5;
		int time_step = 3;
		double delta_t = 1./3;
		for (int k = 1;k < 6; k++){
			N *=2;
//			time_step *= 2;
//			delta_t /= 2.;
			Grid_1D grid(.0,1.0, N);
			Vector<double> init_state(N);
			for (int i = 0 ; i < N ; i++){
				init_state(i) = sin(2*PI*i/(N));
			}
			// periodic boundary
			WENO SOLVER;
			std::cout<< std::setw(8) << N << std::setw(10) << time_step <<
					std::setw(15) << delta_t << std::setw(20);
			SOLVER.CWENO1D(init_state,grid,3,time_step,delta_t,1.0);
			SOLVER.CWENO1D(init_state,grid,5,time_step,delta_t,1.0);
			SOLVER.CWENO1D(init_state,grid,7,time_step,delta_t,1.0);
			SOLVER.CWENO1D(init_state,grid,9,time_step,delta_t,1.0);
			std::cout << '\n';
		}

		return 0;
}

int test2D() {
	std::cout<< std::setw(8) << "N" << std::setw(10) << "steps" << std::setw(15)
		<< "delta_t" << std::setw(20) << "ORDER3" << std::setw(12) << "ORDER5" <<
		std::setw(15) << "ORDER7" << std::setw(15)<< "ORDER9" << '\n';
		int N = 5;
		int time_step = 3;
		double delta_t = 1./3;
		for (int k = 1;k < 6; k++){
			N *=2;
//			time_step *= 2;
//			delta_t /= 2.;
			Grid_2D grid(.0,1.0, N,0.,1.,N);
			Matrix<double> init_state(N,N);
			for (int i = 0 ; i < N ; i++){
				for(int j=0; j<N;j++){
					init_state(i,j) = sin(2*PI*i/(N))*sin(2*PI*j/(N));
				}

			}
			// periodic boundary
			WENO SOLVER;
			std::cout<< std::setw(8) << N << std::setw(10) << time_step <<
					std::setw(15) << delta_t << std::setw(20);
			SOLVER.CWENO2D(init_state,grid,3,time_step,delta_t,1.0,1.0);
			SOLVER.CWENO2D(init_state,grid,5,time_step,delta_t,1.0,1.0);
			SOLVER.CWENO2D(init_state,grid,7,time_step,delta_t,1.0,1.0);
			SOLVER.CWENO2D(init_state,grid,9,time_step,delta_t,1.0,1.0);
			std::cout << '\n';
		}

		return 0;
}

int test3D() {
	std::cout<< std::setw(8) << "N" << std::setw(10) << "steps" << std::setw(15)
		<< "delta_t" << std::setw(20) << "ORDER3" << std::setw(12) << "ORDER5" <<
		std::setw(15) << "ORDER7" << std::setw(15)<< "ORDER9" << '\n';
		int N = 5;
		int time_step = 3;
		double delta_t = 1./3;
		for (int k = 1;k < 5; k++){
			N *=2;
//			time_step *= 2;
//			delta_t /= 2.;
			Grid_3D grid(.0,1.0, N,0.,1.,N , 0., 1., N);
			Tensor<double> init_state(N,N, N);
			for (int i = 0 ; i < N ; i++){
				for(int j=0; j<N;j++){
					for (int k = 0; k<N; k++){
						init_state(i,j,k) = sin(2*PI*i/(N))*sin(2*PI*j/(N))*sin(2*PI*k/(N));
					}
				}

			}
			// periodic boundary
			WENO SOLVER;
			std::cout<< std::setw(8) << N << std::setw(10) << time_step <<
					std::setw(15) << delta_t << std::setw(20);
			SOLVER.CWENO3D(init_state,grid,3,time_step,delta_t,1.0,1.0,1.0);
			SOLVER.CWENO3D(init_state,grid,5,time_step,delta_t,1.0,1.0,1.0);
			SOLVER.CWENO3D(init_state,grid,7,time_step,delta_t,1.0,1.0,1.0);
			SOLVER.CWENO3D(init_state,grid,9,time_step,delta_t,1.0,1.0,1.0);
			std::cout << '\n';
		}

		return 0;
}

int testTensor(){
	Tensor<int> A(3,3,3);
	Tensor<int> B(3,3,3);
	Vector<int> C(3);
	A(0,0,0) = 1;
	A(1,1,1) = 2;

	A += 1;

	cout << A(0,0,0) << A(1,1,1);

	B += 1;

	A += B;

	cout << A(0,0,0) << A(1,1,1);

	C = A.slicerc(1,1);

	cout << C(0) << C(1) << C(2) ;

	cout << A.sqnorm() << endl;
	return 0;
}
int main(){
	test1D();
	test2D();
	test3D();



}
