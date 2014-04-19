/*
 * WENO.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: lurker
 */

#include "WENO.h"

WENO::WENO() {
	// TODO Auto-generated constructor stub

}

WENO::~WENO() {
	// TODO Auto-generated destructor stub
}


void WENO::CWENO1D(Vector<double> init_state, Grid_1D grid, int Order, int time_step, double delta_t, Const Vel_X){
	// periodical boundary, the end point on RHS is useless
	Vector<double> Solution(grid.size_x);
	Vector<double> Aux(grid.size_x);

	double xshift = Vel_X * delta_t / grid.delta_x;
	// backward length

	int rotate = floor(xshift + 0.5);
	double xi = xshift - rotate;

	Vector<double> Xi(Order);

	MakeXi(xi, Xi,Order);

	Matrix<double> CL(Order,Order);
	Matrix<double> CR(Order,Order);
	Assign(CL,CR,Order);

	Solution = init_state;

	int i,j;

	if (xi > 0){
		for (i = 0 ; i < time_step; i++){
#pragma omp parallel for private(j) schedule(static)
			for (j = 0; j < grid.size_x; j++){
				Aux(j+rotate) = Solution(j) - xi*((Solution(j-(Order-1)/2, Order) - Solution(j-(Order+1)/2, Order))*(CL*Xi));
			}
//#pragma omp barrier
			Solution = Aux;
		}
	}
	else{
		for (i =0; i < time_step; i++){
#pragma omp parallel for private(j) schedule(static)
			for (j = 0 ; j < grid.size_x; j++){
				Aux(j+rotate) = Solution(j)  - xi*((Solution(j-(Order-1)/2 + 1,Order) - Solution(j-(Order-1)/2,Order))*(CR*Xi));
			}
//#pragma omp barrier
			Solution = Aux;
		}

	}

	cout << sqrt((Solution - init_state).sqnorm()/grid.size_x) <<"    ";

}

void WENO::CWENO2D(Matrix<double> init_state, Grid_2D grid, int Order, int time_step, double delta_t, double Vel_X, double Vel_Y){
	Matrix<double> Solution(grid.size_x,grid.size_y);
	Matrix<double> Aux(grid.size_x,grid.size_y);
	double xshift = Vel_X * delta_t / grid.delta_x;
	double yshift = Vel_Y * delta_t / grid.delta_y;
	// backward length

	int rotate_x = floor(xshift + 0.5);
	int rotate_y = floor(yshift + 0.5);
	double xi_x = xshift - rotate_x;
	double xi_y = yshift - rotate_y;

//	Vector<double> Sub_1(Order);
//	Vector<double> Sub_2(Order);
	Vector<double> Xi_x(Order);
	Vector<double> Xi_y(Order);

	MakeXi(xi_x, Xi_x, Order);
	MakeXi(xi_y, Xi_y, Order);

	Matrix<double> CL(Order,Order);
	Matrix<double> CR(Order,Order);
	Assign(CL,CR,Order);

	Solution = init_state;

	int i,j,k;
	for (i = 0 ; i < time_step; i++){
		// rows
		for (j = 0; j < grid.size_x; j++){
#pragma omp parallel for private(k) schedule(static)
			for (k = 0; k < grid.size_y; k++){
				if (xi_x > 0){
					Aux(j+rotate_x, k) = Solution(j,k) - xi_x*((Solution.slicec(k)(j-(Order-1)/2, Order) -
							Solution.slicec(k)(j-(Order+1)/2, Order))*(CL*Xi_x));
				}
				else{
					Aux(j+rotate_x, k) = Solution(j,k)  - xi_x*((Solution.slicec(k)(j-(Order-1)/2 + 1,Order) -
							Solution.slicec(k)(j-(Order-1)/2,Order))*(CR*Xi_x));
				}
			}
		}
//#pragma omp barrier
		Solution = Aux;
		// columns
		for (j = 0; j < grid.size_x; j++){
#pragma omp parallel for private(k) schedule(static)
			for (k = 0; k < grid.size_y; k++){
				if (xi_y > 0){
					Aux(j, k+rotate_y) = Solution(j,k) - xi_y*((Solution.slicer(j)(k-(Order-1)/2,Order) -
							Solution.slicer(j)(k-(Order+1)/2, Order))*(CL*Xi_y));
				}
				else{
					Aux(j, k+rotate_y) = Solution(j,k) - xi_y*((Solution.slicer(j)(k-(Order-1)/2 + 1,Order) -
							Solution.slicer(j)(k-(Order-1)/2,   Order))*(CR*Xi_y));
				}
			}
		}
//#pragma omp barrier
		Solution = Aux;

	}
	cout << sqrt((Solution - init_state).sqnorm()/grid.size_x/grid.size_y) << "    ";
}

void WENO::CWENO3D(Tensor<double> init_state, Grid_3D grid, int Order, int time_step, double delta_t, double Vel_X, double Vel_Y, double Vel_Z){
	Tensor<double> Solution(grid.size_x,grid.size_y, grid.size_z);
	Tensor<double> Aux(grid.size_x,grid.size_y, grid.size_z);
	double xshift = Vel_X * delta_t / grid.delta_x;
	double yshift = Vel_Y * delta_t / grid.delta_y;
	double zshift = Vel_Z * delta_t / grid.delta_z;

	int rotate_x = floor(xshift + 0.5);
	int rotate_y = floor(yshift + 0.5);
	int rotate_z = floor(zshift + 0.5);
	double xi_x = xshift - rotate_x;
	double xi_y = yshift - rotate_y;
	double xi_z = zshift - rotate_z;

	Vector<double> Xi_x(Order);
	Vector<double> Xi_y(Order);
	Vector<double> Xi_z(Order);

	MakeXi(xi_x, Xi_x, Order);
	MakeXi(xi_y, Xi_y, Order);
	MakeXi(xi_z, Xi_z, Order);

	Matrix<double> CL(Order,Order);
	Matrix<double> CR(Order,Order);
	Assign(CL,CR,Order);

	Solution = init_state;

	int i,j,k,l;
	for (i = 0 ; i < time_step; i++){
		// rows
		for (j = 0; j < grid.size_x; j++){
			for (k = 0; k < grid.size_y; k++){
#pragma omp parallel for private(l) schedule(static)
				for (l = 0; l < grid.size_z; l++){
					if (xi_x > 0){
						Aux(j + rotate_x, k, l) = Solution(j,k,l) - xi_x*((Solution.slicech(k,l)(j-(Order-1)/2, Order) -
								Solution.slicech(k,l)(j-(Order+1)/2, Order))*(CL*Xi_x));
					}
					else{
						Aux(j+rotate_x, k,l) = Solution(j,k,l)  - xi_x*((Solution.slicech(k,l)(j-(Order-1)/2 + 1,Order) -
								Solution.slicech(k,l)(j-(Order-1)/2,Order))*(CR*Xi_x));
					}
				}
			}
		}
//#pragma omp barrier
		Solution = Aux;
		// columns
		for (j = 0; j < grid.size_x; j++){
			for (k = 0; k < grid.size_y; k++){
#pragma omp parallel for private(l) schedule(static)
				for (l = 0; l < grid.size_z; l++){
					if (xi_y > 0){
						Aux(j, k+rotate_y,l) = Solution(j,k,l) - xi_y*((Solution.slicerh(j,l)(k-(Order-1)/2,Order) -
								Solution.slicerh(j,l)(k-(Order+1)/2, Order))*(CL*Xi_y));
					}
					else{
						Aux(j, k+rotate_y,l) = Solution(j,k,l) - xi_y*((Solution.slicerh(j,l)(k-(Order-1)/2 + 1,Order) -
								Solution.slicerh(j,l)(k-(Order-1)/2,   Order))*(CR*Xi_y));
					}
				}
			}
		}
//#pragma omp barrier
		Solution = Aux;
		// heights
		for (j = 0; j < grid.size_x; j++){
			for (k = 0; k < grid.size_y; k++){
#pragma omp parallel for private(l) schedule(static)
				for (l = 0; l < grid.size_z; l++){
					if (xi_z > 0){
						Aux(j, k,l + rotate_z) = Solution(j,k,l) - xi_z*((Solution.slicerc(j,k)(l-(Order-1)/2,Order) -
								Solution.slicerc(j,k)(l-(Order+1)/2, Order))*(CL*Xi_z));
					}
					else{
						Aux(j, k,l + rotate_z) = Solution(j,k,l) - xi_z*((Solution.slicerc(j,k)(l-(Order-1)/2 + 1,Order) -
								Solution.slicerc(j,k)(l-(Order-1)/2,   Order))*(CR*Xi_z));
					}
				}
			}
		}
//#pragma omp barrier
		Solution = Aux;

	}
	cout << sqrt((Solution - init_state).sqnorm()/grid.size_x/grid.size_y/grid.size_z) << "    ";
}

void WENO::Assign(Matrix<double>& CL, Matrix<double>& CR, int Order){
	if (Order == 3){
		CL(0,0) = -1./6; CL(0,1) = 0.   ; CL(0,2) = 1./6;
		CL(1,0) = 5./6 ; CL(1,1) = 1./2 ; CL(1,2) = -1./3;
		CL(2,0) = 1./3 ; CL(2,1) = -1./2; CL(2,2) = 1./6;

		CR(0,0) = 1./3 ; CR(0,1) = -1./2; CR(0,2) = 1./6;
		CR(1,0) = 5./6 ; CR(1,1) = 1./2 ; CR(1,2) = -1./3;
		CR(2,0) = -1./6; CR(2,1) = 0.   ; CR(2,2) = 1./6;
	}

	if (Order == 5){
		CL(0,0) = 1./30; CL(0,1) = 0. ;CL(0,2) = -1./24; CL(0,3) = 0.;CL(0,4) = 1./120;
		CL(1,0) = -13./60; CL(1,1) = -1./24; CL(1,2) = 1./4; CL(1,3) = 1./24; CL(1,4) = -1./30;
		CL(2,0) = 47./60; CL(2,1) = 5./8; CL(2,2) = -1./3; CL(2,3) = -1./8; CL(2,4) = 1./20;
		CL(3,0) = 9./20; CL(3,1) = -5./8; CL(3,2) = 1./12; CL(3,3) = 1./8; CL(3,4) = -1./30;
		CL(4,0) = -1./20; CL(4,1) = 1./24; CL(4,2) = 1./24; CL(4,3) = -1./24; CL(4,4) = 1./120;

		for (int i = 0; i < Order; i++){
			for (int j = 0; j < Order; j++){
				CR(i,j) = CL(Order-1-i, j);
			}
		}

	}


	if (Order == 7){
		CL(0,0) = -1. / 140 ;CL(0,1) = 0. ;CL(0,2) = 7. / 720;CL(0,3) = 0;CL(0,4) = -1. / 360;CL(0,5) = 0.;CL(0,6) = 1. / 5040;
		CL(1,0) =5. / 84 ;CL(1,1) =1./ 180 ;CL(1,2) = -19. / 240 ;CL(1,3) = -1. / 144;CL(1,4) =  1. / 48 ;CL(1,5) = 1. / 720 ;CL(1,6) =  -1. / 840;
		CL(2,0) =-101./ 420 ;CL(2,1) = -5. / 72;CL(2,2) = 7. / 24;CL(2,3) = 11. / 144 ;CL(2,4) = -13. / 240;CL(2,5) = -1. / 144;CL(2,6) =1./ 336 ;
		CL(3,0) = 319. / 420;CL(3,1) = 49. / 72 ;CL(3,2) = -23. / 72;CL(3,3) = -7. / 36;CL(3,4) = 23. / 360;CL(3,5) = 1./ 72;CL(3,6) = -1. / 252;
		CL(4,0) = 107. / 210 ;CL(4,1) = -49. / 72;CL(4,2) = 1. / 48;CL(4,3) = 7. / 36;CL(4,4) = -1./ 30;CL(4,5) =-1. / 72 ;CL(4,6) = 1. / 336;
		CL(5,0) =-19. / 210 ;CL(5,1) =  5. / 72;CL(5,2) = 7. / 80;CL(5,3) = -11./ 144;CL(5,4) =1. / 240 ;CL(5,5) = 1. / 144;CL(5,6) = -1. / 840;
		CL(6,0) = 1. / 105;CL(6,1) =  -1. / 180;CL(6,2) = -1./ 90;CL(6,3) =  1. / 144;CL(6,4) = 1. / 720;CL(6,5) = -1. / 720;CL(6,6) =1. / 5040 ;

		for (int i = 0; i < Order; i++){
			for (int j = 0; j < Order; j++){
				CR(i,j) = CL(Order-1-i, j);
			}
		}
	}

	if (Order == 9){
		CL(0,0)=1. / 630      ; CL(0,1)= 0.          ; CL(0,2)= -41. / 18144    ; CL(0,3)= 0.           ; CL(0,4)= 13. / 17280  ; CL(0,5)= 0.          ; CL(0,6)= -1. / 12096; CL(0,7)= 0.          ; CL(0,8)= 1./ 362880;
		CL(1,0)=-41. / 2520   ; CL(1,1)= -1. / 1120  ; CL(1,2)=  2081. / 90720  ; CL(1,3)= 7. / 5760    ; CL(1,4)= -1./ 135     ; CL(1,5)= -1. / 2880  ; CL(1,6)= 23. / 30240; CL(1,7)= 1. / 40320  ; CL(1,8)= -1. / 45360;
		CL(2,0)=199./ 2520    ; CL(2,1)= 17. / 1440  ; CL(2,2)= -281. / 2592    ; CL(2,3)= -89. / 5760  ; CL(2,4)= 139. / 4320  ; CL(2,5)= 11./ 2880   ; CL(2,6)= -17. / 6048; CL(2,7)= -1. / 5760  ; CL(2,8)= 1. / 12960;
		CL(3,0)=-641. / 2520  ; CL(3,1)= -127./ 1440 ; CL(3,2)= 4097. / 12960   ; CL(3,3)= 587. / 5760  ; CL(3,4)= -29. / 432   ; CL(3,5)=  -41. / 2880; CL(3,6)= 167./ 30240; CL(3,7)= 1. / 1920   ; CL(3,8)= -1. / 6480;
		CL(4,0)= 1879. / 2520 ; CL(4,1)= 205. / 288  ; CL(4,2)= -797./ 2592     ; CL(4,3)= -91. / 384   ; CL(4,4)= 587. / 8640  ; CL(4,5)= 5. / 192    ; CL(4,6)= -19. / 3024; CL(4,7)= -1./ 1152   ; CL(4,8)= 1. / 5184 ;
		CL(5,0)= 275. / 504 ; CL(5,1)= -205. / 288   ; CL(5,2)= -59. / 2592  ; CL(5,3)= 91./ 384     ; CL(5,4)= -29. / 1080 ; CL(5,5)= -5. / 192  ; CL(5,6)= 25. / 6048  ; CL(5,7)= 1. / 1152; CL(5,8)=    -1./ 6480 ;
		CL(6,0)= -61. / 504  ; CL(6,1)= 127. / 1440  ; CL(6,2)= 1637. / 12960; CL(6,3)= -587./ 5760  ; CL(6,4)= -17. / 4320 ; CL(6,5)= 41. / 2880 ; CL(6,6)= -43. / 30240; CL(6,7)= -1. / 1920; CL(6,8)= 1./ 12960 ;
		CL(7,0)= 11. / 504   ; CL(7,1)= -17. / 1440  ; CL(7,2)= -491. / 18144; CL(7,3)= 89. / 5760   ; CL(7,4)= 11./ 2160   ; CL(7,5)= -11. / 2880; CL(7,6)= 1. / 6048   ; CL(7,7)= 1. / 5760; CL(7,8)= -1. / 45360;
		CL(8,0)= -1./ 504    ; CL(8,1)=  1. / 1120   ; CL(8,2)= 59. / 22680  ; CL(8,3)= -7. / 5760   ; CL(8,4)= -11. / 17280; CL(8,5)= 1./ 2880   ; CL(8,6)= 1. / 60480  ; CL(8,7)= -1. / 40320; CL(8,8)= 1. / 362880;
		for (int i = 0; i < Order; i++){
			for (int j = 0; j < Order; j++){
				CR(i,j) = CL(Order-1-i, j);
			}
		}

	}
}

void WENO::MakeXi(double _xi, Vector<double>& _Xi, int Order){
	// make Xi
	_Xi(0) = 1.0;
	for (int i = 1; i < Order; i++){
		_Xi(i) = _Xi(i-1)*fabs(_xi);
	}
}
