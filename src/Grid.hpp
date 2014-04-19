/*
 * Grid.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: lurker
 */

#ifndef GRID_HPP_
#define GRID_HPP_


class Grid_1D{
public:
	double start_x, final_x,delta_x;
	int size_x;
	Grid_1D(double a, double b, int s){
		start_x = a;
		final_x = b;
		size_x = s;
		delta_x = (b-a)/s;
	}
	~Grid_1D(){
	}

};

class Grid_2D{
public:
	double start_x, final_x , start_y, final_y, delta_x, delta_y;
	int size_x , size_y;
	Grid_2D(double a, double b , int s, double c, double d,int t){
		start_x = a;start_y = c;
		final_x = b;final_y = d;
		size_x = s; size_y = t;
		delta_x = (b-a)/s; delta_y = (d-c)/t;
	}
	~Grid_2D(){

	}
};

class Grid_3D{
public:
	double start_x, final_x , start_y, final_y, start_z, final_z,delta_x, delta_y,delta_z;
	int size_x , size_y, size_z;
	Grid_3D(double a, double b , int s, double c, double d,int t, double e, double f, int r){
		start_x = a;start_y = c; start_z = e;
		final_x = b;final_y = d; final_z = f;
		size_x = s; size_y = t; size_z = r;
		delta_x = (b-a)/s; delta_y = (d-c)/t; delta_z = (f-e)/r;
	}
	~Grid_3D(){

	}
};

#endif /* GRID_HPP_ */
