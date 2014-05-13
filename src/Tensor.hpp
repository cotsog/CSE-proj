/*
 * Tensor.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: lurker
 */

#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "Matrix.hpp"
#include "Vector.hpp"

template <class T> class Tensor{
private:

public:
	int dim_x;
	int dim_y;
	int dim_z;
	valarray<valarray<valarray<T>>> storage;
	Tensor(int arg_1, int arg_2 , int arg_3){
		storage.resize(arg_1);
		for (int i = 0; i < arg_1; i++){
			storage[i].resize(arg_2);
			for (int j = 0; j < arg_3; j++){
				storage[i][j].resize(arg_3);
			}
		}
		dim_x = arg_1;
		dim_y = arg_2;
		dim_z = arg_3;
	}
	~Tensor(){

	}
	T& operator()(const int location_1, const int location_2, const int location_3);

	Tensor& operator+=(const Tensor& rhs){
		(*this).storage += rhs.storage;
		return *this;
	}

	Tensor& operator-=(const Tensor& rhs){
		(*this).storage -= rhs.storage;
		return *this;
	}
	Tensor& operator+=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				(*this).storage[i][j] += scalar;
			}
		}
		return *this;
	}
	Tensor& operator-=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				(*this).storage[i][j] += scalar;
			}
		}
		return *this;
	}

	Tensor operator+(const Tensor& rhs){
		Tensor lhs(rhs.dim_x, rhs.dim_y, rhs.dim_z);
		for(int i = 0 ;i < rhs.dim_x; i++){
			for (int j = 0; j < rhs.dim_y; j++){
				lhs.storage[i][j] = (*this).storage[i][j] + rhs.storage[i][j];
			}
		}
		return lhs;
	}
	Tensor operator-(const Tensor& rhs){
		Tensor lhs(rhs.dim_x, rhs.dim_y, rhs.dim_z);
		for(int i = 0 ;i < rhs.dim_x; i++){
			for (int j = 0; j < rhs.dim_y; j++){
				lhs.storage[i][j] = (*this).storage[i][j] - rhs.storage[i][j];
			}
		}
		return lhs;
	}
	Tensor& operator*=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				(*this).storage[i][j] *= scalar;
			}
		}
		return *this;
	}
	Tensor& operator/=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				(*this).storage[i][j] /= scalar;
			}
		}
		return *this;
	}

	Tensor operator*(double scalar){
		Tensor lhs((*this).dim_x, (*this).dim_y, (*this).dim_z);
		lhs.storage = (*this).storage;
		for (int i = 0; i < lhs.dim_x; i++){
			for (int j = 0; j < lhs.dim_y; j++){
				lhs.storage[i][j] *= scalar;
			}
		}
		return lhs;
	}
	Tensor operator/(double scalar){
		Tensor lhs((*this).dim_x, (*this).dim_y, (*this).dim_z);
		lhs.storage = (*this).storage;
		for (int i = 0; i < lhs.dim_x; i++){
			for (int j = 0; j < lhs.dim_y; j++){
				lhs.storage[i][j] /= scalar;
			}
		}
		return lhs;
	}
	Tensor& operator= (const Tensor rhs){
		(*this).dim_x = rhs.dim_x;
		(*this).dim_y = rhs.dim_y;
		(*this).dim_z = rhs.dim_z;
		(*this).storage = rhs.storage;
		return *this;
	}
	Vector<T> slicerc(const int location1,const int location2, const int start, const int length){
		Vector<T> sub(length);
		for (int i= 0 ; i< length; i++){
			sub(i) = (*this)(location1,location2,start+i);
		}
		return sub;
	}
	Vector<T> slicerh(const int location1,const int location2 ,const int start, const int length){
		Vector<T> sub(length);
		for (int i= 0 ; i< length; i++){
			sub(i) = (*this)(location1,start+i,location2);
		}
		return sub;
	}
	Vector<T> slicech(const int location1,const int location2, const int start, const int length){
		Vector<T> sub(length);
		for (int i= 0 ; i< length; i++){
			sub(i) = (*this)(i+start,location1,location2);
		}
		return sub;
	}

	Tensor operator()(const int lox, const int loy, const int loz, const int num_x, const int num_y, const int num_z){
		Tensor sub(num_x, num_y, num_z);
		for (int i = 0; i < num_x ; i++){
			for (int j = 0; j < num_y; j++){
				for (int k = 0; k < num_z ;k++){
					sub(i,j,k) = (*this)(lox+i,loy+j,loz + k);
				}
			}
		}
		return sub;
	}

	double sqnorm(){
		double ret = 0;
		for (int i=0 ; i < (*this).dim_x; i++){
			for (int j=0 ; j < (*this).dim_y; j++){
				for (int k = 0; k < (*this).dim_z; k++){
					ret += (*this)(i,j,k)*(*this)(i,j,k);
				}
			}
		}
		return ret;
	}

	double L1norm(){
		double ret = 0;
		for (int i=0 ; i < (*this).dim_x; i++){
			for (int j=0 ; j < (*this).dim_y; j++){
				for (int k = 0; k < (*this).dim_z; k++){
					ret += fabs((*this)(i,j,k));
				}
			}
		}
		return ret;
	}

};

template <class T> T& Tensor<T>::operator ()(const int location_1, const int location_2, const int location_3){
	return storage[((location_1) % dim_x + dim_x) % dim_x][((location_2) % dim_y + dim_y) % dim_y][((location_3) % dim_z + dim_z)%dim_z];
}

#endif /* TENSOR_HPP_ */
