/*
 * Matrix.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: lurker
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "Vector.hpp"

using namespace std;

template <class T> class Matrix{
private:

public:
	int dim_x;
	int dim_y;
	valarray<valarray<T>> storage;
	Matrix(){
	}
	Matrix(int arg_1, int arg_2){
		storage.resize(arg_1);
		int i;
//#pragma omp parallel for private(i) schedule(static)
		for (i = 0; i < arg_1; i++){
			storage[i].resize(arg_2);
		}
		dim_x = arg_1;
		dim_y = arg_2;
	}
	~Matrix(){

	}
	void Setup(int arg_1, int arg_2){
		storage.resize(arg_1);
		int i;
//#pragma omp parallel for private(i) schedule(static)
		for (i = 0; i < arg_1; i++){
			storage[i].resize(arg_2);
		}
		dim_x = arg_1;
		dim_y = arg_2;
	}
	Matrix& operator+=(const Matrix& rhs){
		(*this).storage += rhs.storage;
		return *this;
	}

	Matrix& operator-=(const Matrix& rhs){
		(*this).storage -= rhs.storage;
		return *this;
	}
	Matrix& operator+=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			(*this).storage[i] += scalar;
		}
		return *this;
	}
	Matrix& operator-=(const double scalar){
		for (int i = 0; i< (*this).dim_x; i++){
			(*this).storage[i] += scalar;
		}
		return *this;
	}

	Matrix operator+(const Matrix& rhs){
		Matrix lhs(rhs.dim_x, rhs.dim_y);
		for(int i = 0 ;i < rhs.dim_x; i++){
			lhs.storage[i] = (*this).storage[i] + rhs.storage[i];
		}

		return lhs;
	}
	Matrix operator-(const Matrix& rhs){
		Matrix lhs(rhs.dim_x, rhs.dim_y);
		for(int i = 0 ;i < rhs.dim_x; i++){
			lhs.storage[i] = (*this).storage[i] - rhs.storage[i];
		}

		return lhs;
	}
	Matrix& operator*=(double scalar){
		for (int i = 0 ; i < (*this).dim_x; i++)
		 (*this).storage[i] *= scalar;
		 return *this;
	}
	Matrix& operator/=(double scalar){
		for (int i = 0 ; i < (*this).dim_x; i++)
		 (*this).storage[i] /= scalar;
		 return *this;
	}
	Matrix operator*(double scalar){
		Matrix lhs((*this).dim_x, (*this).dim_y);
		lhs.storage = (*this).storage;
		for (int i = 0; i < lhs.dim_x; i++){
			lhs.storage[i] *= scalar;
		}
		return lhs;
	}
	Matrix operator/(double scalar){
		Matrix lhs((*this).dim_x, (*this).dim_y);
		lhs.storage = (*this).storage;
		for (int i = 0; i < lhs.dim_x; i++){
			lhs.storage[i] /= scalar;
		}
		return lhs;
	}

	Matrix& operator= (const Matrix rhs){
		(*this).dim_x = rhs.dim_x;
		(*this).dim_y = rhs.dim_y;
		(*this).storage = rhs.storage;
		return *this;
	}
	T& operator()(const int location_1, const int location_2);
	Matrix operator()(const int lox, const int loy, const int num_x, const int num_y){
		Matrix sub(num_x, num_y);
		for (int i = 0; i < num_x ; i++){
			for (int j = 0; j < num_y; j++){
				sub(i,j) = (*this)(lox+i,loy+j);
			}
		}
		return sub;
	}

	Vector<T> slicer(const int location){
		Vector<T> sub((*this).dim_y);
		for (int i= 0 ; i< (*this).dim_y; i++){
			sub(i) = (*this)(location,i);
		}
		return sub;
	}

	Vector<T> slicec(const int location){
		Vector<T> sub((*this).dim_x);
		for (int i= 0 ; i< (*this).dim_x; i++){
			sub(i) = (*this)(i,location);
		}
		return sub;
	}

	Vector<T> slicer(const int location, const int start, const int length){
		Vector<T> sub(length);
		int i;
//#pragma omp parallel for
		for (i = 0 ; i < length ; i++){
			sub(i) = (*this)(location, start + i);
		}
		return sub;
	}

	Vector<T> slicec(const int location, const int start, const int length){
		Vector<T> sub(length);
		int i;
//#pragma omp parallel for
		for (i=0; i < length; i++){
			sub(i) = (*this)(start + i, location);
		}
		return sub;
	}

	Vector<T> operator*(const Vector<T> lhs){
		Vector<T> sub((*this).dim_x);
		for (int i = 0; i < (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				sub(i) = sub(i) + (*this).storage[i][j]*lhs.storage[j];
			}
		}
		return sub;
	}

	double sqnorm(){
			double ret = 0;
			for (int i=0 ; i < (*this).dim_x; i++){
				for (int j=0 ; j < (*this).dim_y; j++){
				ret += (*this)(i,j)*(*this)(i,j);
				}
			}
			return ret;
		}
		
	double L1norm(){
			double ret = 0;
			for (int i=0 ; i < (*this).dim_x; i++){
				for (int j=0 ; j < (*this).dim_y; j++){
				ret += fabs((*this)(i,j));
				}
			}
			return ret;
		}
		
	double L2norm(){
			return sqrt((*this).sqnorm());
		}

	void Write(string filename){
		ofstream sol;
		sol.open(filename.c_str());
		for (int i=0; i < (*this).dim_x; i++){
			for (int j = 0; j < (*this).dim_y; j++){
				sol << (*this)(i,j) << "\t";
			}
			sol <<"\n";
		}
		sol.close();
	}
};

template <class T> T& Matrix<T>::operator ()(const int location_1, const int location_2){
	return storage[((location_1) % dim_x + dim_x) % dim_x][((location_2 ) % dim_y + dim_y) % dim_y];
}


#endif /* MATRIX_HPP_ */
