/*
 * Vector.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: lurker
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

using namespace std;

template <class T> class Vector{
private:

public:
	int dim;
	valarray<T> storage;
	Vector(int arg){
		storage.resize(arg);
		dim = arg;
	}
	~Vector(){

	}

	Vector& operator*=(double scalar){
		 (*this).storage *= scalar;
		 return *this;
	}
	Vector& operator/=(double scalar){
		 (*this).storage /= scalar;
		 return *this;
	}
	Vector operator*(double scalar){
		Vector lhs((*this).dim);
		lhs.storage = (*this).storage;
		lhs.storage *= scalar;
		return lhs;
	}

	Vector operator/(double scalar){
		Vector lhs((*this).dim);
		lhs.storage = (*this).storage;
		lhs.storage /= scalar;
		return lhs;
	}

	double operator*(const Vector& rhs){
		double lhs = .0;
		for (int i = 0; i < min((*this).dim,rhs.dim); i++){
			lhs += (*this).storage[i]*rhs.storage[i];
		}
		return lhs;
	}

	Vector& operator+=(const Vector& rhs){
		(*this).storage += rhs.storage;
		return *this;
	}

	Vector& operator-=(const Vector& rhs){
		(*this).storage -= rhs.storage;
		return *this;
	}
	Vector& operator+=(const double scalar){
		(*this).storage += scalar;
		return *this;
	}
	Vector& operator-=(const double scalar){
		(*this).storage -= scalar;
		return *this;
	}

	Vector operator+(const Vector& rhs){
		Vector lhs(rhs.dim);
		lhs.storage = (*this).storage + rhs.storage;
		return lhs;
	}

	Vector operator-(const Vector& rhs){
		Vector lhs(rhs.dim);
		lhs.storage = (*this).storage - rhs.storage;
		return lhs;
	}

	Vector& operator= (Vector rhs){
		(*this).dim = rhs.dim;
		(*this).storage = rhs.storage;
		return *this;
	}
	T& operator[](const int location);
	T& operator()(const int location);
	Vector operator()(const int lo, const int num){
		Vector sub(num);
		for (int i = 0; i < num ; i++){
			sub.storage[i] = (*this)(lo+i);
		}
		return sub;
	}
	double sqnorm(){
		double ret = 0;
		for (int i=0 ; i < (*this).dim; i++){
			ret += (*this)(i)*(*this)(i);
		}
		return ret;
	}
};

template <class T> T& Vector<T>::operator ()(const int location){
	return storage[((location) % dim + dim) % dim];
}
// vector access only for 1D vector
template <class T> T& Vector<T>::operator [](const int location){
	if ((location < 0) || (location >= dim)){cout << "subscript exceeds limit" << endl;exit(1);}
	return storage[location];
}


#endif /* VECTOR_HPP_ */
