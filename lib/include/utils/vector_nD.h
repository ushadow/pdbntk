/*
 * vector_nD.h
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MVECTOR_ND_H
#define MVECTOR_ND_H

#include "utils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdarg>
#include "mdarray.h"

namespace mocapy {
	inline int min(int a, int b) { return ((a<=b) ? a : b ); }
	
	// n-dimensional vector class
	class Vector_nD {
	private:
		 double *data;
	public:
		 int size;
	
		 Vector_nD(): data(NULL){
			  size = -1;
		 }
	
		 // Constructor: Only allocate memory
		 Vector_nD(int size): data(NULL) {
			  init(size);
		 }
	
		 // THIS METHOD HAS SOME ISSUES WHEN USED WITH NON DIGIT FLOAT CONSTANTS. e.g. Vector_nD v(2,10,10) - always call using Vector_nD v(2, 10.0, 10.0)
		 // Constructor: Initialize based on variable argument list
		 // First element is special to avoid confusion with previous constructor
		 Vector_nD(int size, double firstElement, ...): data(NULL) {
			  init(size);
	
			  va_list ap;
			  this->data[0] = firstElement;
			  va_start(ap,firstElement);
			  for (int i=1; i<size; i++) {
				   this->data[i] = va_arg(ap, double);
			  }
			  va_end(ap);
		 }
	
	
		 // Constructor: From array
		 Vector_nD(int size, double *data): data(NULL) {
			  init(size);
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = data[i];
			  }
		 }
	
		 // Constructor: From vector
		 Vector_nD(std::vector<double> data): data(NULL) {
			  init(data.size());
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = data[i];
			  }
		 }
	
		 // Constructor: From pointer vector
		 Vector_nD(std::vector<double *> data): data(NULL) {
			  init(data.size());
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = *data[i];
			  }
		 }
	
		 // Constructor: From 1D MDArray
		 Vector_nD(mocapy::MDArray<double> &array): data(NULL) {
			  init(array.get_shape()[0]);
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = array[i];
			  }
		 }
	
		 // Copy constructor
		 Vector_nD(const Vector_nD &v): data(NULL) {
			  init(v.size);
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = v.data[i];
			  }
		 }
	
		  // Persistence
		  friend class boost::serialization::access;
		  template<class Archive>
		  void load(Archive & ar, const unsigned int version);
		  template<class Archive>
		  void save(Archive & ar, const unsigned int version) const;
		  BOOST_SERIALIZATION_SPLIT_MEMBER()
	
		 // Initialize
		 void init(int size) {
			  this->size = size;
		  if (data!=NULL)
			   delete[] data;
			  this->data = new double[this->size];
		 }
	
	
		 // Clear contents
		 void clear() {
			  this->size = -1;
			  if (this->data)
				   delete[] this->data;
		  this->data = NULL;
		 }
	
		 // Destructor
		 ~Vector_nD() {
			  if (data != NULL) {
				   delete[] data;
			  }
		  this->data = NULL;
		 }
	
		 // Check whether vector is initialized
		 bool initialized() {
			  return (size!=-1);
		 }
	
		 //Overload indexing operator
		 double operator[](const int index) const {
			  return data[index];
		 }
	
		 // Overload indexing operator (non-const)
		 double& operator[](const int index) {
			  return data[index];
		 }
	
		 // Overload + operator (with vector)
		 Vector_nD operator+(const Vector_nD& v2) const {
			  Vector_nD res(min(this->size, v2.size));
			  for (int i=0; i<res.size; i++)  {
				   res.data[i] = data[i] + v2.data[i];
			  }
			  return res;
		 }
	
		 // Overload + operator (with value)
		 Vector_nD operator+(const double value) const {
			  Vector_nD res(this->size);
			  for (int i=0; i<res.size; i++)  {
				   res.data[i] = data[i] + value;
			  }
			  return res;
		 }
	
		 // Overload - operator (with vector)
		 Vector_nD operator-(const Vector_nD& v2) const {
			  Vector_nD res(min(this->size, v2.size));
			  for (int i=0; i<res.size; i++)  {
				   res.data[i] = this->data[i] - v2.data[i];
			  }
			  return res;
		 }
	
		 // Overload - operator (with value)
		 Vector_nD operator-(const double value) const {
			  Vector_nD res(this->size);
			  for (int i=0; i<res.size; i++)  {
				   res.data[i] = this->data[i] - value;
			  }
			  return res;
		 }
	
		 // Overload - operator (unary)
		 Vector_nD operator-() const {
			  Vector_nD res(size);
			  for (int i=0; i<res.size; i++)  {
				   res.data[i] = -data[i];
			  }
			  return res;
		 }
	
		 // Overload * operator (dot product)
		 double operator*(const Vector_nD& v2) const {
			  double res = 0.0;
			  assert(this->size == v2.size);
			  for (int i=0; i<size; i++)  {
				   res += data[i] * v2.data[i];
			  }
			  return res;
		 }
	
		 // Overload * operator (with value)
		 Vector_nD operator*(const double value) const {
			  Vector_nD res(this->size);
			  for (int i=0; i<res.size; i++)  {
				   res[i] = data[i] * value;
			  }
			  return res;
		 }
	
		 // Overload / operator
		 Vector_nD operator/(const double value) const {
			  Vector_nD res(this->size);
			  for (int i=0; i<res.size; i++)  {
				   res[i] = data[i] / value;
			  }
			  return res;
		 }
	
		 // Overload assignment operator (vector)
		 const Vector_nD& operator=(const Vector_nD& v2) {
			  if (this->size != v2.size) {
				   if (this->size != -1) {
						if (data) {
							 delete[] data;
						}
				data = NULL;
				   }
				   init(v2.size);
			  }
			  for (int i=0; i<size; i++)  {
				   data[i] = v2.data[i];
			  }
			  return *this;
		 }
	
		 // Overload assignment operator (array)
		 const Vector_nD& operator=(const double *array) {
			  assert(size != -1);
			  for (int i=0; i<size; i++)  {
				   data[i] = array[i];
			  }
			  return *this;
		 }
	
		 // Overload assignment operator (value)
		 const Vector_nD& operator=(const double value) {
			  assert(size != -1);
			  for (int i=0; i<size; i++)  {
				   data[i] = value;
			  }
			  return *this;
		 }
	
		 // Overload in-place += operator (with vector)
		 const Vector_nD& operator+=(const Vector_nD& v2) {
			  for (int i=0; i<v2.size; i++)  {
				   data[i] += v2.data[i];
			  }
			  return *this;
		 }
	
		 // Overload in-place += operator (with value)
		 const Vector_nD& operator+=(const double value) {
			  for (int i=0; i<size; i++)  {
				   data[i] += value;
			  }
			  return *this;
		 }
	
		 // Overload in-place -= operator (with vector)
		 const Vector_nD& operator-=(const Vector_nD& v2) {
			  for (int i=0; i<v2.size; i++)  {
				   data[i] -= v2.data[i];
			  }
			  return *this;
		 }
	
		 // Overload in-place -= operator (with value)
		 const Vector_nD& operator-=(const double value) {
			  for (int i=0; i<size; i++)  {
				   data[i] -= value;
			  }
			  return *this;
		 }
	
		 // Overload in-place *= operator
		 const Vector_nD& operator*=(const double value) {
			  for (int i=0; i<size; i++)  {
				   data[i] *= value;
			  }
			  return *this;
		 }
	
		 // Overload in-place /= operator
		 const Vector_nD& operator/=(const double value) {
			  for (int i=0; i<size; i++)  {
				   data[i] /= value;
			  }
			  return *this;
		 }
	
		 // Overload > operator
		 bool operator>(const Vector_nD& v2) {
			  return (this->norm() > v2.norm());
		 }
	
		 // Overload >= operator
		 bool operator>=(const Vector_nD& v2) {
			  return (this->norm() >= v2.norm());
		 }
	
		 // Overload < operator
		 bool operator<(const Vector_nD& v2) {
			  return (this->norm() < v2.norm());
		 }
	
		 // Overload <= operator
		 bool operator<=(const Vector_nD& v2) {
			  return (this->norm() <= v2.norm());
		 }
	
		 // Calculate length
		 double norm() const {
			  double res = 0.0;
			  for (int i=0; i<size; i++)  {
				   res += data[i]*data[i];
			  }
			  return sqrt(res);
		 }
	
		 // Normalize vector
		 Vector_nD normalize() const {
			  double norm = this->norm();
			  if (norm != 0) {
				   return *this/norm;
			  } else {
				   return *this;
			  }
		 }
	
		 // Return array used internally
		 double *getArray() {
			  return this->data;
		 }
	
		 // Return array as vector
		 std::vector<double> getVector() {
			  return (std::vector<double>(this->data, this->data+this->size));
		 }
	
		 friend Vector_nD range(double firstElement, double lastElement, double stepSize);
		 friend double max(const Vector_nD &v);
		 friend double min(const Vector_nD &v);
		 friend double sum(const Vector_nD &v);
		 friend double mean(const Vector_nD &v);
		 friend Vector_nD cumSum(const Vector_nD &v);
		 friend Vector_nD fmap(double(*f)(double), const Vector_nD &v);
	
		 // Output
		 friend std::ostream& operator<<(std::ostream& Out, const Vector_nD &v) {
			  Out << "(";
			  for (int i=0; i<v.size; i++) {
				   if (i>0) {
						Out << ", ";
				   }
				   Out << v.data[i];
			  }
			  Out << ")";
			  return Out;
		 }
	};
	
	
	// Create vector initialized with range of numbers
	inline Vector_nD range(double firstElement, double lastElement, double stepSize) {
		 int size = (int)fabs(ceil((lastElement-firstElement)/stepSize));
		 Vector_nD res(size);
		 int direction;
		 if (firstElement < lastElement) {
		  direction = 1;
		 } else {
		  direction = -1;
		 }
		 for (int i=0; i<size; i++) {
		  res.data[i] = firstElement + stepSize * i * direction;
		 }
		 return res;
	}
	
	// Find maximum
	inline double max(const Vector_nD &v) {
		 double maxValue = v.data[0];
		 for (int i=1; i<v.size; i++) {
		  if (v.data[i] > maxValue) {
			   maxValue = v.data[i];
		  }
		 }
		 return maxValue;
	}
	
	// Find minimum
	inline double min(const Vector_nD &v) {
		 double minValue = v.data[0];
		 for (int i=1; i<v.size; i++) {
		  if (v.data[i] < minValue) {
			   minValue = v.data[i];
		  }
		 }
		 return minValue;
	}
	
	// return sum of elements
	inline double sum(const Vector_nD &v) {
		 double sum = 0.0;
		 for (int i=0; i<v.size; i++) {
		  sum += v.data[i];
		 }
		 return sum;
	}
	
	// return sum of elements
	inline double mean(const Vector_nD &v) {
		 double sum = 0.0;
		 for (int i=0; i<v.size; i++) {
		  sum += v.data[i];
		 }
		 return sum/v.size;
	}
	
	// Return new vector containin cumulative sum of the elements so far
	inline Vector_nD cumSum(const Vector_nD &v) {
		 Vector_nD res(v.size);
		 double sum = 0.0;
		 for (int i=0; i<v.size; i++) {
		  sum += v.data[i];
		  res.data[i] = sum;
		 }
		 return res;
	}
	
	
	// Apply function to all elements
	inline Vector_nD fmap(double(*f)(double), const Vector_nD &v) {
		 Vector_nD res(v.size);
		 for (int i=0; i<v.size; i++) {
		  res.data[i] = (*f)(v.data[i]);
		 }
		 return res;
	}
	
	template<class Archive>
	 void Vector_nD::save(Archive & ar, const unsigned int version) const {
		ar & size;
		int i;
		for(i = 0; i < size; ++i)
			 ar & data[i];
	  }
	
	template<class Archive>
	 void Vector_nD::load(Archive & ar, const unsigned int version) {
		ar & size;
		data = new double[size];
		int i;
		for(i = 0; i < size; ++i)
			 ar & data[i];
	
	 }

} // end namespace

#endif
