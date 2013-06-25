/*
 * mdarray.h
 *
 *  Copyright (C) 2008, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef MDARRAY_H_
#define MDARRAY_H_

// For Boost serialization
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <assert.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include "../framework/mocapyexceptions.h"
#include "randomgen.h"
#include "utils.h"


namespace mocapy {

// Forwards
template<typename T>
class MDArray;




template<typename T> std::ostream& operator<<(std::ostream&, const MDArray<T>&);

 template<typename T> T sumvect(std::vector<T> & v);

#define __CLPK_integer int
#define __CLPK_doublereal double

extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
extern "C" void dgetrf_(const int*, const int*, double*, const int*, int*, int*);
extern "C" void dgetri_(const int*, double*, const int*, const int*, double*, const int*, int*);

template<typename T>
class MDArray {
public:
	// Constructors
	MDArray() {};
	MDArray(const MDArray& a) {
		copy(a);
	}
	MDArray(const std::vector<uint> & shape) {
		set_shape(shape);
	}

	MDArray(const std::vector<int> & shape) {
	  std::vector<uint> sh;
		for (uint i=0; i<shape.size(); i++) sh.push_back(shape[i]);
		set_shape(sh);
	}

	// Destructor
	virtual ~MDArray() {
		clear();
	}

	// Getting/setting shape
	const std::vector<uint> get_shape() const {
		return shape;
	}
	void set_shape(const std::vector<uint> & new_shape, bool reset = true);
	void set_shape(uint a, uint b = INF, uint c = INF, uint d = INF);

	// Getters
	std::vector<T>& get_values() {
		return values;
	}
	std::vector<T> get_values_flat();

	const std::vector<uint> get_own_index() {
		return ownIndex;
	}

	MDArray<T> & get_view(const std::vector<uint> & indices);
	MDArray<T>& get_view(uint a);
	inline T & get(uint a, uint b);
	inline T & get(uint a, uint b, uint c);
	T get_max();
	void get_max(T & m);
	std::vector<T> get_slice(uint from, uint to);

	// Setters
	void set(uint a, T new_value);
	void set(uint a, uint b, T new_value);
	void set(uint a, uint b, uint c, T new_value);
	void set(uint a, MDArray<T> & new_value);
	void set_values(const std::vector<T> & a);
	void set_values(T * a);

	// set_wildcard assigns specified value to the elements in the index vector.
	// The index vector can have -1 valued indices corresponding to wildcards.
	// e.g. set_wildcard with index=[3,-1,0], value=5, assigns integer value 5
	// to [3,0,0], [3,1,0], [3,2,0], ...
	// There can be an arbitrary number of wildcards.
	void set_wildcard(const std::vector<int> & index, T value, int depth=0);

	// In a 2-dim matrix M, repeat copies 'a' to M[i] for all i
	void repeat(uint s, const std::vector<T> & a);

	void div_inplace(T a);
	void div_inplace(MDArray<T> & a);
	void add_inplace(T a);
	void add_inplace(MDArray<T> & a);
	void add_inplace(std::vector<std::vector<double> > & a);
	void add_inplace(std::vector<double> & a);
	void sub_inplace(MDArray<T> & a);
	void log_all(); //natural log
	void exp_all(); //natural exp
	void sqr_all(); //square the elements (x*x)
	void sqrt_all(); //square root of the elements (x^.5)
	uint bisect(double b);
	MDArray<T> sumlast();
	void cumsum(bool safe=true);
	void normalize();
	void transpose();
	MDArray<T> swapaxes(uint axis1, uint axis2);
	MDArray<T> moveaxis(uint from, uint to);
	void rnd_cov(uint dim, RandomGen * rg);
	void multiply_inplace(T m);
	void multiply_inplace(MDArray<T> & m);
	inline T dot(MDArray<T> & m);
//	inline vector<T> dot(vector<T> & v);
	inline void mult3(MDArray<T> & m);
	double det();
	double det_old();
//	double determinant2();

	void inv();
	void inv_old();

	//Operators
	friend std::ostream& operator<< <T> (std::ostream& output, const MDArray<T>& a);
	MDArray<T>& operator=(const MDArray<T>& a);
	inline T & operator[](std::vector<uint> & indices);
	inline T & operator[](uint a);
	MDArray<T> operator/(T a);
	MDArray<T> operator+(MDArray<T> & a);
	MDArray<T> operator-(MDArray<T> & a);
	std::vector<T> operator*(std::vector<T> & v);
	MDArray<T> operator*(MDArray<T> & v);

	// Misc.
	void clear();
	uint size();
	bool empty() const;
	void randomize( RandomGen * rg, double maxrnd=1.0, bool thin=false);
	std::vector<T*> flat_iterator();
	void flat_array(T* v);
	void clip(double minimum, double maximum);
	void copy(const MDArray<T> & array);
	void copy_cast(const MDArray<double>& array);

	void eye(uint dim, T r = 1);
	void keep_eye();
	void makeRotationMatrix(double angle, MDArray<double> & v);
	void makeRefMatrix(MDArray<double> & u, MDArray<double> & v);
	void makeRotationMatrix(MDArray<double> & u, MDArray<double> & v);

	// Calculate eigen values and vectors of a symmetric matrix
	std::pair<std::vector<double>, std::vector<MDArray<double> > > eigen();

    // String representation
	std::string tostring(
        const uint precision=2,
        const uint width=0,
        const std::string sep_field=" ",
        const std::string sep_record="\n",
        const std::string sep_record_outer="\n"
        ) const;

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    friend class MDArray<eMISMASK>;

protected:
	void initialize(bool reset = true);
	void flat_iterator(std::vector<T*> & v);
	void set(std::vector<uint> & indices, T new_value, uint index);
	inline T & get(std::vector<uint> & indices, uint index);
	MDArray<T>& get_view(const std::vector<uint> & indices, uint index);
	void sumlast(MDArray<T> & sumArray, uint sumDim, uint currentDim);
	void cumsum(uint currentDim, uint dimension, bool safe);
	void transpose(double **a, int n);
	void permuteaxes(const std::vector<uint> & permutation, std::vector<uint> & new_index, MDArray<T> & new_array);
	void CoFactor(double **a, int n, double **b);
	double determinant(double **a, int n);

	MDArray LU(bool & neg, __CLPK_integer *pivots=NULL);

	std::vector<MDArray<T>*> A;
	std::vector<T> values;
	std::vector<uint> ownIndex;
	std::vector<uint> shape;
};

// String representation

template<typename T>
  std::string MDArray<T>::tostring(
        const uint precision,
        const uint width,
        const std::string sep_field,
        const std::string sep_record,
        const std::string sep_record_outer
    ) const
{
    // Create local recursive print function
    struct Local {
         static void recursive_print(
                MDArray<T> a,
                std::ostringstream & os,
                std::string sep_field,
                std::string sep_record,
                std::string sep_record_outer
                )
         {
	   std::vector<uint> shape = a.get_shape();
	   int len = shape.size();
	   if (len == 1)
	     {
                // We've reached the last dimension
                int dim = shape[0];
                for (int i=0; i<dim; i++)
                {
                    os << a[i];
                    if (i<dim-1) { os << sep_field; }
                }
                os << sep_record;
            }
            else
            {
                int dim = shape[0];
                for(int i=0; i<dim; i++) {
                    if (i==dim-1) { sep_record=""; }
                    recursive_print(a.get_view(i), os, sep_field, sep_record, sep_record_outer);
                }
                // ...n x m x p array will be printed in m x p blocks
                os << sep_record_outer;
            }
         }
    };

    std::ostringstream os;
    if (width>0) { os << std::setw(width); }
    os << std::setprecision(precision); // max precision digits after .
    os << std::fixed; // force use of precision digits after .

    if (empty())
        os << "Empty" << sep_record;
    else
        Local::recursive_print(*this, os, sep_field, sep_record, sep_record_outer);

    return os.str();
}

// *** 	Getting/setting shape ***

template<typename T>
  void MDArray<T>::set_shape(const std::vector<uint> & new_shape, bool reset) {
	clear();
	shape = new_shape;
	initialize(reset);
}

// Syntactic sugar
template<typename T>
void MDArray<T>::set_shape(uint a, uint b, uint c, uint d) {
  std::vector<uint> sh;
  sh.push_back(a);
  if (b < INF) {
    sh.push_back(b);
    if (c < INF)
      sh.push_back(c);
    if (d < INF)
      sh.push_back(d);
    
  }
  
  set_shape(sh);
 }

//*** Getters ***

template<typename T>
  std::vector<T> MDArray<T>::get_slice(uint from, uint to) {
  std::vector<T> ret;
	for(uint i = from; i != to; i++) {
		ret.push_back(values[i]);
	}
	return ret;
}

template<typename T>
  std::vector<T> MDArray<T>::get_values_flat() {
  std::vector<T*> b = flat_iterator();
  std::vector<T> v;
	for (uint i=0; i<b.size(); i++) {
	    v.push_back( *(b[i]) );
	}
	return v;
}

template<typename T>
inline T & MDArray<T>::operator[](uint a) {
	assert(get_shape().size() == 1); // Otherwise use the other operators or get_view
	assert(a < values.size());
	return values[a];
}

template<typename T>
inline T & MDArray<T>::get(uint a, uint b) {
	assert(a < A.size());
	assert(b < A[a]->values.size());
	T & ref = A[a]->values[b];
	return ref;
}

template<typename T>
inline T & MDArray<T>::get(uint a, uint b, uint c) {
	return this->get_view(a).get(b,c);
}

template<typename T>
  inline T & MDArray<T>::operator[](std::vector<uint> & indices) {
	if (indices.size() == 2) {
		return A[indices.front()]->values[indices.back()];
	}
	else if (indices.empty()) {
		assert(!values.empty());
		return values[0];
	}
	else if (indices.size() == 1) {
		assert(indices.front() < values.size() );
		return values[indices.front()];
	}
	else {
		return get(indices, 0);
	}
}

template<typename T>
  inline T & MDArray<T>::get(std::vector<uint> & indices, uint index) {
	if ((index + 1) == indices.size()) {
		assert(index < indices.size());
		assert(indices[index] < values.size());
		return values[indices[index]];
	} else {
		assert(index<indices.size());
		assert(indices[index] < A.size());
		return A[indices[index]]->get(indices, index + 1);
	}
}

template<typename T>
  MDArray<T>& MDArray<T>::get_view(const std::vector<uint> & indices) {
	assert(indices.size() <= shape.size());
	if(indices.empty())
		return *this;
	else
		return get_view(indices, 0);
}


template<typename T>
  MDArray<T>& MDArray<T>::get_view(const std::vector<uint> & indices, uint index) {
	if (index == indices.size()) {
		return *this;
	} else {
		assert(index < indices.size());
		assert(indices[index] < A.size());
		return A[indices[index]]->get_view(indices, index + 1);
	}
}

template<typename T>
MDArray<T>& MDArray<T>::get_view(uint a) {
	assert(a<A.size());
	return *A[a];
}

// *** Setters ***

template<typename T>
void MDArray<T>::set(uint a, T new_value) {
	assert(a < values.size());
	values[a] = new_value;
}


template<typename T>
void MDArray<T>::set(uint a, uint b, T new_value) {
	assert(a < A.size());
	assert(b < A[a]->size());
	A[a]->values[b] = new_value;
}

template<typename T>
void MDArray<T>::set(uint a, uint b, uint c, T new_value) {
	this->get_view(a).set(b,c,new_value);
}

template<typename T>
void MDArray<T>::set(uint a, MDArray<T> & new_value) {
	assert(a < A.size());
	A[a]->copy(new_value);
}

template<typename T>
  void MDArray<T>::set(std::vector<uint> & indices, T new_value, uint index) {
	if (index + 1 == indices.size()) {
		assert(!values.empty());
		assert(index < indices.size());
		assert(values.size()> indices[index]);
		values[indices[index]] = new_value;
	} else {
		assert(index < indices.size());
		assert(indices[index] < A.size());
		A[indices[index]]->set(indices, new_value, index + 1);
	}
}

template<typename T>
  void MDArray<T>::set_values(const std::vector<T> & a) {
  std::vector<T*> b = flat_iterator();
	assert(b.size() == a.size());

	for (uint i=0; i<a.size(); i++) {
		*(b[i]) = a[i];
	}
}


template<typename T>
  void MDArray<T>::set_values(T * a) {
  std::vector<T*> b = flat_iterator();
  for (uint i=0; i<b.size(); i++) {
    *(b[i]) = a[i];
  }
}



// *** Arithmetics


// Compute the cummulative sum over the last dimension of an array
// ex: |4 6 1|     |4 10 11|
//     |3 1 0|  =  |3  4  4|
//     |7 7 2|     |7 14 16|
// If safe is true, then last column is set to INF to avoid rounding error problems in bisect
template<typename T>
void MDArray<T>::cumsum(bool safe) {
	uint d = shape.size();
	cumsum(0, d, safe);
}

template<typename T>
void MDArray<T>::cumsum(uint currentDim, uint dimension, bool safe) {
	if (dimension == currentDim + 1) {
		for (uint i = 1; i < values.size(); i++) {
			values[i] += values[i - 1];
		}
		if (safe) {
			assert(!values.empty());
			uint i=values.size()-1;
			values[i] = INF;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->cumsum(currentDim + 1, dimension, safe);
		}
	}
}


// I'm not sure that this works...
// I need to test this! /MP
template<typename T>
MDArray<T> MDArray<T>::sumlast() {
	assert(!shape.empty());
	uint sumDim = shape.size() - 1;
	MDArray<T> sumArray;
	std::vector<uint> sh = shape;
	assert(!sh.empty());
	sh.pop_back();
	sumArray.set_shape(sh);
	sumlast(sumArray, sumDim, 0);
	return sumArray;
}

template<typename T>
void MDArray<T>::sumlast(MDArray<T> & sumArray, uint sumDim, uint currentDim) {
	if (sumDim == currentDim) {
		assert(!values.empty());
		T sum(0);
		for (uint i = 0; i < values.size(); i++) {
			sum += values[i];
		}
		std::vector<uint> indx = ownIndex;
		sumArray.set(indx, sum);
	}

	for (uint i = 0; i < A.size(); i++) {
		A[i]->sumlast(sumArray, sumDim, currentDim + 1);
	}
}

template<typename T>
void MDArray<T>::log_all() {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			double v = values[i];
			if (v>0)
				values[i] = log(v);
			else
				values[i] = -INF;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->log_all();
		}
	}
}

template<typename T>
void MDArray<T>::exp_all() {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] = exp(values[i]);
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->exp_all();
		}
	}
}

template<typename T>
void MDArray<T>::sqr_all() {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] = values[i] * values[i];
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->sqr_all();
		}
	}
}

template<typename T>
void MDArray<T>::sqrt_all() {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] = srqt(values[i]);
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->sqrt_all();
		}
	}
}

template<typename T>
void MDArray<T>::normalize() {
	if (shape.size() == 1) {
		T sum(0);
		for (uint i = 0; i < values.size(); i++) {
			sum += values[i];
		}

		for (uint i = 0; i < values.size(); i++) {
			if (sum != 0)
				values[i] /= sum;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->normalize();
		}
	}
}

template<typename T>
void MDArray<T>::div_inplace(T m) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] /= m;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->div_inplace(m);
		}
	}
}



template<typename T>
void MDArray<T>::div_inplace(MDArray<T> & denum) {
	if (shape.size() == 1) {
	  std::vector<uint> & indx = ownIndex;
		for (uint i = 0; i < values.size(); i++) {
			T d = denum.get(indx);
			values[i] /= d;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->div_inplace(denum);
		}
	}
}

template<typename T>
void MDArray<T>::add_inplace(MDArray<T> & a) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] += a.values[i];
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			MDArray<T> * new_a = a.A[i];
			A[i]->add_inplace(*new_a);
		}
	}
}

template<typename T>
  void MDArray<T>::add_inplace(std::vector<std::vector<double> > & a) {
	for (uint i = 0; i < A.size(); i++) {
		MDArray<T> * ref = A[i];
		for (uint j = 0; j < ref->values.size(); j++) {
			ref->values[j] += a[i][j];
		}
	}
}

template<typename T>
  void MDArray<T>::add_inplace(std::vector<double> & a) {
	for (uint i = 0; i < values.size(); i++) {
		values[i] += a[i];
	}
}


template<typename T>
void MDArray<T>::sub_inplace(MDArray<T> & a) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] -= a.values[i];
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			MDArray<T> * new_a = a.A[i];
			A[i]->sub_inplace(*new_a);
		}
	}
}

template<typename T>
MDArray<T> MDArray<T>::operator+(MDArray<T> & a) {
	MDArray r(*this);
	r.add_inplace(a);
	return r;
}

template<typename T>
MDArray<T> MDArray<T>::operator-(MDArray<T> & a) {
	MDArray r(*this);
	r.sub_inplace(a);
	return r;
}

template<typename T>
MDArray<T> MDArray<T>::operator/(T a) {
	MDArray r(*this);
	r.div_inplace(a);
	return r;
}

template<typename T>
void MDArray<T>::add_inplace(T s) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] += s;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->add_inplace(s);
		}
	}
}

template<typename T>
  void MDArray<T>::rnd_cov(uint dim, RandomGen * rg) {
	/*
	 Return a random covariance matrix with shape (dim x dim).
	 From U{http://www.sci.wsu.edu/math/faculty/genz/papers/mvn/node5.html}:

	 These random covariance matrices were generated with a method described in an
	 article by Marsaglia and Olkin, for m = 3, 4, ..., 10. This method
	 consists of the following steps. First, a random lower triangular matrix is
	 generated with entries chosen uniformly from [-1,1]. Then, the rows of this
	 matrix are normalized so that the sum of the squares of the elements in each
	 row add to one. Finally this lower triangular matrix is multiplied by it's
	 transpose, to produce a random covariance matrix.

	 Note:
	 Lower triangular matrix:
	 A matrix which is only defined at (i,j) when i greater than or equal to j.

	 Reference:
	 Marsaglia, G. and Olkin, I. (1984), Generating Correlation Matrices, SIAM
	 Journal of Scientific and Statistical Computing 5, pp. 470-475.
	 */

  std::vector<uint> sh_c;
	sh_c.push_back(dim);
	sh_c.push_back(dim);
	MDArray<double> c(sh_c);

	for (uint i = 0; i < dim; i++) {
		for (uint j = 0; j < dim; j++) {
			if (i >= j) {
				double d = rg->get_rand();
				double r = 2.0 * (d - 0.5);
				c.set(i, j, r);
			}
		}
	}
	MDArray<double> d(c);
	d.sqr_all();
	for (uint i = 0; i < dim; i++) {
		MDArray<double> & values = d.get_view(i);
		std::vector<double> & v_values = values.get_values();
		double s = sumvect(v_values);
		MDArray<double> view = c.get_view(i);
		view.multiply_inplace(1.0 / sqrt(s));
		c.set(i, view);
	}
	MDArray<double> cov(c);
	c.transpose();
	cov.multiply_inplace(c);
	copy(cov);
}

template<typename T>
void MDArray<T>::multiply_inplace(T m) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] *= m;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->multiply_inplace(m);
		}
	}
}

template<typename T>
void MDArray<T>::multiply_inplace(MDArray<T> & m) {
	assert(shape.size() == 2);
	assert(m.get_shape().size() == 2);

	uint dim = shape.front();

	MDArray<T> cp(shape);

	for (uint i = 0; i < dim; i++) {
		for (uint j = 0; j < dim; j++) {
			T v = 0;
			for (uint k = 0; k < dim; k++) {
				v += get(i, k) * m.get(k, j);
			}
			cp.set(i, j, v);
		}
	}
	copy(cp);
}

template<typename T>
void MDArray<T>::transpose() {
	assert(shape.size() == 2);
	uint dim = shape.front();

	for (uint i = 0; i < dim; i++) {
		for (uint j = i + 1; j < dim; j++) {
			T temp = get(i, j);
			set(i, j, get(j, i));
			set(j, i, temp);
		}
	}
}



template<typename T>
  void printv(std::vector<T> v) {
  std::cout << "{";
  for(uint i=0; i<v.size(); i++)
    std::cout << v[i] << " ";
  std::cout << "}" << std::endl;
}

template<typename T>
MDArray<T> MDArray<T>::swapaxes(uint axis1, uint axis2) {
	const uint dims = shape.size();
	assert(axis1 < dims && axis2 < dims);

	if(axis1==axis2) {
		MDArray<T> new_array(*this);
		return new_array;
	}

	// Calculate the new shape vector
	std::vector<uint> new_shape(shape);
	new_shape[axis1] = shape[axis2];
	new_shape[axis2] = shape[axis1];

	// Construct the new MDArray
	MDArray<T> new_array(new_shape);

	// Calculate the permutation of the axis
	std::vector<uint> permutation;
	for(uint i=0; i<dims; i++) {permutation.push_back(i);}
	permutation[axis1] = axis2;
	permutation[axis2] = axis1;

	// Swap the axis
	std::vector<uint> new_index(dims);
	permuteaxes(permutation, new_index, new_array);

	return new_array;
}

template<typename T>
MDArray<T> MDArray<T>::moveaxis(uint from, uint to) {
	const uint dims = shape.size();
	assert(from < dims && to < dims);

	if(from==to) {
		MDArray<T> new_array(*this);
		return new_array;
	}

	// Calculate the inverse permutation of the axis
	std::vector<uint> inv_permutation;
	for(uint i=0; i<dims; i++) {inv_permutation.push_back(i);}
	inv_permutation.erase(inv_permutation.begin()+from);
	inv_permutation.insert(inv_permutation.begin()+to, from);

	// Calculate the true permutation
	std::vector<uint> permutation(dims);
	for(uint i=0; i<dims; i++) {permutation[inv_permutation[i]] = i;}

	// Calculate the new shape vector
	std::vector<uint> new_shape;
	for(std::vector<uint>::iterator axis=inv_permutation.begin(); axis!=inv_permutation.end(); axis++) {
		new_shape.push_back(shape[*axis]);
	}

	// Construct the new MDArray
	MDArray<T> new_array(new_shape);

	// Swap the axis
	std::vector<uint> new_index(dims);
	permuteaxes(permutation, new_index, new_array);

	return new_array;
}


// Recursive method to permute the axes of this and return the result in new_array
// It assumes that permutation that 'permutation' is a correct permutation!
template<typename T>
  void MDArray<T>::permuteaxes(const std::vector<uint> & permutation, std::vector<uint> & new_index, MDArray<T> & new_array) {
	uint depth = ownIndex.size();
	uint new_depth = permutation[depth];

	// If we are at full depth, copy the value to the new index in the new array
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			new_index[new_depth] = i;
			new_array[new_index] = values[i];
		}
	}
	// If we are not at a full depth, set the new index correct and call permute on all sub arrays
	else {
		for (uint i = 0; i < A.size(); i++) {
			new_index[new_depth] = i;
			A[i]->permuteaxes(permutation, new_index, new_array);
		}
	}
}


template<typename T>
double MDArray<T>::det_old() {
	assert(!shape.empty());
	int n = shape.front();

	double **m = NULL;
	m = (double**) malloc((n) * sizeof(double*));
	for (int i = 0; i < n; i++) {
		m[i] = (double*) malloc((n) * sizeof(double));
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double v = get(i, j);
			m[i][j] = v;
		}
	}

	double d = determinant(m, n);

	for (int i = 0; i < n; i++)
		free(m[i]);
	free(m);

	return d;
}



// Recursive definition of determinate using expansion by minors.
template<typename T>
double MDArray<T>::determinant(double **a, int n) {
	int i, j, j1, j2;
	double det = 0;
	double **m = NULL;

	if (n < 1) { /* Error */
		assert(false);

	} else if (n == 1) { /* Shouldn't get used */
		det = a[0][0];
	} else if (n == 2) {
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	} else {
		det = 0;
		for (j1 = 0; j1 < n; j1++) {
			m = (double**) malloc((n - 1) * sizeof(double *));
			for (i = 0; i < n - 1; i++)
				m[i] = (double*) malloc((n - 1) * sizeof(double));
			for (i = 1; i < n; i++) {
				j2 = 0;
				for (j = 0; j < n; j++) {
					if (j == j1)
						continue;
					m[i - 1][j2] = a[i][j];
					j2++;
				}
			}
			det += pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * determinant(m, n - 1);
			for (i = 0; i < n - 1; i++)
				free(m[i]);
			free(m);
		}
	}
	return (det);
}


template<typename T>
void MDArray<T>::inv_old() {
	assert(!shape.empty());
	int n = shape.front();

	double **m = NULL;
	double **m2 = NULL;
	m = (double**) malloc((n) * sizeof(double*));
	m2 = (double**) malloc((n) * sizeof(double*));
	for (int i = 0; i < n; i++) {
		m[i] = (double*) malloc((n) * sizeof(double));
		m2[i] = (double*) malloc((n) * sizeof(double));
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double v = get(i, j);
			m[i][j] = v;
		}
	}

	CoFactor(m, n, m2);
	transpose(m2, n);
	double d = determinant(m, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double v = m2[i][j] / d;
			set(i, j, v);
		}
	}
	for (int i = 0; i < n; i++) {
		free(m[i]);
		free(m2[i]);

	}
	free(m);
	free(m2);
}


// Find the cofactor matrix of a square matrix
template<typename T>
void MDArray<T>::CoFactor(double **a, int n, double **b) {
	int i, j, ii, jj, i1, j1;
	double det;
	double **c;

	c = (double**) malloc((n - 1) * sizeof(double *));
	for (i = 0; i < n - 1; i++)
		c[i] = (double*) malloc((n - 1) * sizeof(double));

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {

			// Form the adjoint a_ij
			i1 = 0;
			for (ii = 0; ii < n; ii++) {
				if (ii == i)
					continue;
				j1 = 0;
				for (jj = 0; jj < n; jj++) {
					if (jj == j)
						continue;
					c[i1][j1] = a[ii][jj];
					j1++;
				}
				i1++;
			}

			// Calculate the determinate
			det = determinant(c, n - 1);

			// Fill in the elements of the cofactor
			b[i][j] = pow(-1.0, i + j + 2.0) * det;
		}
	}
	for (i = 0; i < n - 1; i++)
		free(c[i]);
	free(c);
}

//   Transpose of a square matrix, do it in place
template<typename T>
void MDArray<T>::transpose(double **a, int n) {
	int i, j;
	double tmp;

	for (i = 1; i < n; i++) {
		for (j = 0; j < i; j++) {
			tmp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = tmp;
		}
	}
}

template<typename T>
  inline std::vector<T> MDArray<T>::operator*(std::vector<T> & m) {
	assert(shape.size() == 2);
	uint dim = shape.front();

	std::vector<T> r(dim);

	for (uint i = 0; i < dim; i++) {
		T v = 0;
		std::vector<T> & vref = A[i]->values;
		for (uint k = 0; k < dim; k++) {
			v += vref[k] * m[k];
		}
		r[i] = v;
	}
	return r;
}

template<typename T>
inline MDArray<T> MDArray<T>::operator*(MDArray<T> & m) {
	assert(shape.size() == 2);
	uint dim = shape.front();

	MDArray<T> r;
	r.set_shape(dim);

	for (uint i = 0; i < dim; i++) {
		T v = 0;
		std::vector<T> & vref = A[i]->values;
		for (uint k = 0; k < dim; k++) {
			v += vref[k] * m[k];
		}
		r[i] = v;
	}
	return r;
}




template<typename T>
inline T MDArray<T>::dot(MDArray<T> & m) {
	T v(0);
	for (uint i = 0; i < m.values.size(); i++) {
		v += values[i] * m.values[i];
	}
	return v;
}

// Multiplication of 3x3 matrices
template<typename T>
inline void MDArray<T>::mult3(MDArray<T> & m) {
     assert(shape.size() == 2);
     uint dim=shape.front();
     assert(dim == 3 || dim == 1);

     // M*v
     if(dim==1) {
       std::vector<T> vref = A[0]->values; //Need a copy

          for(uint i=0; i<m.shape.front();i++) {
               set(0,i,vref[0]*m.get(i,0)+vref[1]*m.get(i,1)+vref[2]*m.get(i,2)); }
     }
     // M*M
     else {
          MDArray<T> cp(shape);

          for (uint i = 0; i < dim; i++) {
	    std::vector<T> & vref = A[i]->values;
               for (uint j = 0; j < dim; j++) {
                    T sum = 0;
                    for (uint k = 0; k < dim; k++) {
                         sum += vref[k] * m.get(k,j);
                    }
                    cp.set(i,j,sum);
               }
          }

          copy(cp);
     }
}

// *** Operators ***

template<typename T>
  std::ostream& operator<<(std::ostream& output, const MDArray<T>& a) {
	if (a.A.empty()) {
		for (uint i = 0; i < a.values.size(); i++) {
			output << a.values[i] << " ";
		}
	} else {
		for (uint i = 0; i < a.A.size(); i++) {
		  output << "[" << *(a.A[i]) << "]" << std::endl;
		}
	}
	return output;
}

template<typename T>
MDArray<T>& MDArray<T>::operator= (const MDArray<T>& a) {
	if (this != &a) {
		copy(a);
	}
	return *this;
}

// *** Misc. ***

template<typename T>
bool MDArray<T>::empty() const {
	return shape.empty();
}

template<typename T>
void MDArray<T>::copy(const MDArray<T> & array) {
	clear();
	shape = array.shape;
	ownIndex = array.ownIndex;
	values = array.values;
	for (uint i = 0; i < array.A.size(); i++) {
		MDArray<T>* newA = new MDArray<T> ;
		newA->copy(*array.A[i]);
		A.push_back(newA);
	}
}

template<typename T>
void MDArray<T>::copy_cast(const MDArray<double>& array) {
	clear();
	shape = array.shape;
	ownIndex = array.ownIndex;

	for (uint i=0; i<array.values.size(); i++) {
		values.push_back((T)array.values[i]);
	}

	for (uint i = 0; i < array.A.size(); i++) {
		MDArray<T>* newA = new MDArray<T> ;
		newA->copy_cast(*array.A[i]);
		A.push_back(newA);
	}
}

template<typename T>
void MDArray<T>::clear() {
	if (shape.size() > 2) {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->clear();
		}
	}

	for (uint i = 0; i < A.size(); i++) {
		delete A[i];
	}
	A.clear();
	values.clear();
}

template<typename T>
void MDArray<T>::clip(double minimum, double maximum) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			values[i] = max(minimum, values[i]);
			values[i] = min(maximum, values[i]);
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->clip(minimum, maximum);
		}
	}
}

 template<typename T>
   void MDArray<T>::initialize(bool reset) {
   assert(!shape.empty());
   
   if (shape.size() > 1) {
     std::vector<uint> new_shape(shape.size() - 1);
     for (uint i = 1; i < shape.size(); i++) {
       new_shape[i - 1] = shape[i];
     }
     
     for (uint i = 0; i < shape.front(); i++) {
       MDArray<T>* a = new MDArray<T> ;
       a->ownIndex = ownIndex;
       a->ownIndex.push_back(i);
       a->set_shape(new_shape); // also calls initialize
       A.push_back(a);
     }
   }
   
   else if (shape.size() == 1) {
     values.resize(shape.front());
   }
   else {
     assert(shape.empty());
     values.resize(1);
   }
 }
 
template<typename T>
uint MDArray<T>::size() {
	if (shape.empty())
		return 0;

	uint sz(1);
	for (uint i = 0; i < shape.size(); i++) {
		sz *= shape[i];
	}
	return sz;
}

template<typename T>
  void MDArray<T>::randomize( RandomGen * rg, double maxrnd, bool thin) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {

		  double d = rg->get_rand()*maxrnd;
			if (thin) {
				double x = rg->get_rand();

				if (x<0.7)
					d=0.01*d;
			}


			values[i] = d;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
		  A[i]->randomize(rg, maxrnd, thin);
		}
	}
}


// Locate the proper insertion point for item in list to maintain sorted order.
// If item is already present in list, the insertion point will be AFTER any existing entries.
template<typename T>
uint MDArray<T>::bisect(double b) {
	// TODO: This should run in log time for long vectors
	assert(!values.empty());

	for (uint i = 0; i < values.size(); i++) {
		if (values[i] > b)
			return i;
	}
	assert(false);
	return 0;
}

template<typename T>
  std::vector<T*> MDArray<T>::flat_iterator() {
  std::vector<T*> v;
	flat_iterator(v);
	return v;
}

template<typename T>
  void MDArray<T>::flat_iterator(std::vector<T*> & v) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			v.push_back(&values[i]);
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->flat_iterator(v);
		}
	}
}


template<typename T>
void MDArray<T>::flat_array(T* v) {
	assert(shape.size() == 2);
	assert(shape[0] = shape[1]);

	uint N = shape[0];
	for (uint i=0; i<N; i++) {
		for (uint j=0; j<N; j++) {
			v[i+N*j] = get(i,j);
		}
	}
}



// set_wildcard assigns specified value to the elements in the index vector.
// The index vector can have -1 valued indices corresponding to wildcards.
// e.g. set_wildcard with index=[3,-1,0], value=5, assigns integer value 5
// to [3,0,0], [3,1,0], [3,2,0], ...
// There can be an arbitrary number of wildcards
template<typename T>
  void MDArray<T>::set_wildcard(const std::vector<int> & index, T value, int depth) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			if (index[depth] == (int)-1 || index[depth] == (int)i)
				values[i] = value;
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			if (index[depth] == (int)-1 || index[depth] == (int)i)
				A[i]->set_wildcard(index, value, depth+1);
		}
	}
}


// In a 2-dim matrix M, repeat copies 'a' to M[i] for all i
template<typename T>
  void MDArray<T>::repeat(uint s, const std::vector<T> & a) {
	set_shape(s, a.size());

	for (uint i=0; i<s; i++) {
		get_view(i).set_values(a);
	}
}


template<typename T>
T MDArray<T>::get_max() {
	T m = -INF;
	if (m>0) { m=0; } // manages the unsigned int case
	get_max(m);
	return m;
}

template<typename T>
void MDArray<T>::get_max(T & m) {
	if (shape.size() == 1) {
		for (uint i = 0; i < values.size(); i++) {
			m = std::max(m, values[i]);
		}
	} else {
		for (uint i = 0; i < A.size(); i++) {
			A[i]->get_max(m);
		}
	}
}

template<typename T>
void MDArray<T>::eye(uint dim, T r) {
	clear();
	std::vector<uint> sh;
	sh.push_back(dim);
	sh.push_back(dim);
	set_shape(sh);

	for (uint i = 0; i < dim; i++) {
		assert(i < A.size());
		assert(i < A[i]->values.size());
		A[i]->values[i] = r;
	}
}

template<typename T>
void MDArray<T>::keep_eye() {
	assert(!empty());
	uint dim = shape.front();
	for (uint i = 0; i < dim; i++) {
		for (uint j = 0; j < dim; j++) {
			assert(i < A.size());
			assert(j < A[i]->values.size());
			if (i != j)
				A[i]->values[j] = 0;
		}
	}
}


// Calculate eigen values and vectors of a general matrix
template<typename T>
  std::pair<std::vector<double>, std::vector<MDArray<double> > > MDArray<T>::eigen() {
	// N': left eigenvectors of A are not computed
	char JOBVL = 'N';

	// 'V': right eigenvectors of A are computed.
	char JOBVR = 'V';

	assert(shape.size() == 2);
	__CLPK_integer N = shape[0];

	__CLPK_doublereal eigenVectorsL[N*N];
	__CLPK_doublereal eigenVectorsR[N*N];
	__CLPK_doublereal eigenValues[N];
	__CLPK_doublereal eigenValues_img[N];

	__CLPK_doublereal A[N*N];
	flat_array(A);
	__CLPK_integer LDA = N;

	__CLPK_doublereal *WR = eigenValues;
	__CLPK_doublereal *WI = eigenValues_img;
	__CLPK_doublereal *VR = eigenVectorsR;
	__CLPK_doublereal *VL = eigenVectorsL;

	__CLPK_integer INFO;
	__CLPK_integer LWORK = -1;

	__CLPK_doublereal WORK_tmp;
	dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &N, VR, &N, &WORK_tmp, &LWORK, &INFO);
	LWORK = (int)WORK_tmp;
	__CLPK_doublereal *WORK = new __CLPK_doublereal[LWORK];
	dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &N, VR, &N, WORK, &LWORK, &INFO);
	std::vector<double> eigenVal;
	std::vector< MDArray<double> > eigenVec;

	for (int i=0; i<N; i++) {
		eigenVal.push_back(eigenValues[i]);
		MDArray<double> vec;
		vec.set_shape(N);
		for (int j=0; j<N; j++) {
			vec[j] = eigenVectorsR[j + i*N];
		}
		eigenVec.push_back(vec);
	}

	delete[] WORK;
	return make_pair(eigenVal, eigenVec);
 }


template<typename T>
void MDArray<T>::makeRotationMatrix(double angle, MDArray<double> & v) {
	set_shape(3, 3);
	assert(v.get_shape().size() == 1 && v.get_shape()[0] == 3);

	double ux = v[0];
	double uy = v[1];
	double uz = v[2];

	double a00 = ux*ux + cos(angle)*(1.0-ux*ux);
	double a10 = ux*uy*(1.0-cos(angle))+uz*sin(angle);
	double a20 = uz*ux*(1.0-cos(angle)) - uy*sin(angle);

	double a01 = ux*uy*(1.0-cos(angle)) - uz*sin(angle);
	double a11 = uy*uy + cos(angle)*(1.0-uy*uy);
	double a21 = uy*uz*(1.0-cos(angle)) + ux*sin(angle);

	double a02 = uz*ux*(1.0-cos(angle)) + uy*sin(angle);
	double a12 = uy*uz*(1.0-cos(angle)) - ux*sin(angle);
	double a22 = uz*uz + cos(angle)*(1.0 - uz*uz);

	set(0,0, a00);
	set(0,1, a01);
	set(0,2, a02);

	set(1,0, a10);
	set(1,1, a11);
	set(1,2, a12);

	set(2,0, a20);
	set(2,1, a21);
	set(2,2, a22);
}


// LU decomposition of a general matrix
template<typename T>
  MDArray<T> MDArray<T>::LU(bool & neg, __CLPK_integer *pivots) {

	uint nRow = shape[0];
	uint nCol = shape[1];
	if (nRow != nCol) {
		throw MocapyExceptions("Error (LU): Cannot compute LU decomposition for non-square matrix");
	}


	__CLPK_integer M = nRow;
	__CLPK_integer N = nCol;
	__CLPK_integer LDA = M;


	__CLPK_doublereal A[N*N];
	flat_array(A);

	__CLPK_integer info;

	int minMN = M;
	if (N<M)
		minMN = N;

	bool newPivots = false;
	if (!pivots) {
		pivots = new __CLPK_integer[minMN];
		newPivots = true;
	}

	// Call LAPACK function
	dgetrf_(&M, &N, A, &LDA, pivots, &info);

	// Check error status
	if (info<0) {
		throw MocapyExceptions("Determinant: Illegal value");
	}
	else if (info>0) {
	  //	  std::cout << *this << "\n";
	  //	  throw MocapyExceptions("Zero determinant");

	  MDArray<T> newMatrix = *this;
	  return newMatrix;
	} else {
	  neg=false;
	  /* Take the product of the diagonal elements */
	  for (int c1 = 0; c1 < minMN; c1++) {
	    if (pivots[c1] != (c1+1)) neg = !neg;
	  }
	  
	  if (newPivots)
	    delete[] pivots;
	  
	  MDArray<T> newMatrix = *this;
	  newMatrix.set_values(A);
	  return newMatrix;	 
	}

}


template<typename T>
double MDArray<T>::det() {
  uint nRow = shape[0];

  // Expansion by minors (faster for small matrices)
  if (nRow <=4)
    return det_old();

  // Compute determinant using LU decomposition
  bool neg;
  MDArray<T> newMatrix = LU(neg);
  
  double d=1;
  
  for (uint i=0;i<nRow;i++) {
    d*=newMatrix.get(i, i);
  }
  
  if (neg)
    d*=-1;
  
  return d;
 
}


template<typename T>
void MDArray<T>::inv() {
  int n = shape[1];

  if (n <= 4) // Faster for small matrices
    return inv_old();
  

  // Use LU decomposition
  int m = n;
  int info, lwork;
  int *ipiv = NULL;
  double *work = NULL;
  int err = 0;
  
  ipiv = (int*) malloc(m * sizeof *ipiv);
  if (ipiv == NULL) {
    fputs("out of memory\n", stderr);
    return ;
  }
  
  work = (double*) malloc(sizeof *work);
  if (work == NULL) {
    fputs("out of memory\n", stderr);
    err = 1;
    goto bailout;
  }   


  __CLPK_doublereal a[n*n];
  flat_array(a);
  
  dgetrf_(&m, &n, a, &m, ipiv, &info);   

  if (info != 0) {
    fprintf(stderr, "dgetrf: info = %d\n", info);
    err = 1;
    goto bailout;
  }
  
  // query for workspace size 
  lwork = -1;
  dgetri_(&n, a, &n, ipiv, work, &lwork, &info);
  
  if (info != 0 || work[0] <= 0.0) {
    fprintf(stderr, "dgetri (workspace): info = %d\n", info);
    err = 1;
    goto bailout;
  }
  
  lwork = (int) work[0];
  
  work = (double*) realloc(work, lwork * sizeof *work);
  if (work == NULL) {
    fputs("out of memory\n", stderr);
    err = 1;
    goto bailout;
  } 
  
  dgetri_(&n, a, &n, ipiv, work, &lwork, &info);
  
  if (info != 0) {
    fprintf(stderr, "dgetri: info = %d\n", info);
    err = 1;
  }
  
 bailout:
  
  free(ipiv);
  free(work);
  set_values(a);
  return ;
}

// Return a (left multiplying) matrix that mirrors p onto q.
template<typename T>
void MDArray<T>::makeRefMatrix(MDArray<double> & u, MDArray<double> & v) {

     assert(v.get_shape().size() == 1 && v.get_shape()[0] == 3);
     assert(u.get_shape().size() == 1 && u.get_shape()[0] == 3);

     set_shape(3,3);
     eye(3);

     MDArray<double> uv = u-v;

     uv.div_inplace(sqrt(uv[0]*uv[0]+uv[1]*uv[1]+uv[2]*uv[2]));

     MDArray<double> B;
     B.set_shape(3,3);
     for (uint i = 0; i < uv.size(); i++) {
          for (uint j = 0; j < uv.size(); j++) {
               double val = uv[i] * uv[j];
               B.set(i, j, val);
          }
     }

     B.multiply_inplace(2.0);

     //A-B;
     set(0,0, 1-B.get(0,0));
     set(0,1,  -B.get(0,1));
     set(0,2,  -B.get(0,2));

     set(1,0,  -B.get(1,0));
     set(1,1, 1-B.get(1,1));
     set(1,2,  -B.get(1,2));

     set(2,0,  -B.get(2,0));
     set(2,1,  -B.get(2,1));
     set(2,2, 1-B.get(2,2));

}

//  Return a left multiplying matrix that rotates u onto v.
template<typename T>
void MDArray<T>::makeRotationMatrix(MDArray<double> & u, MDArray<double> & v) {

     assert(v.get_shape().size() == 1 && v.get_shape()[0] == 3);
     assert(u.get_shape().size() == 1 && u.get_shape()[0] == 3);

     MDArray<double> mu = u;
     mu.multiply_inplace(-1.0);

     MDArray<double> refMat1;
     refMat1.set_shape(3,3);
     refMat1.makeRefMatrix(v,mu);

     MDArray<double> refMat2;
     refMat2.set_shape(3,3);
     refMat2.makeRefMatrix(u,mu);
     refMat1.mult3(refMat2);

     copy(refMat1);

}

template<typename T> template< class Archive>
void MDArray<T>::serialize(Archive & ar, const unsigned int version) {
	ar & A;
	ar & values;
	ar & ownIndex;
	ar & shape;
}
}

#endif /* MDARRAY_H_ */
