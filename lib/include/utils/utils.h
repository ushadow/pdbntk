/*
 * utils.h
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

#ifndef UTILS_H_
#define UTILS_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <string>
#include <stdlib.h>

#include <vector>

namespace mocapy {
#define _MIN_TRANSITION 1e-50
#define uint unsigned int
#define INF 100000000
#define MOCAPY_ALL -1

enum eMISMASK {
	MOCAPY_OBSERVED, MOCAPY_HIDDEN, MOCAPY_MISSING
};

}


#include "mdarray.h"



namespace mocapy {

typedef MDArray<double> Sequence ;

#define FlatSequence std::vector<double*>

#define CPD MDArray<double>

#define TOINT(x) (int)(x+0.5)

// Forward declaration
class Node;

void printvec(std::vector<uint> & vec);
std::vector<Node*> vec_concNode(std::vector<Node*> & v1, std::vector<Node*> & v2);
template<typename T>
std::vector<T> vec_conc(std::vector<T> & v1, std::vector<T> & v2) {
	std::vector<T> v3;
	v3 = v1;
	for (uint i = 0; i < v2.size(); i++) {
		v3.push_back(v2[i]);
	}
	return v3;
}

template<typename T>
std::vector<T> vec_conc(std::vector<T> & v1, T & v2) {
	std::vector<T> v3;
	v3 = v1;
	v3.push_back(v2);
	return v3;
}

template<typename T>
std::vector<T> operator+(std::vector<T> & v1, std::vector<T> & v2) {
	std::vector<T> v3;
	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		v3.push_back(v1[i] + v2[i]);
	}
	return v3;
}

template<typename T>
std::vector<T> operator-(std::vector<T> & v1, std::vector<T> & v2) {
	assert(v1.size() == v2.size());
	std::vector<T> v3;
	for (uint i = 0; i < v1.size(); i++) {
		v3.push_back(v1[i] - v2[i]);
	}
	return v3;
}

template<typename T>
MDArray<T> toMDArray(std::vector<std::vector<T> > & v) {
	MDArray<T> a;
	uint dim1 = v.size();
	assert(dim1>0);
	uint dim2 = v.front().size();

	a.set_shape(dim1, dim2);
	for (uint i = 0; i < dim1; i++) {
	  if (v[i].size() != dim2)
	    std::cout << "Wrong number of data values in mismask or training data" << std::endl;
		assert(v[i].size() == dim2);
		for (uint j = 0; j < dim2; j++) {
			a.set(i, j, v[i][j]);
		}
	}
	return a;
}

template<typename T, typename Y>
std::vector<T> operator/(std::vector<T> & v1, Y v2) {
	std::vector<T> v3;
	for (uint i = 0; i < v1.size(); i++) {
		v3.push_back(v1[i] / v2);
	}
	return v3;
}

template<typename T>
void add_inplace(std::vector<T> & v1, std::vector<T> & v2) {
	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		v1[i] = v1[i] + v2[i];
	}
}

template<typename T>
inline void vec_sub(std::vector<T> & v1, std::vector<T> & v2, std::vector<T> & v3) {
	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		v3.push_back(v1[i] - v2[i]);
	}
}

template<typename T>
std::vector<T> vec_add(std::vector<T> & v1, T v2) {
	std::vector<T> v3;
	for (uint i = 0; i < v1.size(); i++) {
		v3.push_back(v1[i] + v2);
	}
	return v3;
}

template<typename T>
inline T dot(std::vector<T> & v1, std::vector<T> & v2) {
	T v(0);
	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		v += v1[i] * v2[i];
	}
	return v;
}

template<typename T>
inline void dot2(std::vector<T> & v1, std::vector<T> & v2, MDArray<double> & v3) {
	std::vector<uint> sh(2);
	sh[0] = v1.size();
	sh[1] = v1.size();

	v3.set_shape(sh, false);

	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		for (uint j = 0; j < v1.size(); j++) {
			T val = v1[i] * v2[j];
			v3.set(i, j, val);
		}
	}
}

template<typename T>
inline void dot2(MDArray<T> & v1, MDArray<T> & v2, MDArray<double> & v3) {
	std::vector<uint> sh(2);
	sh[0] = v1.size();
	sh[1] = v1.size();

	v3.set_shape(sh, false);

	assert(v1.size() == v2.size());
	for (uint i = 0; i < v1.size(); i++) {
		for (uint j = 0; j < v1.size(); j++) {
			T val = v1[i] * v2[j];
			v3.set(i, j, val);
		}
	}
}


template<typename T>
inline void dot2(std::vector<T> & v1, std::vector<T> & v2, std::vector<std::vector<T> > & v3) {

	assert(v1.size() == v2.size());

	v3.resize(v1.size());
	for (uint i = 0; i < v1.size(); i++) {
		v3[i].resize(v1.size());
		for (uint j = 0; j < v1.size(); j++) {
			T val = v1[i] * v2[j];
			v3[i][j] = val;
		}
	}
}

template<typename T>
T sumvect(std::vector<T> & v) {
	T s(0);
	for (uint i = 0; i < v.size(); i++) {
		s += v[i];
	}
	return s;
}

template<typename T>
std::vector<T> take(std::vector<T> & seq, std::vector<uint> & index) {
	// Return the values of seq with indices in index
	std::vector<T> r;
	for (uint i = 0; i < index.size(); i++) {
		assert( i < index.size() );
		assert( seq.size()> index[i] );
		r.push_back(seq[index[i]]);
	}
	return r;
}

template<typename T>
inline void take_value(std::vector<T*> & seq, std::vector<uint> & index, std::vector<T> & r) {
	// Return the values of seq with indices in index
	for (uint i = 0; i < index.size(); i++) {
		assert( i < index.size() );
		assert( seq.size()> index[i] );
		r.push_back(*seq[index[i]]);
	}
}

template<typename T>
inline void take_value_toint(std::vector<T*> & seq, std::vector<uint> & index, std::vector<
		uint> & r) {
	// Return the values of seq with indices in index
	for (uint i = 0; i < index.size(); i++) {
		assert( i < index.size() );
		assert( seq.size()> index[i] );
		r.push_back((uint) (*seq[index[i]]));
	}
}

template<typename T>
uint argmax(std::vector<T> & v) {
	T m = -INF;
	uint index = INF;
	for (uint i = 0; i < v.size(); i++) {
		if (v[i] > m) {
			m = v[i];
			index = i;
		}
	}

	assert(index < INF );
	return index;
}

template<typename T>
inline void toint(std::vector<T> & v, std::vector<uint> & v2) {
	v2.resize(v.size());
	for (uint i = 0; i < v.size(); i++) {
		v2[i] = (TOINT(v[i]));
	}
}

std::vector<Node*> vec_conc(std::vector<Node*> & v1, std::vector<Node*> & v2);


// Initialize the random number generators that are used by Mocapy.
void mocapy_seed(uint s);
extern uint moc_seed;
//extern uint moc_seed1;
//extern int moc_seed2;

class NodeID {
public:
	NodeID();
	NodeID(int i) {
		integer = i;
		isInteger = true;
	}
	NodeID(const char* n) {
		name = n;
		isInteger = false;
	}
	bool isInteger;
	int integer;
	std::string name;

	// Persistence
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
};

template<class Archive>
void NodeID::serialize(Archive & ar, const unsigned int version) {
	ar & isInteger;
	ar & integer;
	ar & name;
}

template<typename T>
  std::ostream& operator<<(std::ostream& output, const std::vector<T>& a) {
	for (uint i = 0; i < a.size(); i++) {
		output << a[i];
		if (i + 1 < a.size()) {
			output << ", ";
		}
	}
	return output;
}

template<typename T>
  std::ostream& operator<<(std::ostream& output, const std::pair<T, T>& a) {
	output << "(" << a.first << ", " << a.second << ")";

	return output;
}


template<typename T>
const std::vector<T> vec(T t) {
	std::vector<T> v;
	v.push_back(t);
	return v;
}

template<typename T>
const std::vector<T> vec(T t1, T t2) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	return v;
}

template<typename T>
const std::vector<T> vec(T t1, T t2, T t3) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	return v;
}

template<typename T>
const std::vector<T> vec(T t1, T t2, T t3, T t4) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	v.push_back(t4);
	return v;
}

template<typename T>
const std::vector<T> vec(T t1, T t2, T t3, T t4, T t5) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	v.push_back(t4);
	v.push_back(t5);
	return v;
}

template<typename T>
const std::vector<T> vec(T & t1, T & t2, T & t3, T & t4, T & t5, T & t6) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	v.push_back(t4);
	v.push_back(t5);
	v.push_back(t6);
	return v;
}

template<typename T>
	const std::vector<T> vec(T & t1, T & t2, T & t3, T & t4, T & t5, T & t6, T & t7) {
	std::vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	v.push_back(t4);
	v.push_back(t5);
	v.push_back(t6);
	v.push_back(t7);
	return v;
}

std::vector<char*> get_tokens(char* line, char sep);
std::vector<MDArray<double> > data_loader(const char* filename, char sep = ' ', std::vector<uint> columns = std::vector<uint>());
std::vector<MDArray<eMISMASK> > toMismask(std::vector<MDArray<double> > data);
void itoa(int n, char s[]);
char *ftoa(float f, int *status);
bool FileExists(std::string strFilename);
}


void test_library();

#endif
