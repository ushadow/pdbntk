/*
 * MultiGauss.h
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

#ifndef MULTIGAUSS_H_
#define MULTIGAUSS_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../utils/utils.h"

namespace mocapy {

// Density function for multidimensional Gaussian distribution
class MultiGauss {
public:
	MultiGauss() {};
	MultiGauss(MDArray<double> & new_mean, MDArray<double> & new_sigma);
	inline double get_lik(std::vector<double> & a, bool log_space = false);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	double c;
	double c2;
	bool one_dim;
	double sq_det;
	double b;
	double logb;
	uint dim;
	MDArray<double> inv_sigma;
	MDArray<double> sigma;
	MDArray<double> mean;
};

inline double MultiGauss::get_lik(std::vector<double> & a, bool log_space) {
	// Return likelihood of a.
	// a: argument to Gaussian density function
	// log_space: if true return log(density)
	if (one_dim) {
		// 1D Gaussian
		assert(a.size() == 1);
		double d = a.front() - mean[0];
		double x = c - d * d / c2;
		if (log_space)
			return x;
		else
			return exp(x);
	} else {
		// Multidimensional Gaussian
		assert(a.size()> 1);
		std::vector<double> p;
		p.reserve(a.size());
		vec_sub(a, mean.get_values(), p);
		std::vector<double> q = inv_sigma*p;

		double r = -0.5 * dot(p, q);
		if (log_space)
			return logb + r;
		else
			return b * exp(r);
	}
}


template<class Archive>
void MultiGauss::serialize(Archive & ar, const unsigned int version) {
	ar & c;
	
	ar & c2;
	
	ar & one_dim;
	ar & sq_det;
	ar & b;


	
	ar & logb;
	ar & dim;
	ar & inv_sigma;
	ar & sigma;
	ar & mean;

}

}

#endif /* MULTIGAUSS_H_ */
