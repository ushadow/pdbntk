/*
 * vonmises.h
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre,
 *  University of Copenhagen.
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

#ifndef VONMISES_H_
#define VONMISES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <math.h>
#include "../utils/utils.h"
#include "../utils/randomgen.h"

namespace mocapy {

class VonMises {
public:
	VonMises() {};
	VonMises(double new_mu, double new_kappa);
	inline double get_lik(double angle, bool log_space=false);
	inline double sample(RandomGen* rg);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	double mu;
	double kappa;

	// Variables for the density function
	double log_denominator;

	// Variables for the sampler
	double a, b, r;
};

inline double VonMises::get_lik(double angle, bool log_space) {
	// Returns the likelihood of angle.
	// angle: argument to the Von Mises density function
	// log_space: if true return log(density)

	double log_nominator = kappa * cos(angle-mu);
	double log_density = log_nominator - log_denominator;

	if(log_space) {
		return log_density;
	}
	else {
		return exp(log_density);
	}
}

inline double VonMises::sample(RandomGen* rg) {
	// Returns a sample from the Von Mises distribution.
	//
	// Based on the algorithm of Best & Fisher (1979), as described in: Jupp PE
	// and Mardia KV,  "Directional Statistics",  John Wiley & Sons, 1999.

	// For kappas close to zero return a sample from the uniform distribution on the circle
	if(kappa < 1e-6) {
		return 2.0 * M_PI * rg->get_rand();  // TODO: This should be a random number in the interval [0,1)
	}

	// For all other use the algorithm of Best & Fisher (1997)
	double U1, U2, U3, z, f, c, theta;

	do {
		U1 = rg->get_rand();
		z = cos(M_PI * U1);
		f = (1.0 + r * z) / (r + z);
		c = kappa * (r - f);

		U2 = rg->get_rand();
	} while ( ( c * (2.0 - c) - U2 <= 0 ) && ( c * exp(1.0 - c) - U2 < 0 ) );

	U3 = rg->get_rand();

	if(U3 - 0.5 > 0) {
		theta = mu + acos(f);
	}
	else {
		theta = mu - acos(f);
	}

	// Make sure that theta is in [0, 2Pi[
	bool less;
	double TWO_M_PI = 2.0 * M_PI;
	while( (less=(theta<0.0)) || (theta>TWO_M_PI)) {
		if (less)
			theta += TWO_M_PI;
		else
			theta -= TWO_M_PI;
	}

	return theta;
}

template<class Archive>
void VonMises::serialize(Archive & ar, const unsigned int version) {
	ar & mu;
	ar & kappa;
	ar & log_denominator;
	ar & a;
	ar & b;
	ar & r;
}

}

#endif /* VONMISES_H_ */
