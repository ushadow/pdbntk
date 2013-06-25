/*
 * bippo.h
 *
 *  Copyright (C) 2008, Thomas Hamelryck, The Bioinformatics Centre,
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

#ifndef BIPPO_H_
#define BIPPO_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <math.h>
#include "../utils/utils.h"

#include "../utils/LogFactorial.h"

// These will be used to make sure moment estimates make sense
// Boundaries for theta (binomial)
#define THETA_MAX 0.9
#define THETA_MIN 0.1
// Minimum for lambda (Poisson)
#define LAMBDA_MIN 0.1

namespace mocapy {

class Bippo {
public:
	Bippo() {};
	Bippo(double lambda, double theta, uint n, uint cache_max=50);
	inline double get_lik(uint c, bool log_space=false);
	double sample(RandomGen* rg);
    // estimate AND update parameters
	void estimate(double m, double m2, double m3);
    std::vector<double> get_parameters(void);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
    // METHODS
    void cache(void); // cache loglik/lik
    // Calculate loglik (only used in caching)
    double calc_loglik(uint c); 
    double poisson_sample(RandomGen* rg);

    // These std::vector cache the lik/loglik
    std::vector<double> lik_cache;
    std::vector<double> loglik_cache;
    // Poisson parameter
	double lambda; 
    // Binomial parameters
	double theta;
    uint n;
    // Loglik/lik values will be cached up to cache_max
    uint cache_max;
    // Cached log factorial (use: log_fac.get(uint) 
    LogFactorial log_fac;
};

inline double Bippo::get_lik(uint c, bool log_space) {
    if (c<cache_max)
    {
        // cached
        if (log_space) return loglik_cache[c];
        else return lik_cache[c];
    }
    else
    {
        // not cached - calculate
        double ll = calc_loglik(c);

        if (log_space) return ll;
        else return exp(ll);
    }
}

template<class Archive>
void Bippo::serialize(Archive & ar, const unsigned int version) {
  
    ar & lambda;
	ar & theta;
	ar & n;
	ar & cache_max;
	ar & log_fac;

	// Redo cache on load (and save)
	cache();
}

}

#endif /* BIPPO_H_ */
