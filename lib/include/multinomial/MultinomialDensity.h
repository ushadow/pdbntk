/*
 *  multinomialdensity.h
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




#ifndef MULTINOMIALDENSITY_H_
#define MULTINOMIALDENSITY_H_

#include <boost/math/special_functions/gamma.hpp>

#include "../utils/randomgen.h"

#include "../utils/LogFactorial.h"

#include "../utils/utils.h"

namespace mocapy
{

// Density function for multidimensional Gaussian distribution
class MultinomialDensity {
public:
	MultinomialDensity() {};
	MultinomialDensity(MDArray<double> & p, uint sample_size = 1);

	inline double get_lik(std::vector<double> & counts, bool log_space = false,
        uint start_index = 0);

	inline std::vector<double> sample(RandomGen * rg, uint sample_size = 0);

    void set_parameters(MDArray<double> &a);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);

private:
	uint dim;
    uint sample_size;
	MDArray<double> p;
	MDArray<double> log_p;
	MDArray<double> cum_p;
    // Cache of the log factorial calculation
    LogFactorial log_fac;
    // Will contain sampled counts
    std::vector<double> sampled;
};

inline double MultinomialDensity::get_lik(std::vector<double> & counts,
    bool log_space, uint start_index)
{
	// Return likelihood of counts.
	// counts: argument to multinomial density function
	// log_space: if true return log(density)
    // start_index: where do the counts start in the std::vector
    uint sum_counts = 0;
    double ll = 0;


    for (uint i = 0; i < dim; i++)
    {
        double count = counts[i + start_index];
        uint int_count = (unsigned int) (count + 0.1);
        ll -= log_fac.get(int_count);
        ll += count * log_p[i];
        sum_counts+=int_count;
    }

    ll += log_fac.get(sum_counts);

    if (log_space)
        return ll;
    else {
    	return exp(ll);
    }
}

 inline std::vector<double> MultinomialDensity::sample(RandomGen * rg, uint sample_size)
{
    if (sample_size == 0)
        sample_size = this->sample_size;

    // clear sample container
    sampled.assign(dim, 0);

    for (uint i = 0; i < sample_size; i++) {
        double r = rg->get_rand();
        uint index = cum_p.bisect(r);
        sampled[index] += 1;
    }

    return sampled;
}

template<class Archive>
    void MultinomialDensity::serialize(Archive & ar, const unsigned int version)
{
    ar & dim;
    ar & sample_size;
	ar & p;
	ar & log_p;
	ar & cum_p;
    ar & log_fac;
    ar & sampled;
}

}

#endif /* MULTINOMIALDENSITY_H_ */
