/*
 *  DirichletDensity.h
 *
 *  Copyright (C) 2010, Peter Kerpedjiev, The Bioinformatics Centre,
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

#ifndef DIRICHLETDENSITY_H_
#define DIRICHLETDENSITY_H_

#include "../utils/utils.h"
#include "../utils/LogFactorial.h"
#include <boost/random.hpp>

namespace mocapy
{

// Density function for Dirichlet distribution
class DirichletDensity {
public:
	DirichletDensity() {};
	DirichletDensity(MDArray<double> & p, uint sample_size = 1);

	inline double get_lik(std::vector<double> & counts, bool log_space = false,
        uint start_index = 0);

    inline std::vector<double> sample(uint sample_size = 0);

    void set_parameters(MDArray<double> &a);
    MDArray<double> get_parameters() { return a; }

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);

private:
	uint dim;
    uint sample_size;
	MDArray<double> a;
    MDArray<double> am;

    double lnc;

    // Will contain sampled counts
    std::vector<double> sampled;
};

inline double DirichletDensity::get_lik(std::vector<double> & ps,
    bool log_space, uint start_index)
{
	// Return likelihood of counts.
	// counts: argument to multinomial density function
	// log_space: if true return log(density)
    // start_index: where do the counts start in the std::vector
    double ll = lnc;

    bool validlog=true;
    for (uint i = 0; i < dim; i++)
    {
        double p = ps[i + start_index];
        ll += am[i] * log(p);
    }

    // test this before removing
    validlog = true;

    if (!validlog)
    	ll = -INF;

    if (log_space)
        return ll;
    else {
    	if (validlog)
    		return exp(ll);
    	else
			return 0;
    }
}

inline std::vector<double> DirichletDensity::sample(uint sample_size)
{
    if (sample_size == 0)
        sample_size = this->sample_size;

    // clear sample container
    sampled.assign(dim, 0);

    static boost::mt19937 rng(static_cast<unsigned> (std::time(0)));

    for (uint i = 0; i < sample_size; i++) {
        double sum = 0.0; 

        for (uint j = 0; j < dim; j++) {
            boost::gamma_distribution<> gamma(a[j]);
            boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > rvg(rng, gamma);
            sampled[j] = rvg();
            sum += sampled[j];
        }

        for (uint j = 0; j < dim; j++) {
            sampled[j] /= sum;
        }
    }

    return sampled;
}

template<class Archive>
    void DirichletDensity::serialize(Archive & ar, const unsigned int version)
{
    ar & dim;
    ar & sample_size;
	ar & a;
    ar & am;
    ar & lnc;
    ar & sampled;
}

}

#endif /* DIRICHLETDENSITY_H_ */
