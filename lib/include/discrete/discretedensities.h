/*
 * discretedensities.h
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

#ifndef DISCRETEDENSITIES_H_
#define DISCRETEDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "discreteess.h"

namespace mocapy {

class DiscreteDensities;
std::ostream& operator<<(std::ostream&, const DiscreteDensities&);

class DiscreteDensities : public DensitiesBase {
public:
	DiscreteDensities();
	DiscreteDensities(uint new_node_size, Prior * new_prior=NULL, bool new_init_random = false);
	DiscreteDensities(uint new_node_size, CPD & new_user_cpd, Prior * new_prior=NULL);
	virtual ~DiscreteDensities() {};

	// Initializes the Density arrays
	void construct(std::vector<uint> & parent_sizes);

	// Parameter estimation based on the ESS
	void estimate(std::vector<MDArray<double> > & ess);

	// Return a sample, based on indicated parent values
	std::vector<double> sample(std::vector<double> & ptv);

	// Return likelihood, that is: P(child|parents)
	double get_lik(std::vector<double> & ptv, bool log=false);

	// Return the distribution's parameters
	std::vector< MDArray<double> > get_parameters();

	void set_user_cpd(CPD & userCPD);
	void set_prior(Prior * new_prior);
	CPD & getCPD() {return cpd;}

	// Persistence
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	friend std::ostream& operator<< (std::ostream& output, const DiscreteDensities& a);

private:
	void test_cpd(CPD & cpd, std::vector<uint> & shape);
	void set_cpd(CPD & cpd);
	void initialize();

	CPD make_uniform_cpd(const std::vector<uint> & shape);
	CPD make_random_cpd(std::vector<uint> & shape, bool no_zeroes = true);
	Prior* prior;
	bool init_random;

	// CPDs
	CPD cpd;
	CPD user_cpd;
	CPD cum_cpd;
	CPD log_cpd;
	std::vector<uint> CPD_shape;
};

template<class Archive>
void DiscreteDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & prior;
    ar & init_random;
    ar & cpd;
    ar & user_cpd;
    ar & cum_cpd;
    ar & log_cpd;
    ar & CPD_shape;
}

}

#endif /* DISCRETEDENSITIES_H_ */
