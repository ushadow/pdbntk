/*
 * bippodensities.h
 *
 *  Copyright (C) 2009, Thomas Hamelryck, The Bioinformatics Centre, 
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

#ifndef BIPPODENSITIES_H_
#define BIPPODENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"

#include "bippo.h"
#include "bippoess.h"

namespace mocapy {

class BippoDensities;
std::ostream& operator<<(std::ostream&, const BippoDensities&);

class BippoDensities : public DensitiesBase {
public:
	BippoDensities(std::vector<double> new_user_lambdas = std::vector<double> (),
        std::vector<double> new_user_thetas = std::vector<double> (),
        std::vector<double> new_user_ns = std::vector<double> ());
	virtual ~BippoDensities() {};

	// Initializes the Density arrays
	void construct(std::vector<uint> & parent_sizes);

	// Parameter estimation based on the ESS
	void estimate(std::vector<MDArray<double> > & ess);

	// Return a sample, based on indicated parent values
	std::vector<double> sample(std::vector<double> & pv);

	// Return likelihood, that is: P(child|parents)
	double get_lik(std::vector<double> & ptv, bool log=false);

	// Return the distribution's parameters
	std::vector< MDArray<double> > get_parameters();

	friend std::ostream& operator<< (std::ostream& output, const BippoDensities& a);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
    // Bippo parameters
    std::vector<double> user_lambdas;
    std::vector<double> lambdas;
    std::vector<double> user_thetas;
    std::vector<double> thetas;
    std::vector<double> user_ns;
    std::vector<double> ns;

    // List of densities: one for each parent state
    std::vector<Bippo> densities;

    void initialize();
};

template<class Archive>
void BippoDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & user_lambdas;
    ar & lambdas;
    ar & user_thetas;
    ar & thetas;
    ar & user_ns;
    ar & ns;
    ar & densities;
}

}

#endif /* BIPPODENSITIES_H_ */
