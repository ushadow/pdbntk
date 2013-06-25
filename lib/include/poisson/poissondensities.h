/*
 * poissondensities.h
 *
 *  Copyright (C) 2008, Mikael Borg, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef POISSONDENSITIES_H_
#define POISSONDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "poissoness.h"
#include "poissonsample.h"

namespace mocapy {

class PoissonDensities;
std::ostream& operator<<(std::ostream&, const PoissonDensities&);

class PoissonDensities : public DensitiesBase {
public:
	PoissonDensities(std::vector<double> new_user_means = std::vector<double> ());
	virtual ~PoissonDensities() {};

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

	friend std::ostream& operator<< (std::ostream& output, const PoissonDensities& a);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:

        std::vector<double> user_means;
        std::vector<double> means;

	void initialize();


};

template<class Archive>
void PoissonDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & user_means;
    ar & means;
}

}

#endif /* POISSONDENSITIES_H_ */
