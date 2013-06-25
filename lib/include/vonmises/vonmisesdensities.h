/*
 * vonmisesdensities.h
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef VONMISESDENSITIES_H_
#define VONMISESDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>


#include "../framework/densitiesbase.h"
#include "vonmises.h"

namespace mocapy {

class VonMisesDensities;
std::ostream& operator<<(std::ostream&, const VonMisesDensities&);

class VonMisesDensities : public DensitiesBase {
public:
	VonMisesDensities(std::vector<double> user_mus = std::vector<double> (),
			std::vector<double> user_kappas = std::vector<double> ());

	virtual ~VonMisesDensities() {};

	// Initializes the Density arrays
	void construct(std::vector<uint> & parent_sizes);

	// Parameter estimation based on the ESS
	void estimate(std::vector<MDArray<double> > & ess);

	// Return a sample, based on indicated parent values
	std::vector<double> sample(std::vector<double> & ptv);

	// Return likelihood, that is: P(child|parents)
	double get_lik(std::vector<double> & ptv, bool log=false);

	// Return the distribution's parameters
	std::vector<MDArray<double> > get_parameters();

	// Discrete node specific
	std::vector<double> & get_mus() {return mus;}
	std::vector<double> & get_kappas() {return kappas;}

	friend std::ostream& operator<< (std::ostream& output, const VonMisesDensities& a);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	std::vector<double> make_uniform_mus(uint node_size);
	std::vector<double> make_uniform_kappas(uint node_size);
	void initialize();
	std::pair<std::pair<double, double>, bool> estimate_mu_kappa(double rx, double ry, int n);
	std::vector<VonMises> make_vm_list(uint node_size, std::vector<double> mus, std::vector<double> kappas);

	std::vector<double> user_mus;
	std::vector<double> user_kappas;
	std::vector<double> mus;
	std::vector<double> kappas;
	std::vector<VonMises> vm_list;
	uint parent_index;
};


template<class Archive>
void VonMisesDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & user_mus;
	ar & user_kappas;
	ar & mus;
	ar & kappas;
	ar & vm_list;
	ar & parent_index;
}

}

#endif /* VONMISESDENSITIES_H_ */
