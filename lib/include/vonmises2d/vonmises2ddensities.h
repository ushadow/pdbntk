/*
 * vonmises2ddensities.h
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
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


#ifndef VONMISES2DDENSITIES_H_
#define VONMISES2DDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "vonmises2d.h"

namespace mocapy {

class VonMises2dDensities;
std::ostream& operator<<(std::ostream&, const VonMises2dDensities&);

class VonMises2dDensities:public DensitiesBase {
private:
     void initialize();

     MDArray<double> make_uniform_kappas(uint node_size);
     MDArray<double> make_uniform_mus(uint node_size);

     // attributes
     MDArray<double> mus;
     MDArray<double> kappas;
     MDArray<double> user_mus;
     MDArray<double> user_kappas;

     uint kappa_size;
     uint parent_index;

     std::vector<VonMises2d> dist_list;


     std::vector<VonMises2d> make_dist_list(uint node_size, MDArray<double> &mus, MDArray<double> &kappas);

public:

     // Constructors
     VonMises2dDensities(MDArray<double> user_mus = MDArray<double> (),
			 MDArray<double> user_kappas = MDArray<double> ());


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

     virtual void setRandomGen(RandomGen* rg);

     // Persistence
     friend class boost::serialization::access;
     template<class Archive>
     void serialize(Archive & ar, const unsigned int version);

     friend std::ostream& operator<< (std::ostream& output, const VonMises2dDensities& a);

};


template<class Archive>
void VonMises2dDensities::serialize(Archive & ar, const unsigned int version) {
	ar & boost::serialization::base_object<DensitiesBase>(*this);

	ar & user_mus;
	ar & user_kappas;
	ar & mus;
	ar & kappas;
	ar & kappa_size;

	ar & parent_index;
	ar & dist_list;
}

}



#endif /* VONMISES2DDENSITIES_H_ */
