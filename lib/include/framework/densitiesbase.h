/*
 * DensitiesBase.h
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

#ifndef DENSITIESBASE_H_
#define DENSITIESBASE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../utils/randomgen.h"
#include "../discrete/discretepriors.h"
#include "node.h"

namespace mocapy {


// This is the base class for density classes (DiscreteDensities, GaussianDensities, ...)
class DensitiesBase {
public:
	DensitiesBase() { node_size = output_size = 0; };
	virtual ~DensitiesBase() {};

	// Called in node.construct and initializes the density arrays
	virtual void construct(std::vector<uint> & parent_sizes) = 0;

	// Parameter estimation based on the ESS
	virtual void estimate(std::vector<MDArray<double> > & ess) = 0;

	// Return a sample, based on indicated parent values
	virtual std::vector<double> sample(std::vector<double> & ptv) = 0;

	// Return likelihood, that is: P(child|parents)
	virtual double get_lik(std::vector<double> & ptv, bool log=false) = 0;

	// Return the distribtion's parameters
	virtual std::vector< MDArray<double> > get_parameters() = 0;

	virtual void setRandomGen(RandomGen* rg) {randomGen = rg; }

	uint get_output_size() {return output_size; }
	uint get_node_size() { return node_size; }
	eNodeType get_node_type() {return type; }

	RandomGen* randomGen;

	// Persistence
	friend class boost::serialization::access;

	/*
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const;

	template<class Archive>
	void load(Archive & ar, const unsigned int version);
	BOOST_SERIALIZATION_SPLIT_MEMBER()
	*/


	  
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
	  

protected:
	// Discrete, Gaussian, ...
	eNodeType type;

	// Dimension of output std::vector
	uint output_size;
	uint node_size;
};

template<class Archive>
void DensitiesBase::serialize(Archive & ar, const unsigned int version) {
  	ar & type;
	ar & output_size;
	ar & node_size;

	//	ar & randomGen;
}


/*
template<class Archive>
void DensitiesBase::load(Archive & ar, const unsigned int version) {
	ar & type;
	ar & output_size;
	ar & node_size;


	randomGen = new RandomGen();
}

template<class Archive>
void DensitiesBase::save(Archive & ar, const unsigned int version) const {
	ar & type;
	ar & output_size;
	ar & node_size;
}
*/


BOOST_SERIALIZATION_ASSUME_ABSTRACT(DensitiesBase)

}

#endif /* DENSITIESBASE_H_ */
