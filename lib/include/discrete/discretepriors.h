/*
 * discretepriors.h
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


#ifndef DISCRETEPRIORS_H_
#define DISCRETEPRIORS_H_

#include <boost/serialization/version.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "discreteess.h"
#include "../framework/essbase.h"

namespace mocapy {

class Prior {
public:
	Prior() {};
	virtual ~Prior() {};
	virtual void apply_prior(MDArray<double> & ess) = 0;

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};




template<class Archive>
void Prior::serialize(Archive & ar, const unsigned int version) {}

class PseudoCountPrior : public Prior {
public:
	PseudoCountPrior() {};
	PseudoCountPrior(MDArray<double> new_pcounts) : pcounts(new_pcounts) {};
	~PseudoCountPrior() {};
	void apply_prior(MDArray<double> & ess);
	Prior* base() { return dynamic_cast<Prior*>(this); }

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	MDArray<double> pcounts;
};

template<class Archive>
void PseudoCountPrior::serialize(Archive & ar, const unsigned int version) {
	if (version > 0) {
		ar & boost::serialization::base_object<Prior>(*this);
		ar & pcounts;
	}

}

}

BOOST_CLASS_VERSION(mocapy::PseudoCountPrior, 1)

#endif /* DISCRETEPRIORS_H_ */
