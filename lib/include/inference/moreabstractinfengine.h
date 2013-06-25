/*
 * moreabstractinfengine.h
 *
 *  Copyright (C) 2009, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef MOREABSTRACTINFENGINE_H_
#define MOREABSTRACTINFENGINE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/dbn.h"

namespace mocapy {

// TODO: At some point MoreAbstractInfEngine, SampleInfEngine and AbstractInfEngine should be combined somehow
class MoreAbstractInfEngine {
public:
	MoreAbstractInfEngine() {};
	MoreAbstractInfEngine(DBN * new_dbn);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	void check_seq_and_mismask(Sequence & seq, MDArray<eMISMASK> & mismask);

	std::vector<ParentMap> make_parentmap_list(Sequence * seq, double weight=1.0);
	void set_parentmap(std::vector<ParentMap> & l);
	void set_parentmap(Sequence * seq, double weight=1.0);

	DBN* dbn;
};


template<class Archive>
void MoreAbstractInfEngine::serialize(Archive & ar, const unsigned int version) {
	ar & dbn;
}
}

#endif /* ABSTRACTINFENGINE_H_ */
