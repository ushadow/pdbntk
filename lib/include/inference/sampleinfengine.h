/*
 * sampleinfengine.h
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

#ifndef SAMPLEINFENGINE_H_
#define SAMPLEINFENGINE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "moreabstractinfengine.h"
#include "../framework/node.h"
#include "../framework/dbn.h"
#include "../utils/randomgen.h"

namespace mocapy {

class SampleInfEngine: public MoreAbstractInfEngine  {
public:
	SampleInfEngine() {};
	SampleInfEngine(DBN* new_dbn, Sequence & new_seq, MDArray<eMISMASK> & new_mismask, double new_weight = 1.0);

	void set_start_end(int new_start=0, int new_end=-1);
	void set_seq_mismask(Sequence & new_seq, MDArray<eMISMASK> & new_mismask);
	void undo();


	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
   
	void store_seq();
	void replace_seq_in_paremtmap_list(Sequence & seq);

	RandomGen* randomGen;

	Sequence seq;
	Sequence prev_seq;
	MDArray<eMISMASK> mismask;
	double weight;

	std::vector<ParentMap> parentmap_list;

	uint seq_len;
	uint start;
	uint end;
};


template<class Archive>
void SampleInfEngine::serialize(Archive & ar, const unsigned int version) {
	ar & seq;
	ar & prev_seq;
	ar & mismask;
	ar & weight;
	ar & parentmap_list;
	ar & seq_len;
	ar & start;
	ar & end;
};

inline void SampleInfEngine::set_start_end(int new_start, int new_end) {
	// Set and check start
	start = new_start;
	if (!(0 <= new_start && start < seq_len))     throw MocapyExceptions("SampleInfEngine: start must be between in [0, sequence length)");;

	// Set and check end
	end = new_end<0 ? seq_len : new_end;
	if (!(start < end && end <= seq_len))         throw MocapyExceptions("SampleInfEngine: end must be between in (start, sequence length]");;
}

inline void SampleInfEngine::store_seq() {
	prev_seq = seq;
}

inline void SampleInfEngine::replace_seq_in_paremtmap_list(Sequence & seq) {
	// Replace the current sequence with the previous in the parent map
	for (std::vector<ParentMap>::iterator pm=parentmap_list.begin(); pm<parentmap_list.end(); pm++) {
		pm->replace_seq(& seq);
	}
}
}

#endif /* SAMPLEINFENGINE_H_ */
