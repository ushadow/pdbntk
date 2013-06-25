/*
 * ParentMap.h
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

#ifndef PARENTMAP_H_
#define PARENTMAP_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>

#include "../utils/utils.h"

// Fast parent node value/node value lookup for discrete nodes.

namespace mocapy {

class ParentMap {
public:
	ParentMap() {};
	ParentMap(Sequence * new_seq, std::vector<uint> & parents_0,
			std::vector<uint> & parents_1, uint new_this_index, double new_weight, uint new_dim);
	virtual ~ParentMap() {};
	inline void get(uint l, std::vector<double> & ptv);
	inline void get(uint l, uint & pv, std::vector<double> & v);
	inline void set(uint l, uint v);
	inline void set(uint l, double v);
	void set(uint l, MDArray<double> & v);
	void set(uint l, std::vector<double> & v);
	double get_weight();
	void replace_seq(Sequence *new_seq);
	uint lng;
	uint dim;

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	FlatSequence seq;
	std::vector<std::vector<uint> > indices;
	uint step;
	double weight;
	uint this_index;
};

inline void ParentMap::set(uint i, uint v) {
	std::vector<uint> & v1 = indices[i];
	uint j = v1.back();
	*(seq[j]) = v;
}


inline void ParentMap::set(uint i, double v) {
	std::vector<uint> & v1 = indices[i];
	uint j = v1.back();
	*(seq[j]) = v;
}

inline void ParentMap::get(uint i, std::vector<double> & ptv) {
	assert(i < indices.size());
	std::vector<uint> & index = indices[i];
	// Return the values of seq with indices in index
	take_value(seq, index, ptv);
	uint j = index.back();
	for (uint k = 1; k < dim; k++) {
		ptv.push_back(*((seq)[j + k]));
	}
}

template<class Archive>
void ParentMap::serialize(Archive & ar, const unsigned int version) {
	ar & lng;
	ar & indices;
	ar & step;
	ar & weight;
	ar & this_index;
}

}

#endif /* PARENTMAP_H_ */
