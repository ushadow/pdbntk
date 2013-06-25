/*
 * abstractinfengine.h
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

#ifndef ABSTRACTINFENGINE_H_
#define ABSTRACTINFENGINE_H_


#include "../framework/dbn.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>


namespace mocapy {

class AbstractInfEngine {
public:
	AbstractInfEngine() {};
	AbstractInfEngine(DBN * new_dbn, Sequence * new_seq, double new_weight);

	std::vector<ParentMap> make_parentmap_list(Sequence * seq, double weight);
	void set_parentmap(std::vector<ParentMap> & l);
	std::pair<Sequence, double> get_viterbi();

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	DBN* dbn;
	double weight;
	Sequence* seq;
	uint seq_len;
	uint output_size;
	uint nr_slices;
	std::vector<Sequence> seq_list;
};


template<class Archive>
void AbstractInfEngine::serialize(Archive & ar, const unsigned int version) {
	ar & dbn;
	ar & weight;
	ar & seq;
	ar & seq_len;
	ar & output_size;
	ar & nr_slices;
	ar & seq_list;
}
}

#endif /* ABSTRACTINFENGINE_H_ */
