/*
 * forwardbacktracker.h
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

#ifndef FORWARDBACKTRACK_H_
#define FORWARDBACKTRACK_H_

#include "../../framework/dbn.h"
#include "../../utils/mdarray.h"
#include "../../utils/utils.h"
#include "../../utils/randomgen.h"
#include "../../framework/mocapyexceptions.h"

namespace mocapy {

/* Implementation of the forward-backtrack algorithm for DBNs that are
 * hidden Markov models. The implementation supports models with both
 * multiple input and output for the hidden node.
 *
 * NOTE: This is the raw implementation of the algorithms. It is
 * recommended to use the interfaces to the algorithm in infenginehmm.h
 */

class ForwardBacktracker {
public:
	ForwardBacktracker() {};
	ForwardBacktracker(DBN* new_dbn, uint new_hidden_node, bool check_dbn);

	void check_dbn();
	void forward_pass(uint start, uint end, uint seq_len,
			MDArray<eMISMASK> & mismask, bool multiply_by_parents);
	void backtrack_pass(uint start, uint end, uint seq_len,
			    MDArray<eMISMASK> & mismask, RandomGen* rg);
	double calc_ll();

	DBN * dbn;
	uint hidden_node;
	DiscreteNode *hd_0, *hd_1;
	MDArray<double> *cpd_0, *cpd_1, cpd_1_swaped;
	std::vector<Node*> nodes_0, nodes_1;
	uint hidden_node_size;
	MDArray<double> forward;
	std::vector<double> scales;
	std::vector<MDArray<double>*> transition_matrices;

	// Persistence
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
};

template<class Archive>
void ForwardBacktracker::serialize(Archive & ar, const unsigned int version) {
	ar & dbn;
	ar & hidden_node;
	ar & hd_0;
	ar & hd_1;
	ar & cpd_0;
	ar & cpd_1;
	ar & cpd_1_swaped;
	ar & nodes_0;
	ar & nodes_1;
	ar & hidden_node_size;
	ar & forward;
	ar & scales;
	ar & transition_matrices;
};

}

#endif /* FORWARDBACKTRACK_H_ */
