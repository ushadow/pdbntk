/*
 * infenginehmm.h
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

#ifndef INFENGINEHMM_H_
#define INFENGINEHMM_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/dbn.h"
#include "../inference/moreabstractinfengine.h"
#include "../inference/sampleinfengine.h"
#include "../inference/detail/forwardbacktracker.h"

namespace mocapy {

/* Interfaces to the implementation of the forward-backtrack algorithm
 * for DBNs that are hidden Markov models (implemented in
 * forwardbacktracker.h).
 *
 * The implementation supports models with both multiple input and
 * output for the hidden node. In the implementation it is assumed
 * that the DBN has the design:
 *
 *           I  I  ...  I  I     I  I  ...  I  I
 *            \  \     /  /       \  \     /  /
 *             \  \   /  /         \  \   /  /
 *              V  V V  V           V  V V  V
 *                /---\               /---\
 *               |  H  |------------>|  H  |
 *                \---/	              \---/
 *              / /   \ \           / /   \ \
 *             / /     \ \         / /     \ \
 *            V V       V V       V V       V V
 *            O O  ...  O O       O O  ...  O O
 *
 * where I are input nodes, H is the hidden node and O are output
 * nodes. Further, it is assumed that all inputs are observed and that
 * the hidden node is unobserved (no matter the contents of the mismask).
 * The output nodes cannot have any children.
 *
 * The classes must be initialized with the *index* of the hidden
 * node. Per default, the design of the network is check for
 * comparability.
 */

// Sampler for HMMs
class SampleInfEngineHMM : public SampleInfEngine  {
public:
	SampleInfEngineHMM() {};
	SampleInfEngineHMM(DBN* new_dbn, Sequence & new_seq, MDArray<eMISMASK> & new_mismask,
			uint new_hidden_node_index, bool check_dbn=true, double new_weight=1.0);

	Sequence & sample_next();
	double calc_ll(MDArray<eMISMASK> & mismask, int start=0, int end=-1, bool multiply_by_input=false);

	// TODO: Implement!
	//double calc_ll_of_output(int start=0, int end=-1, bool multiply_by_parents=false);

     	friend class FwbtRandom;

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	ForwardBacktracker fwbt;
};

template<class Archive>
void SampleInfEngineHMM::serialize(Archive & ar, const unsigned int version) {
	ar & fwbt;
};


// Likelihood calculator for HMMs
class LikelihoodInfEngineHMM : public MoreAbstractInfEngine {
public:
	LikelihoodInfEngineHMM() {};
	LikelihoodInfEngineHMM(DBN* new_dbn, uint new_hidden_node, bool check_dbn=true);

	double calc_ll(Sequence & seq, MDArray<eMISMASK> & mismask, int start=0, int end=-1, bool multiply_by_input=false);

	// Likelihood calculation of entire dataset
    double calc_ll(std::vector<Sequence *> &seqs, std::vector<MDArray<eMISMASK> > &mismasks, bool multiply_by_input=false);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	ForwardBacktracker fwbt;
};


template<class Archive>
void LikelihoodInfEngineHMM::serialize(Archive & ar, const unsigned int version) {
	ar & fwbt;
};

}

#endif /* INFENGINEHMM_H_ */
