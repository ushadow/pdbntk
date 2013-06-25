/*
 * infenginemm.h
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

#ifndef INFENGINEMM_H_
#define INFENGINEMM_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/dbn.h"
#include "moreabstractinfengine.h"
#include "sampleinfengine.h"
#include "../inference/detail/generalinfenginemm.h"
#include "../utils/randomgen.h"

namespace mocapy {

/* Classes for sampling from and calculating likehood for DBNs that are
 * mixture models.
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
 *               |  H  |             |  H  |
 *                \---/	              \---/
 *              / /   \ \           / /   \ \
 *             / /     \ \         / /     \ \
 *            V V       V V       V V       V V
 *            O O  ...  O O       O O  ...  O O
 *
 * where I are input nodes, H is the hidden node and O are output
 * nodes. Further, it is assumed that all inputs are observed and that
 * the hidden node unobserved (no matter the contents of the mismask).
 * The output nodes cannot have any children.
 *
 * The classes must be initialized with the *index* of the hidden
 * node. Per default, the design of the network is check for
 * comparability.
 */

// Sampler for MMs
class SampleInfEngineMM : public SampleInfEngine  {
public:
	SampleInfEngineMM() {};
	SampleInfEngineMM(DBN* new_dbn, Sequence & new_seq, MDArray<eMISMASK> & new_mismask,
			uint new_hidden_node_index, bool check_dbn=true, double new_weight=1.0);

	Sequence & sample_next();
	double calc_ll(MDArray<eMISMASK> & mismask, int start=0, int end=-1, bool multiply_by_input=false);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	GeneralInfEngineMM generalinfengine;
};

template<class Archive>
void SampleInfEngineMM::serialize(Archive & ar, const unsigned int version) {
	ar & generalinfengine;
};


// Likelihood calculator for MMs
class LikelihoodInfEngineMM : public MoreAbstractInfEngine {
public:
	LikelihoodInfEngineMM() {};
	LikelihoodInfEngineMM(DBN* new_dbn, uint new_hidden_node, bool check_dbn=true);

	double calc_ll(Sequence & seq, MDArray<eMISMASK> & mismask, int start=0, int end=-1, bool multiply_by_input=false);

	// Likelihood calculation of entire dataset
    double calc_ll(std::vector<Sequence *> &seqs, std::vector<MDArray<eMISMASK> > &mismasks, bool multiply_by_input=false);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	GeneralInfEngineMM generalinfengine;
};


template<class Archive>
void LikelihoodInfEngineMM::serialize(Archive & ar, const unsigned int version) {
	ar & generalinfengine;
};

}

#endif /* INFENGINEMM_H_ */
