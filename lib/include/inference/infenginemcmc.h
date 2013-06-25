/*
 * infenginemcmc.h
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

#ifndef INFENGINEMCMC_H_
#define INFENGINEMCMC_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "abstractinfengine.h"
#include "mcmc.h"

namespace mocapy {

class Sample {
public:
	Sample() {};
	Sample(Sequence & new_seq, double new_ll, uint new_slice_count) :
		seq(new_seq), ll(new_ll), slice_count(new_slice_count) {
	}

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	Sequence seq;
	double ll;
	uint slice_count;
};

template<class Archive>
void Sample::serialize(Archive & ar, const unsigned int version) {
	ar & seq;
	ar & ll;
	ar & slice_count;
}


class InfEngineMCMC: public AbstractInfEngine {
public:
	InfEngineMCMC() {};
	InfEngineMCMC(DBN* new_dbn, MCMC* new_sampler, Sequence * new_seq,
			MDArray<eMISMASK> & new_mismask, double new_weight = 1.0);

	virtual ~InfEngineMCMC() {};

	void initialize_sample_generator(uint burn_in_steps = 0, bool init_random =	true,
			MDArray<eMISMASK> mismask = MDArray<eMISMASK> (), uint start = 0,
			int end = -1);
	void initialize_viterbi_generator(uint new_mcmc_steps, uint burn_in_steps = 0,
			bool init_random = false, MDArray<eMISMASK> mismask = MDArray<eMISMASK> (),
			uint start = 0, int end = -1, bool restart = false);

	Sample viterbi_next();
	Sample sample_next();

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

private:
	MCMC* sampler;
	MDArray<eMISMASK> mismask;

	std::vector<ParentMap> parentmap_list;
	bool is_initialized;
	bool burn_in_flag;

	// For the generators
	MDArray<int> g_mismask;
	uint g_start;
	uint g_end;
	Sequence prev_seq;
	uint mcmc_steps;
};

template<class Archive>
void InfEngineMCMC::serialize(Archive & ar, const unsigned int version) {
	ar & boost::serialization::base_object<AbstractInfEngine>(*this);
	ar & sampler;
	ar & mismask;
	ar & parentmap_list;
	ar & is_initialized;
	ar & burn_in_flag;
	ar & g_mismask;
	ar & g_start;
	ar & g_end;
	ar & prev_seq;
	ar & mcmc_steps;
}

}

#endif /* INFENGINEMCMC_H_ */
