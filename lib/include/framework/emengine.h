/*
 * EMEngine.h
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

#ifndef EMENGINE_H_
#define EMENGINE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "dbn.h"
#include "../inference/infenginemcmc.h"
#include "../inference/mcmc.h"

namespace mocapy {

// Expectation Maximization (EM) engine, using Monte Carlo EM.
class EMEngine {
public:
  EMEngine() {dataIsVerified=false;}

	EMEngine(DBN* new_dbn, MCMC* new_sampler, std::vector<Sequence>* new_seq_list = NULL,
			std::vector< MDArray<eMISMASK> > * new_mismask_list = NULL, std::vector<double> * new_weight_list=NULL);

	void initialize(DBN* new_dbn, MCMC* new_sampler,
			std::vector<Sequence*> * new_seq_list,
			std::vector<MDArray<eMISMASK> > * new_mismask_list, std::vector<double> * new_weight_list);

	void do_E_step(uint mcmc_steps, uint burn_in_steps, bool init_random=false);
	void do_M_step();
	double get_loglik();
	std::vector<InfEngineMCMC> make_inf_list();

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	// File format:
    // A line contains values in a slice separated by "sep" (default space).
    // A block of lines therefore represents a sequence.
    // An empty line separates two sequences.
    // The columns vector holds the numbers of the columns that should be read from the file (default is all)
    // Lines starting with # are ignored
    //
    // Example:
    //---------
    // # This line is ignored
    // 6 5 0 1
    // 4 3 6 1
    // 6 4 1 8
    //
    // 9 8 7 5
    // 2 3 4 0
    // 5 5 5 5
    //---------
    // Reading the lines above, results in two sequences of length three having four values in a slice.
	bool load_sequences(const char* filename, char sep = ' ', std::vector<uint> columns=std::vector<uint>());

	// Mismask files have the same structure as sequence files (with 0 for observed and 1 for hidden nodes)
	bool load_mismask(const char* filename, char sep = ' ', std::vector<uint> columns=std::vector<uint>());

	// A list of weights (must have the same number of lines as the number of sequences)
	bool load_weights(const char* filename);

	// Made public in order to allow users to access data after reading it in via load_sequences
	// Perhaps the load sequences should be moves out of the emengine class?
	std::vector<Sequence*> seq_list;
	std::vector< MDArray<eMISMASK> >  mismask_list;

private:
	void verifyData();

	DBN* dbn;
	MCMC* sampler;
	std::vector<InfEngineMCMC> inf_list;
	std::vector<double> weight_list;
	bool E_done;
	bool keep_sequences;
	double loglik;
	uint slice_count;
	uint nr_seqs;
	bool dataIsVerified;
	
};

template<class Archive>
void EMEngine::serialize(Archive & ar, const unsigned int version) {
	ar & dbn;
	ar & sampler;
	ar & inf_list;
	ar & seq_list;
	ar & mismask_list;
	ar & weight_list;
	ar & E_done;
	ar & keep_sequences;
	ar & loglik;
	ar & slice_count;
	ar & nr_seqs;
}

}

#endif /* EMENGINE_H_ */
