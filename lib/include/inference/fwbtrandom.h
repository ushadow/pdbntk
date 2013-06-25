/*
 * fwbtrandom.h
 *
 *  Copyright (C) 2009, Wouter Boomsma, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef FWBTRANDOM_H_
#define FWBTRANDOM_H_

#include <boost/serialization/base_object.hpp>
#include <fstream>

#include "../inference/mcmc.h"

namespace mocapy {

// Forward-backtrack sampler for use in E-step
class FwbtRandom: public MCMC {
public:
	FwbtRandom() {};
	FwbtRandom(DBN* new_dbn);

	std::vector<std::pair<uint, uint> > traverse(uint nr_nodes, uint start, uint end);
	std::pair<double, uint> sweep(uint mcmc_steps, bool burn_in_flag,
                                 Sequence &seq, MDArray<eMISMASK> & mismask, uint start, uint end);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

template<class Archive>
void FwbtRandom::serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::base_object<MCMC>(*this);
}

}

#endif /* FWBTRANDOM_H_ */
