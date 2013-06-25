/*
 * mcmc.h
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

#ifndef MCMC_H_
#define MCMC_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/dbn.h"

namespace mocapy {

class MCMC {
public:
	MCMC() {};
	MCMC(DBN * new_dbn);
	virtual ~MCMC() {};

	virtual std::pair<double, uint> sweep(uint mcmc_steps, bool burn_in_flag,
                                         Sequence &seq, MDArray<eMISMASK> & mismask, uint start, uint end)=0;
	void initialize(MDArray<eMISMASK> & mismask, uint start, uint end);

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

protected:
	DBN* dbn;
	uint total_output_size;
};

template<class Archive>
void MCMC::serialize(Archive & ar, const unsigned int version) {
	ar & dbn;
	ar & total_output_size;
}

}

#endif /* MCMC_H_ */
