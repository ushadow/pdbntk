/*
 * essbase.h
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

#ifndef ESSBASE_H_
#define ESSBASE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../utils/utils.h"

namespace mocapy {

class ESSBase {
public:
	ESSBase() {};
	virtual ~ESSBase() {};

	// Initializes the ESS arrays. Called in Node::construct.
    virtual void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size) = 0;

    // Add another sample to the ESS
    virtual void add_ptv(std::vector<double> ptv) = 0;

    // End of sequence during Gibbs sampling.
    virtual void save_ess(double weight=1.0);

    // ESS has been used in M-step - get ready for new cycle
    virtual void clear();

    // Return a std::vector of MDArrays containing the ESS
    // This will be used in the density object for estimation
    virtual std::vector< MDArray<double> > get_array();

    // Combine ESS - for use with MPI
    virtual void combine(ESSBase & other_ess);

protected:
	std::vector< MDArray<double> > ess;
	std::vector< MDArray<double> > saved_ess;
	std::vector< MDArray<double> > ess_shape;
	uint ess_size;
};

}

#endif /* ESS_H_ */
