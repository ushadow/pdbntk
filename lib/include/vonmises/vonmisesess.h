/*
 * vonmisesess.h
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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


#ifndef VONMISESESS_H_
#define VONMISESESS_H_

#include "../framework/essbase.h"

namespace mocapy {

// Indices for the ess std::vector
enum VONMISES_ESS_INDEX {VM_RX, VM_RY, VM_N, VM_ESSSIZE};

class VonMisesESS: public ESSBase {
public:
	VonMisesESS() {};
	virtual ~VonMisesESS() {};

	// Initializes the ESS arrays. Called in Node::construct.
	void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

	// Add another sample to the ESS
	void add_ptv(std::vector<double> ptv);

private:
	uint output_dim;
	uint node_size;
};

}


#endif /* VONMISESESS_H_ */
