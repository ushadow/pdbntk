/*
 * vonmises2dess.h --- Bivariate von Mises distribution, cosine variant
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
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


#ifndef VONMISES2DESS_H_
#define VONMISES2DESS_H_

#include "../framework/essbase.h"

namespace mocapy {

class VonMises2dESS: public ESSBase {
public:

     // Indices for the ess std::vector
     enum index {_2COS1=0, _2SIN1, _2COS2, _2SIN2,
		 COS1MIN2, SIN1MIN2, COS1PLUS2, SIN1PLUS2,
		 C1, S1, C2, S2, COUNT, ESSSIZE};

     VonMises2dESS() {};
     virtual ~VonMises2dESS() {};

     // Initializes the ESS arrays. Called in Node::construct.
     void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

     // Add another sample to the ESS
     void add_ptv(std::vector<double> ptv);

     // Save ESS (called after each sequence)
     void save_ess(double weight=1.0);

private:
     uint output_dim;
     uint node_size;
};

}


#endif /* VONMISES2DESS_H_ */
