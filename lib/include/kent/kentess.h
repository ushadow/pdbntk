/*
 *  kentess.h
 *
 *  Copyright (C) 2008, Kasper Stovgaard, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef KENTESS_H_
#define KENTESS_H_

#include "../framework/essbase.h"


namespace mocapy {

// Indices for the ess std::vector, containing the average unit std::vector (ybar) and dispersion (S)
enum KENT_ESS_INDEX {K_YBAR, K_S, K_N, K_DATALENGTH, K_ESSSIZE};

class KentESS: public ESSBase {
public:
     KentESS() {};
     virtual ~KentESS() {};

     // Initializes the ESS arrays. Called in Node::construct.
     void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

     // Add another sample to the ESS
     void add_ptv(std::vector<double> ptv);
     
     // Save ESS (called after each sequence)  
     //void save_ess(double weight=1.0);

private:
     uint output_dim;
     uint node_size;
};

}

#endif /*KENTESS_H_*/
