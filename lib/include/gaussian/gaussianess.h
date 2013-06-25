/*
 * GaussianESS.cpp
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

#ifndef GAUSSIANESS_H_
#define GAUSSIANESS_H_

#include "../framework/essbase.h"

namespace mocapy {

// Indices for the ess std::vector
enum GAUSSIAN_ESS_INDEX {G_S, G_SY, G_SYTY, G_STYY, G_DATALENGTH, G_ESSSIZE};

class GaussianESS: public ESSBase {
 public:
  GaussianESS() {useShrinkage=false;}
  virtual ~GaussianESS() {};
  
  // Initializes the ESS arrays. Called in Node::construct.
  void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);
  
  // Add another sample to the ESS
  void add_ptv(std::vector<double> ptv);
  
  void save_ess(double weight);
  void clear();

  bool useShrinkage;

  
 private:
  uint output_dim;
  uint node_size;
};
 
}

#endif /* GAUSSIANESS_H_ */
