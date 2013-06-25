/*
 * DirichletESS.h
 *
 *  Copyright (C) 2010, Peter Kerpedjiev, The Bioinformatics Centre,
 *  University of Copenhagen.
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

#ifndef DIRICHLETESS_H_
#define DIRICHLETESS_H_

#include <cmath>
#include "../framework/essbase.h"

namespace mocapy
{

enum DIRICHLET_ESS_INDEX {G_SUMOFDATA, G_SUMOFLOG, G_DATAPOINTS, G_DIRESSSIZE};

class DirichletESS : public ESSBase {
    public:
        DirichletESS() {};
        virtual ~DirichletESS() {};

        void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

        void add_ptv(std::vector<double> ptv);
    private:
        uint dim;
};

}

#endif // DIRICHLETLESS_H_
