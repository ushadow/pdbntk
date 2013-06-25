/*
 * multinomialess.h
 *
 *  Copyright (C) 2008, Thomas Hamelryck, The Bioinformatics Centre,
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

#ifndef MULTINOMIALESS_H_
#define MULTINOMIALESS_H_

#include "../framework/essbase.h"

namespace mocapy
{

class MultinomialESS : public ESSBase {
    public:
        MultinomialESS() {};
        virtual ~MultinomialESS() {};

        void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

        void add_ptv(std::vector<double> ptv);
    private:
        uint dim;
};

}

#endif // MULTINOMIALESS_H_
