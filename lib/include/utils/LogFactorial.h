/*
 *  logfactorial.h
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


#ifndef LOGFACTORIAL_H_
#define LOGFACTORIAL_H_

#ifdef __MACH__
#define uint unsigned int
#endif

#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#define LOGFACTORIAL_MAXARG 1000

namespace mocapy
{

class LogFactorial
{
    public:
        LogFactorial();
        double get(uint i);

        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive & ar, const unsigned int version);

    private:
        static std::vector<double> *cache;
};

template<class Archive>
    void LogFactorial::serialize(Archive & ar, const unsigned int version)
{
    ar & cache;
}

}


#endif
