/*
 * randomgen.h
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




#ifndef RANDOMGEN_H_
#define RANDOMGEN_H_

#ifdef __MACH__
#define uint unsigned int
#endif

#include <boost/random.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

namespace mocapy {

    extern uint moc_seed;

    class RandomGen {
    public:
        RandomGen(){
            r_lagged_fibonacci607 = new boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_real<> >(rng_lagged_fibonacci607, unf);
            r = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(rng, unf);
            norm = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(rng, normal);

            mocapy_seed(moc_seed);
        }

        ~RandomGen() {
            delete r;
            delete norm;
            delete r_lagged_fibonacci607;
        }
        double get_rand(const bool use_lagged_fibonacci607=false);
        void mocapy_seed(uint s);
        double get_rand_normal(double mean=0.0, double sigma=1.0, const bool use_lagged_fibonacci607=false);

        // For Twister
        uint moc_seed1;
        // for multivariate normal
        int moc_seed2;
        
        // Random number generators methods
        boost::lagged_fibonacci607 rng_lagged_fibonacci607;
        boost::mt19937 rng; //(moc_seed1);
        // distributions
        boost::uniform_real<> unf; //(0,1);
        boost::normal_distribution<> normal; //(0,1);
        
        // Boost generators
        boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_real<> >* r_lagged_fibonacci607; //(rng, unf);
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> >* r; //(rng, unf);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >* norm; //(rng, normal);


        // Persistence
        template<class Archive>
        void save(Archive & ar, const unsigned int version) const;
        
        template<class Archive>
        void load(Archive & ar, const unsigned int version);
        BOOST_SERIALIZATION_SPLIT_MEMBER()

    };



    template<class Archive>
    void RandomGen::load(Archive & ar, const unsigned int version) {
        ar & moc_seed1;

        /*
        r = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(rng, unf); 
        norm = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(rng, normal); 
        
        mocapy_seed(moc_seed);
        */
    }

    template<class Archive>
    void RandomGen::save(Archive & ar, const unsigned int version) const {
        ar & moc_seed1;

        /*
        r = new boost::variate_generator<boost::mt19937&, boost::uniform_real<> >(rng, unf); 
        norm = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(rng, normal); 
        
        mocapy_seed(moc_seed);
        */
    }

}

#endif
