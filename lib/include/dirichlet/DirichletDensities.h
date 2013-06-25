/*
 *  DirichletDensities.h
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

#ifndef DIRICHLETDENSITIES_H_
#define DIRICHLETDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "../framework/mocapyexceptions.h"

#include "DirichletDensity.h"
#include "DirichletESS.h"
#include "../utils/LogFactorial.h"

namespace mocapy
{

class DirichletDensities;
std::ostream& operator<<(std::ostream&, const DirichletDensities&);

// Density function for multidimensional Gaussian distribution
class DirichletDensities : public DensitiesBase {
public:
	DirichletDensities();
	DirichletDensities(uint dim,
        double pseudo_count=0,
        uint sample_size=1,
        bool init_random=false);

    virtual ~DirichletDensities() {};

    // Initializes the density arrays
    void construct(std::vector<uint> & parent_sizes);

    void initialize();

    // Parameter estimation based on the ESS
    void estimate(std::vector< MDArray<double> > & ess);

    // Return a sample, based on indicated parent values
    std::vector<double> sample(std::vector<double> & pv);

    // Return likelihood, that is: P(child|parents)
    double get_lik(std::vector<double> & ptv, bool log = false);

    // Return the distribution's parameters
    std::vector<MDArray<double> > get_parameters();
    void set_parameters(MDArray<double> & cpd);

    // Persistence
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    friend std::ostream& operator<< (std::ostream& output, const DirichletDensities& a);

    std::vector<DirichletDensity> densities;

private:
    void _create_densities(void);

    //maybe this should be called cad (conditional alpha distribution)
    MDArray<double> cpd;
    double pseudo_count;
	uint dim;
    uint sample_size;
    bool init_random;

};

template<class Archive>
    void DirichletDensities::serialize(Archive & ar,
        const unsigned int version)
{
	ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & cpd;
    ar & pseudo_count;
    ar & dim;
    ar & sample_size;
    ar & init_random;
    ar & densities;
}

}

#endif /* DIRICHLETDENSITIES_H_ */
