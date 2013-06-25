/*
 *  GDMdensities.h
 *
 *  Copyright (C) 2010, Simon H. Rasmussen, The Bioinformatics Centre,
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
#ifndef GDMDENSITIES_H_
#define GDMDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "../framework/mocapyexceptions.h"

#include "GDMdensity.h"
#include "../utils/LogFactorial.h"

namespace mocapy
{

class GDMdensities;
std::ostream& operator<<(std::ostream&, const GDMdensities&);

// Density function for multidimensional Gaussian distribution
class GDMdensities : public DensitiesBase {
public:
	GDMdensities();
	GDMdensities(uint dim,uint sample_size=1);

	virtual ~GDMdensities(){};

    // Initializes the density arrays
    void construct(vector<uint> & parent_sizes);

    void initialize();

    // Parameter estimation based on the ESS
    void estimate(vector< MDArray<double> > & ess);

    // Return a sample, based on indicated parent values
    vector<double> sample(vector<double> & pv);

    // Return likelihood, that is: P(child|parents)
    double get_lik(vector<double> & ptv, bool log = false);

    void update_pars();

    vector<int> get_idxVec();

    void set_SumOfCounts(int s);
    void set_metro_samples(int m);

    // Return the distribution's parameters
    vector<MDArray<double> > get_parameters();
    void set_parameters(MDArray<double> & cpd);

    // Persistence
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
    friend ostream& operator<< (ostream& output, const GDMdensities& a);

    vector<GDMdensity> densities;

private:
    void _create_densities(void);
    vector< MDArray<double> > cpd;
    double pseudo_count;
    uint dim;
    uint sample_size;
    vector<bool> init;
    bool first;
    int metro_samples;
    int sumOfCounts;
};
/* */
template<class Archive>
    void GDMdensities::serialize(Archive & ar,
        const unsigned int version)
{
  	ar & boost::serialization::base_object<DensitiesBase>(*this);
/*     ar & cpd;
    ar & pseudo_count;
    ar & dim;
    ar & sample_size;
    ar & init;
    ar & densities;
 */  
}
/* */
}

#endif /* GDMDENSITIES_H_ */
