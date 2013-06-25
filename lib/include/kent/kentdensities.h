/*
 *  kentdensities.h
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


#ifndef KENTDENSITIES_H_
#define KENTDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "kentess.h"



namespace mocapy {

class KentDensities;
std::ostream& operator<<(std::ostream&, const KentDensities&);

class KentDensities : public DensitiesBase {
public:
     KentDensities( std::vector<double> new_user_kappas = std::vector<double> (),
                    std::vector<double> new_user_betas = std::vector<double> (),
                    std::vector<MDArray<double> > new_user_axes = std::vector<MDArray<double> > (),
                    double new_kappa_max=5000,
                    double new_beta_max=100,
                    double new_beta_factor=2.0001,
                    double new_kappa_init=50,
                    double new_beta_init=30 );

     virtual ~KentDensities() {};

     // Initializes the Density arrays
     void construct(std::vector<uint> & parent_sizes);

     // Parameter estimation based on the ESS
     void estimate(std::vector<MDArray<double> > & ess);

     // Return a sample, based on indicated parent values
     std::vector<double> sample(std::vector<double> & ptv);

     // Return likelihood, that is: P(child|parents)
     double get_lik(std::vector<double> & ptv, bool log=true);

     // Return the distribution's parameters
     std::vector<MDArray<double> > get_parameters();

     // Return string representation of parameters
     std::string get_string() const;

     MDArray<double> sphere_rand(uint dim=3);
     void create_rand_axes(MDArray<double> &a);


     // Persistence
     friend class boost::serialization::access;
     template<class Archive>
     void serialize(Archive & ar, const unsigned int version);

     friend std::ostream& operator<< (std::ostream& output, const KentDensities& a);

private:

     std::vector<double> cart_to_polar(MDArray<double> &v);
     void polar_to_cart(double theta, double tau, MDArray<double> &v, uint dim);
     void initialize();

     std::vector<double> user_kappas;
     std::vector<double> user_betas;
     std::vector<MDArray<double> > user_es;
     std::vector<double> kappas;
     std::vector<double> betas;
     std::vector<MDArray<double> > es;

     double kappa_max, beta_max, beta_factor, kappa_init, beta_init;

     // Variables for the density function
     std::vector<double> logc;

     // Variables for the sampler
     std::vector<double> c2, gamma, lam1, lam2, a, b;
     //std::vector<MDArray<double>* > rho;
     std::vector<MDArray<double> > rho;

};


template<class Archive>
void KentDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);

    ar & user_kappas;
    ar & user_betas;
    ar & user_es;
    ar & kappas;
    ar & betas;
    ar & es;

    ar & kappa_max;
    ar & beta_max;
    ar & beta_factor;
    ar & kappa_init;
    ar & beta_init;
    ar & logc;
    ar & c2;
    ar & gamma;
    ar & lam1;
    ar & lam2;
    ar & a;
    ar & b;
    ar & rho;
}

}


#endif /*KENTDENSITIES_H_*/
