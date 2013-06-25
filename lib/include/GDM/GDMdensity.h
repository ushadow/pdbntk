/*
 *  GDMdensity.h
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

#ifndef GDMDENSITY_H_
#define GDMDENSITY_H_

#include <vector>
#include <iostream>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../utils/utils.h"
#include <boost/random/mersenne_twister.hpp>

namespace mocapy
{
  class GDMdensity{

  public:
    GDMdensity();
    GDMdensity(uint dim, ulong s, int nn);
    void add_countv(vector<double> & c);
    void set_parameters(vector<double> & as,vector<double> & bs);
    void init_pars();
    vector<double> Sample_mult(vector<double> & alphas,vector<double> & betas, int N);
    double getMult_loglik(vector<double> & countv, vector<double> & f);
    vector< vector<double> > metro_hast(int samples,bool init,bool sort);
    vector<int> get_idxVec();
    double loglik(vector<double> & countv);
    double sum_loglik();
    void reset_cvec();
  private:
    vector<double> Sample_Gendir(vector<double> & alphas, vector<double> & betas);    
    vector< vector<double> > reject(vector< vector<double> > & p);
    vector< vector<double> > estimate_pars(bool init);
    double sample_beta(double a,double b);
    void scanData();
    double lngammai(int c);
    double mean(int dim, vector< vector<double> > & vec);
    double var(double mean, int dim,vector< vector<double> > & pvec);
    bool sortdim(vector< vector<double> > & vec);
    double covsum(int d1,double *means,vector< vector<double> > & vec);
    void swap_dims(vector< int > & idx,vector< vector<double> > & vec);


    boost::mt19937 rng1;                 
    boost::mt19937 rng2;                 
    vector<int> idxvec;
    vector< vector<double> > pvector;
    vector< vector<double> > cvector;
    vector<double> alphas;
    vector<double> betas;
    int n;
    int Nmax;
    vector<double> lnGammaCache;
    int dim;
    bool cache_initialized;

    // Persistence
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version);
    
  };

  template<class Archive>
    void GDMdensity::serialize(Archive & ar, const unsigned int version)
    {
      ar & idxvec;
      ar & pvector;
      ar & cvector;
      ar & alphas;
      ar & betas;
      ar & n;
      ar & Nmax;
      ar & lnGammaCache;
      ar & dim;
      ar & rng1;
      ar & rng2;
      ar & cache_initialized;
    }
}
#endif /* GDMDENSITY_H_ */
