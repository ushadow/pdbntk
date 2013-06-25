/*
 *  GDdensity.h
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
#ifndef GDDENSITY_H_
#define GDDENSITY_H_

#include <vector>
#include <iostream>
#include <sstream>
#include "../utils/utils.h"
#include <boost/random/mersenne_twister.hpp>
namespace mocapy
{
  class GDdensity{

  public:
    GDdensity();
    GDdensity(uint dim, ulong s, int nn);
    void GDadd_pvec(vector<double> & c);
    void GDset_parameters(vector<double> & as,vector<double> & bs);
    void GDinit_pars();
    void reset_pvec();
    vector<double> GDSample_mult(vector<double> & alphas,vector<double> & betas, int N);
    double GDloglik(vector<double> & countv,bool init);
    double GDsum_loglik();
    vector<double> GDSample_Gendir(vector<double> & alphas, vector<double> & betas);    
    vector< vector<double> > GDestimate_pars(bool init);

  private:
    vector< vector<double> > GDreject(vector< vector<double> > & p);
    double GDsample_beta(double a,double b);
    void GDscanData();
    double lngammaf(double c);
    double GDlngammai(int c);
    double mean(int dim, vector< vector<double> > & vec);
    double var(double mean, int dim,vector< vector<double> > & pvec);
    int place_in_sequence(vector<double> & a, double item);
    bool GDsortdim(vector< vector<double> > & vec);
    double GDcovsum(int d1,double *means,vector< vector<double> > & vec);
    void GDswap_dims(vector< int > & idx,vector< vector<double> > & vec);


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
    void GDdensity::serialize(Archive & ar, const unsigned int version)
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
#endif /* GDDENSITY_H_ */
