/*
 * vonmises2d.h --- Bivariate von Mises distribution, cosine variant
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef VONMISES2D_H_
#define VONMISES2D_H_

#include <math.h>
#include "../utils/utils.h"
#include "../utils/vector_nD.h"
#include "../utils/optimize.h"
#include "../utils/mdarray.h"

namespace mocapy {

     class VM2CosMarginalSampler;

     double computeLogNormConst(Vector_nD &k, double lower=-M_PI, double upper=M_PI);


     // Mixture of univariate von Mises. Used as comparison function in rejection sampling.
     class VM1CosMixtureDensity {
        public:
             Vector_nD k;
             Vector_nD mu;
             Vector_nD proportion;
             Vector_nD cumProportion;

             VM1CosMixtureDensity() {}
             VM1CosMixtureDensity(Vector_nD k, Vector_nD mu, Vector_nD proportion);
             double density(double phi, double scaleTerm=0.0);

             double logDensity(double phi);
             Vector_nD logDensity(Vector_nD phi);

       	  // Persistence
         	  friend class boost::serialization::access;
         	  template<class Archive>
         	  void serialize(Archive & ar, const unsigned int version);

        };

  	  template<class Archive>
  	  void VM1CosMixtureDensity::serialize(Archive & ar, const unsigned int version) {

  		  ar & k;
  		  ar & mu;
  		  ar & proportion;
  		  ar & cumProportion;

  	  }




        // Mixture of univariate von Mises. Used as comparison function in rejection sampling..
        class VM1CosMixtureSampler {
        private:
             int n;                      // Number of mixture components
             Vector_nD k;
             Vector_nD mu;
             Vector_nD proportion;
             Vector_nD cumProportion;
	     RandomGen* rg;

        public:
             // Constructor
             VM1CosMixtureSampler() {}
             VM1CosMixtureSampler(Vector_nD &k, Vector_nD &mu, Vector_nD &proportion, RandomGen* rg);

             // Set random number generator
             void setRandomGen(RandomGen *rg) {
                  this->rg = rg;
             }

             // Return single sample
             double operator()();

             // Return n samples
             Vector_nD operator()(int n);

             // Persistence
             friend class boost::serialization::access;
             template<class Archive>
             void serialize(Archive & ar, const unsigned int version);
        };

   	  template<class Archive>
   	  void VM1CosMixtureSampler::serialize(Archive & ar, const unsigned int version) {

    		  ar & k;
    		  ar & mu;
    		  ar & proportion;
    		  ar & cumProportion;

    	  }



     // Density class
     class VM2cosDensity {
     private:
     public:
	  Vector_nD k;
	  Vector_nD mu;
	  double logNormConst;

	  VM2cosDensity() {};
	  VM2cosDensity(Vector_nD &k, Vector_nD &mu);
	  VM2cosDensity(Vector_nD &k, Vector_nD &mu, double logNormConst);
	  ~VM2cosDensity();
	  void init(Vector_nD &k, Vector_nD &mu);
	  void update(Vector_nD &k, Vector_nD &mu);
	  void update(Vector_nD &k, Vector_nD &mu, double logNormConst);
	  double get_log_likelihood(double *anglePair);
	  double get_likelihood(double *anglePair);

	  // Persistence
  	  friend class boost::serialization::access;
  	  template<class Archive>
  	  void serialize(Archive & ar, const unsigned int version);

     };

	  template<class Archive>
	  void VM2cosDensity::serialize(Archive & ar, const unsigned int version) {

		  ar & k;
		  ar & mu;
		  ar & logNormConst;

	  }


     // Sampler class
     class VM2cosSampler {
     private:
	  Vector_nD k;
	  Vector_nD mu;
	  double logNormConst;
	  VM2CosMarginalSampler *marginalPhiSampler;
	 
     public:
    	 VM2cosSampler() {}
	 VM2cosSampler(Vector_nD &k, Vector_nD &mu, RandomGen* rg  );
	 VM2cosSampler(Vector_nD &k, Vector_nD &mu, double logNormConst , RandomGen* rg);
	  VM2cosSampler(const VM2cosSampler &sampler, RandomGen* rg);
	  ~VM2cosSampler();
	  void update(Vector_nD &k, Vector_nD &mu);
	  void update(Vector_nD &k, Vector_nD &mu, double logNormConst);
	  void sample(std::vector<double> &rVal);
	  RandomGen* rg;

          // Set random number generator
          void setRandomGen(RandomGen *rg);

	  // Persistence
  	  friend class boost::serialization::access;
  	  template<class Archive>
  	  void serialize(Archive & ar, const unsigned int version);
     };

	  template<class Archive>
 	  void VM2cosSampler::serialize(Archive & ar, const unsigned int version) {
		  ar & k;
		  ar & mu;
		  ar & logNormConst;
		  ar & marginalPhiSampler;
	  }


     // VonMises2d component (sampling and density functionality gathered in single object
     class VonMises2d {
     public:

	  VM2cosDensity *density;
	  VM2cosSampler *sampler;
	  RandomGen* rg;

	  // Constructor
	  VonMises2d():density(NULL), sampler(NULL) {}

	  // Constructor
	  VonMises2d(MDArray<double> &k, MDArray<double> &mu, RandomGen* rg) {
	       Vector_nD k_vec(k);
	       Vector_nD mu_vec(mu);
	       this->rg = rg;
	       density = new VM2cosDensity(k_vec, mu_vec);
	       sampler = new VM2cosSampler(k_vec, mu_vec, density->logNormConst, rg);

	  }


	  // Copy constructor
	  VonMises2d(const VonMises2d &other):density(NULL), sampler(NULL) {
	    this->rg = other.rg;
	       if (other.density)
		    density = new VM2cosDensity(*other.density);
	       if (other.sampler)
		    sampler = new VM2cosSampler(*other.sampler);
	  }


	  // Destructor
	  ~VonMises2d() {
	       delete density;
	       delete sampler;
	  }


          // Set random number generator
          void setRandomGen(RandomGen *rg) {
               this->rg = rg;
               if (sampler)
                    sampler->setRandomGen(rg);
          }

	  // Assignment operator
	  VonMises2d &operator=(const VonMises2d &other) {
	       if (other.density)
		    density = new VM2cosDensity(*other.density);
	       else
		    density = NULL;

	       if (other.sampler)
		    sampler = new VM2cosSampler(*other.sampler);
	       else
		    sampler = NULL;
	       return *this;
	  }

	  // Persistence
  	  friend class boost::serialization::access;
  	  template<class Archive>
  	  void serialize(Archive & ar, const unsigned int version);

	  // Update k and mu parameters
	  void update(MDArray<double> &k, MDArray<double> &mu) {
	       Vector_nD k_vec(k);
	       Vector_nD mu_vec(mu);
	       density->update(k_vec, mu_vec);
	       sampler->update(k_vec, mu_vec, density->logNormConst);
	  }


     };

	  template<class Archive>
	  void VonMises2d::serialize(Archive & ar, const unsigned int version) {
		  ar & sampler;
		  ar & density;
	  }


     // Estimator class (Maximum likelihood)
     std::vector<double> VonMises2dEstimate_ML(MDArray<double> &ess);

     // Estimator class (Moment estimation)
     std::vector<double> VonMises2dEstimate_ME(MDArray<double> &ess);

     // // Estimator class (Maximum likelihood)
     // class VM2cosEstimator_ML {
     // public:
     // 	  Vector_nD getParameters(MDArray<double> &ess);
     // };

     // // Estimator class (Moment estimation)
     // class VM2cosEstimator_ME {
     // public:
     // 	  Vector_nD getParameters(MDArray<double> &ess);
     // };

     // Bivariate von Mises - cosine model. Marginal of phi. Sampler
     class VM2CosMarginalSampler {
     private:
          Vector_nD k;
          double logNormConst;
          double root;
          Brent<double> brentOptimizer;

          VM1CosMixtureSampler *envelopeSampler;
          VM1CosMixtureDensity *envelopeDensity;
          double K;
	  RandomGen* rg;

     public:

          // Initializer
          void init(Vector_nD &k);

          // Constructor
          VM2CosMarginalSampler() {}
          // Constructor
          VM2CosMarginalSampler(Vector_nD &k, RandomGen* rg);

          // Constructor
          VM2CosMarginalSampler(Vector_nD &k, double logNormConst, RandomGen* rg);

          // Copy constructor
          VM2CosMarginalSampler(const VM2CosMarginalSampler &sampler, RandomGen* rg);

          // destructor
          ~VM2CosMarginalSampler();
          // Generate a random sample from the marginal distribution.
          // Rejection sampling: 1) Sample random value r from comparison distribution
          //                     2) Sample random value u from uniform distribution
          //                     3) Accept r if u < f_target(r)/f_comparison(r) """

          // Set random number generator
          void setRandomGen(RandomGen *rg) {
               this->rg = rg;
               if (envelopeSampler)
                    envelopeSampler->setRandomGen(rg);
          }

          double operator()(bool *status);
         // Generate n samples
         Vector_nD operator()(int n);

          // Factor b from derivative of density. Defined in proof of theorem 4
          double derivativeFactorB(double phi);

      	// Persistence
          friend class boost::serialization::access;
          template<class Archive>
          void serialize(Archive & ar, const unsigned int version);

          // Find parameter k and scale factor K so the the target marginal density is at
          // no point above the comparison function.
          // Comparison function: Mixture of two univariate von Mises distributions
          class MaxRatio {
          private:
               Vector_nD k;
               double root,logNormConst;
          public:
     	  // Constructor
               MaxRatio(Vector_nD &k, double root, double logNormConst);

     	  // Reparameterize k parameters
               void reparameterisation(double *k);

     	  // Calculate parameters
               double operator()(Vector_nD x);
     	  // Calculate parameters
               double operator()(double x);
          };


          // Find parameter k and scale factor K so the the target marginal density is at
          // no point above the comparison function.
          // Comparison function: Mixture of four univariate von Mises distributions
          // This method is used for large values of k."""
          class MaxRatioLargerMixture {
          private:
               Vector_nD k;
               double root,logNormConst;
          public:
               MaxRatioLargerMixture(Vector_nD &k, double root, double logNormConst);

               // Reparameterisation of parameters.
               // k values are forced to be positive
               // proportion values are forced to be between 0 and 1"""
               void reparameterisation(Vector_nD *x, int kSize, bool inverse=false);

               // Calculate parameters
               double operator()(Vector_nD x);
          };

          // Minimize the necessary scale factor K
          void findOptimalComparisonFunction();
     };

     template<class Archive>
     void VM2CosMarginalSampler::serialize(Archive & ar, const unsigned int version) {

         ar & k;
         ar & logNormConst;
         ar & root;

         ar & envelopeSampler;
         /*
         ar & envelopeDensity;
         ar & K;
         */

     }


}

#endif /* VONMISES2D_H_ */





