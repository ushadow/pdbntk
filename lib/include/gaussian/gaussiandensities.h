/*
 * gaussiandensities.h
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

#ifndef GAUSSIANDENSITIES_H_
#define GAUSSIANDENSITIES_H_


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>


#include "../framework/densitiesbase.h"
#include "../framework/mocapyexceptions.h"
#include "multigauss.h"

namespace mocapy {


enum eCovType {
	FULL, SPHE, DIAG, WRONG_COV_TYPE
};

class GaussianDensities;
std::ostream& operator<<(std::ostream&, const GaussianDensities&);

class GaussianDensities: public DensitiesBase {
public:
  GaussianDensities() { initialize(); useShrinkage=false; }
	GaussianDensities(uint new_dim, bool new_init_random=false, eCovType new_cov_type = FULL, bool new_cov_tied = false,
			MDArray<double> new_user_means = MDArray<double> (),
			MDArray<double> new_user_covs = MDArray<double> ());

	virtual ~GaussianDensities() {};

	// Called in node.construct, and initializes the density arrays
	void construct(std::vector<uint> & parent_sizes);

	// Parameter estimation based on the ESS
	void estimate(std::vector< MDArray<double> > & ess);

	// Return a sample, based on indicated parent values
	std::vector<double> sample(std::vector<double> & ptv);

	// Return likelihood, that is: P(child|parents)
	double get_lik(std::vector<double> & ptv, bool log=false);

	// Return the distribtion's parameters
	std::vector< MDArray<double> > get_parameters();

	bool useShrinkage;

	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    friend std::ostream& operator<< (std::ostream& output, const GaussianDensities& a);

private:
	// Covariance tied?
	bool cov_tied;
	uint par_count;
	uint parent_index;
	bool init_random;
	// Covariance SPHE, DIAG or FULL
	eCovType cov_type;
	MDArray<double> user_means;
	MDArray<double> user_covs;
	MDArray<double> means;
	MDArray<double> covs;

	std::vector<uint> mean_shape;
	std::vector<uint> cov_shape;

	std::vector<MultiGauss> mgauss_list;

	MDArray<double> calc_means(std::vector<MDArray<double> > & ess);
	MDArray<double> calc_cov_full(std::vector<MDArray<double> > & ess);
	MDArray<double> calc_cov_full_tied(std::vector<MDArray<double> > & ess);
	MDArray<double> multivariate_normal(MDArray<double> & m, MDArray<double> & c);
	MDArray<double> multivariate_normal_DUMMY(MDArray<double> & m, MDArray<double> & c);
	MDArray<double> make_rnd_covs(uint node_size, uint dim, std::vector<uint> & cov_shape, bool tied);
	MDArray<double> make_uniform_covs(uint node_size, uint dim, std::vector<uint> & cov_shape);
	std::vector<MultiGauss> make_mgauss_list(MDArray<double> & means, MDArray<double> & covs);
	void initialize();
};


template<class Archive>
void GaussianDensities::serialize(Archive & ar, const unsigned int version) {
  
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    
    ar & cov_tied;
	ar & par_count;
	ar & parent_index;
	ar & cov_type;

	
	ar & user_means;
	ar & user_covs;
	ar & means;
	ar & covs;
	ar & mean_shape;
	ar & cov_shape;

	
	ar & mgauss_list;
    
	
}

}

#endif /* GAUSSIANDENSITIES_H_ */
