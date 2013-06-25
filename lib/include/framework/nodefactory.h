/*
 * nodefactory.h
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

// The NodeFactory creates the desired objects for you. You can also do this yourself, but it should be easier using the NodeFactory.

#ifndef NODEFACTORY_H_
#define NODEFACTORY_H_

#include "dbn.h"

namespace mocapy {

class NodeFactory {
public:
	NodeFactory();
	virtual ~NodeFactory();

	static DiscreteNode* new_discrete_node(uint node_size, const char* name = NULL, bool init_random=false, CPD new_cpd=CPD(), Node* discrete_node=NULL, bool fixed=false);
	static GaussianNode* new_gaussian_node(uint dimension, const char* name = NULL, bool init_random=false, bool shrinkage=false, MDArray<double> means=MDArray<double>(), MDArray<double> covs=MDArray<double>(), eCovType new_cov_type = FULL);
	static MultinomialNode* new_multinomial_node(uint dimension, const char* name = NULL, bool init_random=false, double pseudocount=0);
    static DirichletNode* new_dirichlet_node(uint dimension, const char* name = NULL, bool init_random=false, double pseudocount=0);
	static VonMises2dNode* new_vonmises2d_node(const char* name = NULL, MDArray<double> user_mus = MDArray<double>(), MDArray<double> user_kappas = MDArray<double>());
	static VonMisesNode* new_vonmises_node(const char* name = NULL, std::vector<double> user_mus = std::vector<double>(), std::vector<double> user_kappas = std::vector<double>());
	static KentNode* new_kent_node(const char* name = NULL, std::vector<double> kappas=std::vector<double>(), std::vector<double> betas=std::vector<double>(), std::vector<MDArray<double> > es=std::vector<MDArray<double> >(), double kappa_max=5000, double beta_max=1000);
	static PoissonNode* new_poisson_node(const char* name = NULL, std::vector<double> new_user_means=std::vector<double>());

	static BippoNode* new_bippo_node(const char* name = NULL, std::vector<double> new_lambdas=std::vector<double>(), std::vector<double> new_thetas=std::vector<double>(),std::vector<double> new_ns=std::vector<double>() );
};
}

#endif /* NODEFACTORY_H_ */
