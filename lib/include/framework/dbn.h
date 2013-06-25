/*
 * DBN.h
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

#ifndef DBN_H_
#define DBN_H_

#include "../multinomial/MultinomialDensities.h"
#include "../multinomial/MultinomialESS.h"

#include <fstream>

#include <map>
#include <string>
#include <vector>
#include "../discrete/discreteess.h"
#include "../discrete/discretedensities.h"
#include "../gaussian/gaussiandensities.h"
#include "../gaussian/gaussianess.h"
#include "../poisson/poissondensities.h"
#include "../poisson/poissoness.h"
#include "childnode.h"
#include "../vonmises/vonmisesess.h"
#include "../vonmises/vonmisesdensities.h"
#include "../vonmises2d/vonmises2dess.h"
#include "../vonmises2d/vonmises2ddensities.h"
#include "../dirichlet/DirichletESS.h"
#include "../dirichlet/DirichletDensities.h"
#include "../bippo/bippoess.h"
#include "../bippo/bippodensities.h"
#include "../kent/kentess.h"
#include "../kent/kentdensities.h"
#include "../utils/utils.h"
#include "../utils/randomgen.h"


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

namespace mocapy {

typedef ChildNode<DiscreteESS, DiscreteDensities> DiscreteNode;
typedef ChildNode<GaussianESS, GaussianDensities> GaussianNode;
typedef ChildNode<VonMisesESS, VonMisesDensities> VonMisesNode;
typedef ChildNode<VonMises2dESS, VonMises2dDensities> VonMises2dNode;
typedef ChildNode<PoissonESS, PoissonDensities> PoissonNode;
typedef ChildNode<MultinomialESS, MultinomialDensities> MultinomialNode;
typedef ChildNode<DirichletESS, DirichletDensities> DirichletNode;
typedef ChildNode<KentESS, KentDensities> KentNode;
typedef ChildNode<BippoESS, BippoDensities> BippoNode;

// Dictionary that maps node names to (index, slice) tuples
#define IndexMap std::map<std::string, std::pair<int,int> >

class DBN {
public:
  DBN(uint seed=0);
	DBN(std::vector<Node*> & new_nodes_0, std::vector<Node*> & new_nodes_1,
	    std::string new_name = "DBN", uint seed=0);

	void load_dbn(const char* fname);
	// Add an edge between slices, from parent to child.
	void add_inter(NodeID & parent_i, NodeID & child_i);
	void add_inter(int parent_i, int child_i);
	void add_inter(const char* parent_i, const char* child_i);

	// Add an edge inside a slice, from parent to child.
	void add_intra(NodeID & parent_i, NodeID & child_i);
	void add_intra(int parent_i, int child_i);
	void add_intra(const char* parent_i, const char* child_i);

	// Initialize the DBN data structures based on the added nodes and edges.
	// After calling this method no edges can be added.
	void construct();

        uint get_total_output_size();
        uint get_nr_nodes() {return nr_nodes;}
	Node* get_node_by_name(std::string name);

	std::pair<std::vector<Node*> , std::vector<Node*> > get_nodes();
	std::pair<Sequence, double> sample_sequence(uint length);
	std::vector<Node*> getNodes0() { return nodes_0; }
	std::vector<Node*> getNodes1() { return nodes_1; }

	// Persistence
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const;

	template<class Archive>
	void load(Archive & ar, const unsigned int version);
	BOOST_SERIALIZATION_SPLIT_MEMBER()


	void save(const char* name);
	void load(const char* name);
	void set_slices(std::vector<Node*> new_nodes_0, std::vector<Node*> new_nodes_1);
	void setRandomGen(RandomGen* rg);
	void seed(uint s);
	
	friend class MCMC;
	friend class GibbsRandom;
	friend class FwbtRandom;
	friend class EMEngine;
	friend class AbstractInfEngine;
	friend class MoreAbstractInfEngine;
	friend std::ostream& operator<< (std::ostream& output, const DBN& a);

	RandomGen* randomGen;


protected:
	void make_index_map(std::vector<Node*> & nodes_0, std::vector<Node*> & nodes_1);
	std::pair<uint, uint> map_to_indices(NodeID & parent_i, NodeID & child_i);

	bool is_constructed;
	uint nr_nodes;
	uint total_output_size;
	std::string name;

	// Map the name of a node to its (index, slice).
	IndexMap index_map;

	// Nodes in slice 0
	std::vector<Node*> nodes_0;

	// Nodes in slices 1-T
	std::vector<Node*> nodes_1;

	std::vector<Node*> unique_nodes;
	std::vector<Node*> all_nodes;

	uint myseed;
};


template<class Archive>
void DBN::save(Archive & ar, const unsigned int version) const {
  
    ar.register_type(static_cast<DiscreteNode *>(NULL));
    ar.register_type(static_cast<GaussianNode *>(NULL));
    ar.register_type(static_cast<VonMisesNode *>(NULL));
    ar.register_type(static_cast<PoissonNode *>(NULL));
    ar.register_type(static_cast<MultinomialNode *>(NULL));
     ar.register_type(static_cast<KentNode *>(NULL));
    ar.register_type(static_cast<VonMises2dNode *>(NULL));


    if (version > 0) {
      ar.register_type(static_cast<PseudoCountPrior *>(NULL));
    }

    if (version > 1) {
      ar.register_type(static_cast<BippoNode *>(NULL));
      ar & myseed;
    }

    if (version > 2)
      ar.register_type(static_cast<DirichletNode *>(NULL));
 

      
    ar & is_constructed;
    ar & nr_nodes;
    ar & total_output_size;
    ar & name;
    ar & index_map;

    
    ar & nodes_0;
    ar & nodes_1;
   
    ar & unique_nodes;
    ar & all_nodes;
 
 }

template<class Archive>
void DBN::load(Archive & ar, const unsigned int version) {
    ar.register_type(static_cast<DiscreteNode *>(NULL));
    ar.register_type(static_cast<GaussianNode *>(NULL));
    ar.register_type(static_cast<VonMisesNode *>(NULL));
    ar.register_type(static_cast<PoissonNode *>(NULL));
    ar.register_type(static_cast<MultinomialNode *>(NULL));
    ar.register_type(static_cast<KentNode *>(NULL));
    ar.register_type(static_cast<VonMises2dNode *>(NULL));

    
    if (version > 0) {
      ar.register_type(static_cast<PseudoCountPrior *>(NULL));
    }

    if (version > 1) {
      ar.register_type(static_cast<BippoNode *>(NULL));
      ar & myseed;
    }

    if (version > 2)
      ar.register_type(static_cast<DirichletNode *>(NULL));
    
    ar & is_constructed;
    ar & nr_nodes;
    ar & total_output_size;
    ar & name;
    ar & index_map;
    ar & nodes_0;
    ar & nodes_1;
    ar & unique_nodes;
    ar & all_nodes;
    
    randomGen = new RandomGen();
    setRandomGen(randomGen);

    if (version > 1) {
      seed(myseed);
    }
 }



 std::ostream& operator <<(std::ostream &os,const DBN &obj);

}


BOOST_CLASS_VERSION(mocapy::DBN, 3)

#endif /* DBN_H_ */
