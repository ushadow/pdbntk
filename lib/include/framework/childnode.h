/*
 * childnode.h
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

#ifndef CHILDNODE_H_
#define CHILDNODE_H_

#include <boost/serialization/base_object.hpp>
#include <fstream>


#include "../discrete/discreteess.h"
#include "../discrete/discretedensities.h"
#include "../gaussian/gaussiandensities.h"
#include "../gaussian/gaussianess.h"
#include "node.h"
#include "../vonmises/vonmisesess.h"
#include "../vonmises/vonmisesdensities.h"

namespace mocapy {

// Forwards
template<class ESS, class Densities>
class ChildNode;

template<class ESS, class Densities> std::ostream& operator<<(std::ostream&, const ChildNode<ESS, Densities>&);

template<class ESS, class Densities>
class ChildNode: public Node {
public:
	ChildNode() {};
	virtual ~ChildNode() {};

	// Return a sample of the i'th slice
	void sample(uint i);

	// Add a new sample to the ESS
	void update_ess(std::vector<double> & ptv);

	// Flag end of sequence sampling
	void save_ess();

	// Blanket sample discrete nodes only
	void blanket_sample(uint i, MDArray<eMISMASK> & mismask);

	// Compute LL of i'th slice
	double get_slice_log_likelihood(uint i);

	// Compute LL of i'th slice and also return ptv
	std::pair<double, std::vector<double> > get_slice_log_likelihood_ptv(uint i);

	std::vector< MDArray<double> > get_ess();

	void do_M_step(std::vector<MDArray<double> > & ess);

	// Initialize ESS and densities objects with the DBN specific parameters
	void construct(eSliceType new_slice);

	// Construct and return parentmap with node-specific parameters
	ParentMap get_parentmap(Sequence * seq, double weight);

	// Return the parent (base) class
	Node* base() { return dynamic_cast<Node*>(this); }

	void set_densities(Densities d) { densities = d; }
	Densities* get_densities() { return &densities; }
	std::vector<MDArray<double> > get_parameters();

	virtual void printNode(std::ostream& output) const;

	void setRandomGen(RandomGen* rg);

	// Dimension
	uint get_output_size();
	uint get_node_size();
	eNodeType get_node_type();


	ESS& get_ess_ref() {return ess;}

	// Persistence
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const;

	template<class Archive>
	void load(Archive & ar, const unsigned int version);
	BOOST_SERIALIZATION_SPLIT_MEMBER()
	  /*
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
	  */

	//Operators
	friend std::ostream& operator<< <ESS, Densities> (std::ostream& output, const ChildNode<ESS, Densities>& a);
	Densities densities;

private:
	ESS ess;
};


template<class ESS, class Densities>
void ChildNode<ESS, Densities>::setRandomGen(RandomGen* rg) {
  randomGen = rg;
  densities.setRandomGen(rg);
}

template<class ESS, class Densities>
inline double ChildNode<ESS, Densities>::get_slice_log_likelihood(uint i) {
	std::vector<double> ptv;
	parentmap.get(i, ptv);
	return densities.get_lik(ptv, true);
}

template<class ESS, class Densities>
void ChildNode<ESS, Densities>::sample(uint i) {
  assert(randomGen != NULL);
	std::vector<double> pv;
	parentmap.get(i, pv);
	pv.pop_back();

	assert(densities.randomGen != NULL);
	std::vector<double> choice = densities.sample(pv);

	parentmap.set(i, choice);
}

template<class ESS, class Densities>
void ChildNode<ESS, Densities>::blanket_sample(uint i, MDArray<eMISMASK> & mismask) {
	// Blanket sample for discrete nodes
	if (densities.get_node_type() != DISCRETE)
		return sample(i);

	// Get parent+this values
	std::vector<double> ptv;
	parentmap.get(i, ptv);

	uint ptv_sz = ptv.size();
	assert(ptv_sz>0);

	// Probabilities
	CPD p;
	std::vector<uint> p_shape;
	p_shape.push_back(get_node_size());
	p.set_shape(p_shape);

	// Loop over all proposed new states
	for (uint state = 0; state < get_node_size(); state++) {
		// Set value of node in slice_1 to state
		ptv[ptv_sz - 1] = (double)state;
		parentmap.set(i, state);
		// Note that we are working in LOG space!
		p.set(state, densities.get_lik(ptv,true) );
		double log_lik_product(0);
		// Children in the current slice
		for (uint j = 0; j < children_1.size(); j++) {
			Node* child = children_1[j];

			uint index = child->node_index;
			eMISMASK mm;
			if (mismask.get_shape()[0] == 1) {
				// Mismask is only specified for a slice
				mm = mismask.get(0,index);
			} else {
				mm = mismask.get(i, index);
			}

			if (mm != MOCAPY_MISSING) {
				double lik = child->get_slice_log_likelihood(i);
				log_lik_product += lik;
			}
		}
		// Children in next slice
		if (i + 1 < seq_len) {
			// Skip this when we reach end of sequence
			// (because there are no children in the next slice)
			for (uint j = 0; j < children_2.size(); j++) {
				Node* child = children_2[j];


				uint index = child->node_index;
				eMISMASK mm;
				if (mismask.get_shape()[0] == 1) {
					// Mismask is only specified for a slice
					mm = mismask.get(0,index);
				} else {
					mm = mismask.get(i, index);
				}

				if (mm != MOCAPY_MISSING) {
					double lik = child->get_slice_log_likelihood(i + 1);
					log_lik_product += lik;
				}
			}
		}

		double a = p[state];
		p.set(state, a + log_lik_product);
	}
	// Get rid of potentially very large values in p
	double m = -p.get_max();
	p.add_inplace(m);

	// Get out of LOG space
	p.exp_all();

	// Normalize probabilities
	p.normalize();

	// Make cumulative cpd
	p.cumsum();

	// Choose a node value
	double r = randomGen->get_rand();

	assert(r<=1 && r>=0);
	uint choice = p.bisect(r);
	parentmap.set(i, choice);
}


template<class ESS, class Densities>
void ChildNode<ESS, Densities>::update_ess(std::vector<double> & ptv) {
	ess.add_ptv(ptv);
}

template<class ESS, class Densities>
void ChildNode<ESS, Densities>::save_ess() {
	ess.save_ess(weight);
}

template<class ESS, class Densities>
std::pair<double, std::vector<double> > ChildNode<ESS, Densities>::get_slice_log_likelihood_ptv(uint i) {
	std::vector<double> ptv;
	parentmap.get(i, ptv);
	return make_pair(densities.get_lik(ptv, true), ptv);
}

template<class ESS, class Densities>
std::vector<MDArray<double> > ChildNode<ESS, Densities>::get_ess() {
	return ess.get_array();
}

template<class ESS, class Densities>
void ChildNode<ESS, Densities>::do_M_step(std::vector<MDArray<double> > & new_ess) {
    // Update the parameters
     densities.estimate(new_ess);
     // Get ready for new cycle
     ess.clear();
}

// Called from DBN
template<class ESS, class Densities>
void ChildNode<ESS, Densities>::construct(eSliceType new_slice) {
  assert(randomGen != NULL);
	std::vector<uint> parent_sizes = vec_conc(parents_0_sizes, parents_1_sizes);
	densities.construct(parent_sizes);
	ess.construct(parent_sizes, get_output_size(), get_node_size());
	slice = new_slice;
	is_constructed = true;
}

template<class ESS, class Densities>
uint ChildNode<ESS, Densities>::get_output_size() {
	return densities.get_output_size();
}

template<class ESS, class Densities>
uint ChildNode<ESS, Densities>::get_node_size() {
	return densities.get_node_size();
}

template<class ESS, class Densities>
eNodeType ChildNode<ESS, Densities>::get_node_type() {
	return densities.get_node_type();
}


template<class ESS, class Densities>
ParentMap ChildNode<ESS, Densities>::get_parentmap(Sequence * seq, double weight) {
	return ParentMap(seq, parents_0, parents_1, data_index, weight, densities.get_output_size());
}


/*
template<class ESS, class Densities> template<class Archive>
void ChildNode<ESS, Densities>::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<Node>(*this);
    ar & densities;
    
    std::vector<uint> parent_sizes = vec_conc(parents_0_sizes, parents_1_sizes);
    ess.construct(parent_sizes, get_output_size(), get_node_size());
}
*/

 template<class ESS, class Densities> template<class Archive>
  void ChildNode<ESS, Densities>::load(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<Node>(*this);
    ar & densities;
    
    std::vector<uint> parent_sizes = vec_conc(parents_0_sizes, parents_1_sizes);
    ess.construct(parent_sizes, get_output_size(), get_node_size());

 }


 template<class ESS, class Densities> template<class Archive>
  void ChildNode<ESS, Densities>::save(Archive & ar, const unsigned int version) const {
   ar & boost::serialization::base_object<Node>(*this);
    ar & densities;
 }

template<class ESS, class Densities>
std::vector<MDArray<double> > ChildNode<ESS, Densities>::get_parameters() {
        return densities.get_parameters();
}

template<class ESS, class Densities>
  void ChildNode<ESS, Densities>::printNode(std::ostream& output) const {
  output << densities << std::endl;
}

}


#endif /* CHILDNODE_H_ */
