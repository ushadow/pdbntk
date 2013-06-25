/*
 * Node.h
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

#ifndef NODE_H_
#define NODE_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <string>
#include <vector>
#include "essbase.h"
#include "../utils/mdarray.h"
#include "parentmap.h"

namespace mocapy {

enum eNodeType {
     DISCRETE, GAUSSIAN, DIRICHLET, KENT, VONMISES, VONMISES2D, POISSON, MULTINOMIAL, BIPPO
};
enum eSliceType {
	TIED, START, END
};

std::ostream& operator<<(std::ostream&, const Node&);

class Node {
public:
	Node() { fixed=false; is_constructed=false; randomGen=NULL; }
	virtual ~Node();
	void set_node_index(uint ni);
	void set_data_index(uint di);

	void add_intra_child(Node* n);
	void add_inter_child(Node* n);
	void add_inter_parent(uint data_index, uint node_size);
	void add_intra_parent(uint data_index, uint node_size);
	virtual void construct(eSliceType new_slice) = 0;
	void fix(bool flag);
	virtual uint get_output_size()=0;
	virtual uint get_node_size()=0;
	virtual eNodeType get_node_type()=0;
	std::string get_name() {return name;}
	void set_name(const char* new_name) {name = new_name;}
	virtual void sample(uint i)=0;
	virtual void update_ess(std::vector<double> & ess)=0;
	virtual void save_ess()=0;
	virtual void blanket_sample(uint i, MDArray<eMISMASK> & mismask)=0;
	virtual std::pair<double, std::vector<double> > get_slice_log_likelihood_ptv(uint i)=0;
	virtual double get_slice_log_likelihood(uint i)=0;
	virtual ParentMap get_parentmap(Sequence * seq, double weight)=0;
	virtual void set_parentmap(ParentMap * pm);
	virtual std::vector<MDArray<double> > get_ess()=0;
	virtual void do_M_step(std::vector<MDArray<double> > & ess)=0;
    virtual std::vector<MDArray<double> > get_parameters() = 0;
    virtual void printNode(std::ostream& output) const = 0;

    virtual void setRandomGen(RandomGen* rg) = 0; 


	// Persistence
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	// Index of node data in slice
	int data_index;
	bool fixed;

	// Index of node in node list
	uint node_index;

	RandomGen* randomGen;


	friend class GibbsRandom;
	friend class FwbtRandom;
	friend class EMEngine;
	friend class ForwardBacktracker;
	friend class GeneralInfEngineMM;
	friend std::ostream& operator<< (std::ostream& output, const Node& a);

protected:
	ParentMap parentmap;

	// Children in same slice
	std::vector<Node*> children_1;

	// Children in next slice
	std::vector<Node*> children_2;

	// Parents in previous slice
	std::vector<uint> parents_0;

	// Parents in same slice
	std::vector<uint> parents_1;

	std::vector<uint> parents_0_sizes;
	std::vector<uint> parents_1_sizes;

	bool is_constructed;



	uint seq_len;

	// 'slice' indicates to which slice the node belongs
	eSliceType slice;

	std::string name;

	double weight;
};


template<class Archive>
void Node::serialize(Archive & ar, const unsigned int version) {
	ar & data_index;
	ar & fixed;

	if (version == 0) {
		ar & parentmap;
	}

	ar & children_1;
	ar & children_2;
	ar & parents_0;
	ar & parents_1;
	ar & parents_0_sizes;
	ar & parents_1_sizes;
	ar & is_constructed;
	ar & node_index;
	ar & seq_len;
	ar & slice;
	ar & name;
	ar & weight;

}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Node)

}

BOOST_CLASS_VERSION(mocapy::Node, 1)

#endif /* NODE_H_ */
