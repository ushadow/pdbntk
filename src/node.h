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

#ifndef PDBNTK_NODE_H_
#define PDBNTK_NODE_H_

#include "mocapy.h"
#include "cond_prob_dist.h"

#include <string>
#include <vector>

namespace pdbntk {

enum eSliceType {
	TIED, START, END
};

class Node;

std::ostream& operator<<(std::ostream&, const Node&);

/// A node in the Bayesian network.
class Node {
public:
	Node(uint node_index, CondProbDist *cpd); 
	virtual ~Node();
	void set_node_index(uint ni);
	void set_data_index(uint di);

	void add_intra_child(Node* n);
	void add_inter_child(Node* n);
	void add_inter_parent(uint data_index, uint node_size);
	void add_intra_parent(uint data_index, uint node_size);
	void fix(bool flag);
	std::string get_name() {return name;}
	void set_name(const char* new_name) {name = new_name;}
  virtual void set_parentmap(mocapy::ParentMap * pm);

  bool operator< (const Node &n) const;
  bool operator> (const Node &n) const;
  bool operator== (const Node &n) const;
  bool operator<= (const Node &n) const;
  bool operator>= (const Node &n) const;
  bool operator!= (const Node &n) const;

  // Persistence
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

  // Index of node data in slice
  int data_index;
  bool fixed;

protected:
  // Index of node in node list
  uint node_index_;

  mocapy::ParentMap parentmap;

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

  /// Conditional probability distribution associated with this node.
  CondProbDist *cpd_;
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
  ar & node_index_;
  ar & seq_len;
  ar & slice;
  ar & name;
}

}
#endif // PDBNTK_NODE_H_ 
