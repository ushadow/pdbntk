#ifndef PDBNTK_NODE_H_
#define PDBNTK_NODE_H_

#include "cpd/cond_prob_dist.h"
#include "framework/parentmap.h"

#include <string>
#include <vector>
#include <iostream>

namespace pdbntk {

class CondProbDist;

enum eSliceType {
	TIED, START, END
};

/// A node in the Bayesian network.
class Node {
public:
	Node(uint node_index, CondProbDist *cpd); 
	virtual ~Node();
	void set_node_index(uint ni);
	void set_data_index(uint di);

/// Accessors
//@{
  uint size() const { return cpd_->node_size(); }
  uint index() const { return index_; }
//@}

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
  uint index_;

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
  ar & index_;
  ar & seq_len;
  ar & slice;
  ar & name;
}

inline std::ostream& operator<<(std::ostream& os, const Node *node) {
  os << node->index();
  return os;
}
}
#endif // PDBNTK_NODE_H_ 
