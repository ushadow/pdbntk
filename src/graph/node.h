#ifndef PDBNTK_NODE_H_
#define PDBNTK_NODE_H_

#include "../cpd/cond_prob_dist.h"
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
  Node() {};
  /// \param cpd cannot be NULL, and this node does not take the ownership of
  //  the pointer.
	Node(uint node_index, CondProbDist *cpd, bool observed = false); 
	virtual ~Node();
	void set_node_index(uint ni);
	void set_data_index(uint di);

/// Accessors
//@{
  /// If the node is discrete, its size is the number of possible values it can
  /// take on; if the node is continuous, it can be a vector and its size is the 
  /// length of this vector. If the node is observed, its size is 1.
  uint size() const; 
  uint index() const { return index_; }
//@}

	void add_intra_child(Node* n);
	void add_inter_child(Node* n);
	void add_inter_parent(uint data_index, uint node_size);
	void add_intra_parent(uint data_index, uint node_size);
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

protected:
  /// Index of node in node list
  uint index_;

  bool observed_;

  /// Conditional probability distribution associated with this node.
  CondProbDist *cpd_;

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
};

template<class Archive>
void Node::serialize(Archive & ar, const unsigned int version) {
  ar & data_index;

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
