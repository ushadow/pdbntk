#ifndef PDBNTK_NODE_H_
#define PDBNTK_NODE_H_

#include "../cpd/cond_prob_dist.h"
#include "factor.h"

#include <string>
#include <vector>
#include <iostream>
#include <memory>

namespace pdbntk {

class CondProbDist;

enum eSliceType {
	TIED, START, END
};

/// A node in the Bayesian network.
class Node {
public:
  Node() : factor_(new Factor(this)) {};
  /// \param cpd cannot be NULL, and this node does not take the ownership of
  //  the pointer.
	Node(uint node_index, uint di, CondProbDist *cpd, bool observed = false); 
	virtual ~Node();

/// Accessors
//@{
  /// If the node is discrete, its size is the number of possible values it can
  /// take on; if the node is continuous, it can be a vector and its size is the 
  /// length of this vector. If the node is observed, its size is 1.
  uint Size() const; 
  uint node_index() const { return node_index_; }
  bool observed() const { return observed_; }
	std::string get_name() const { return name_; }
//@}

/// Mutators
//@{
	void set_node_index(uint ni);
	void set_data_index(uint di);
	void add_parent(uint node_index, uint node_size);
	void set_name(const char* new_name) {name_ = new_name;}
  /// Updates the CPD factor this node belongs to.
  const Factor& UpdateFactor(const std::vector<Real> &ev);

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


protected:
  // Index of node data in slice
  int data_index_;
  /// Index of node in node list
  uint node_index_;

  bool observed_;

  /// Conditional probability distribution associated with this node.
  CondProbDist *cpd_;

  std::vector<uint> parents_;

  std::vector<uint> parents_sizes_;

  // 'slice' indicates to which slice the node belongs
  eSliceType slice;

  std::string name_;

  std::unique_ptr<Factor> factor_;
};

template<class Archive>
void Node::serialize(Archive & ar, const unsigned int version) {
  ar & data_index_;
  ar & parents_;
  ar & parents_sizes_;
  ar & node_index_;
  ar & slice;
  ar & name_;
}

inline std::ostream& operator<<(std::ostream& os, const Node *node) {
  os << node->node_index();
  return os;
}

}
#endif // PDBNTK_NODE_H_ 
