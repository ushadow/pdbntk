#include "node.h"

#include "mocapy.h"

using namespace std;

namespace pdbntk {

Node::Node(uint ni, uint di, CondProbDist *cpd, bool observed) : 
    node_index_(ni), data_index_(di), cpd_(cpd), observed_(observed) {}

Node::~Node() {} 

uint Node::Size() const {
  if (observed_)
    return 1;
  else
    return cpd_->NodeSize();
}

void Node::set_data_index(uint di) {
	data_index_ = di;
}

void Node::add_parent(uint data_index, uint node_size) {
	// Add an intra-slice parent.
	// node_size: output size of parent
	parents_.push_back(data_index);
	parents_sizes_.push_back(node_size);
}

bool Node::operator< (const Node &n) const {
  return node_index_ < n.node_index_;
}

bool Node::operator> (const Node &n) const { 
  return node_index_ > n.node_index_; 
}

bool Node::operator>= (const Node &n) const { 
  return node_index_ >= n.node_index_; 
}

bool Node::operator<= (const Node &n) const { 
  return node_index_ <= n.node_index_; 
}

bool Node::operator!= (const Node &n) const { 
  return node_index_ > n.node_index_; 
}

const Factor& Node::UpdateFactor(const std::vector<Real> &ev) {
  
  return *factor_; 
}
}
