#include "node.h"

#include "mocapy.h"

using namespace std;

namespace pdbntk {

Node::Node(uint ni, CondProbDist *cpd) : 
    index_(ni), cpd_(cpd), fixed(false), is_constructed(false) {}

Node::~Node() {} 

void Node::set_data_index(uint di) {
	data_index = di;
}

void Node::add_intra_child(Node* n) {
	// Add an intra-slice child.
	// child node (same slice) of current node
	children_1.push_back(n);
}

void Node::add_inter_child(Node* n) {
	// Add a child in the previous slice.
	// child node (previous slice) of current node
	children_2.push_back(n);
}

void Node::add_inter_parent(uint data_index, uint node_size) {
	// Add a parent in the previous slice.
	// node_size: output size of parent
	parents_0.push_back(data_index);
	parents_0_sizes.push_back(node_size);
}

void Node::add_intra_parent(uint data_index, uint node_size) {
	// Add an intra-slice parent.
	// node_size: output size of parent
	parents_1.push_back(data_index);
	parents_1_sizes.push_back(node_size);
}

void Node::fix(bool flag) {
	// Fix the node (flag=1), or unfix the node (flag=0).
	fixed = flag;
}

void Node::set_parentmap(mocapy::ParentMap * pm) {
	parentmap = *pm;
	seq_len = parentmap.lng;
}

bool Node::operator< (const Node &n) const {
  return index_ < n.index_;
}

bool Node::operator> (const Node &n) const { 
  return index_ > n.index_; 
}

bool Node::operator>= (const Node &n) const { 
  return index_ >= n.index_; 
}

bool Node::operator<= (const Node &n) const { 
  return index_ <= n.index_; 
}

bool Node::operator!= (const Node &n) const { 
  return index_ > n.index_; 
}
}
