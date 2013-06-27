/*
 * Node.cpp
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

#include "mocapyexceptions.h"
#include "node.h"

using namespace std;

namespace mocapy {

Node::~Node() {

}

void Node::set_node_index(uint ni) {
	node_index = ni;
}

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

void Node::set_parentmap(ParentMap * pm) {
	parentmap = *pm;
	weight = parentmap.get_weight();
	seq_len = parentmap.lng;
}

ostream& operator<<(ostream& output, const Node& a)  {
	output << a.name << endl;
	a.printNode(output);
	return output;
}

}
