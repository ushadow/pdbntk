#include "node_set.h"
#include "node.h"
#include "../util.h"
#include "dai/util.h"

namespace pdbntk {
using std::vector;

bool NodeComparator::operator()(const Node *n1, const Node *n2) const {
  return *n1 < *n2;
}

size_t NodeSet::NStates() const {
  size_t s = 1;
  bforeach(Node* n, _elements) {
    s *= n->Size();
  }
  return s;
}

std::ostream& operator<<(std::ostream& os, const vector<NodeSet>& v) {
  os << "(";
  for (vector<NodeSet>::const_iterator it = v.begin(); it != v.end(); it++) {
    os << (it != v.begin() ? ", " : "") << *it;
  }
  os << ")";
  return os;
}

NodeSet::NodeSet(Node *t1, Node *t2) {
  if(*t1 < *t2) {
    _elements.push_back(t1);
    _elements.push_back(t2);
  } else if (*t2 < *t1) {
    _elements.push_back(t2);
    _elements.push_back(t1);
  } else
    _elements.push_back(t1);
}


}
