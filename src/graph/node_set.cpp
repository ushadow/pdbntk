#include "node_set.h"
#include "../util.h"
#include "dai/util.h"

namespace pdbntk {
using std::vector;

size_t NodeSet::NStates() const {
  size_t s = 1;
  bforeach(Node* n, _elements) {
    s *= n->size();
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
}
