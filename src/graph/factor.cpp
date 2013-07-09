#include "factor.h"
#include "node.h"

#include "dai/util.h"

namespace pdbntk {

Factor createFactorIsing(Node *n) {
  return Factor(n);
}

/// Writes a factor to an output stream
/** \relates TFactor
*/
std::ostream& operator<< (std::ostream& os, const Factor& f) {
  os << "(" << f.nodes() << ")";
  return os;
}

Factor createFactorIsing(Node *n1, Node *n2) {
  DAI_ASSERT(n1 != n2);
  return Factor(NodeSet(n1, n2));
}


Factor createFactorExpGauss(const NodeSet &ns) {
  Factor fac(ns);
  return fac;
}

Factor createFactorPotts(Node *n1, Node *n2) {
  Factor fac(NodeSet(n1, n2));
  return fac;
}

Factor createFactorDelta(Node *v, size_t state ) {
  Factor fac(v);
  return fac;
}

Factor createFactorDelta( const NodeSet& vs, size_t state ) {
  Factor fac(vs);
  return fac;
}

Factor::Factor(Node *n) : ns_(n), states_(n->Size()) {}

Factor Factor::operator*(Real x) const {
 return Factor();
} 

Factor Factor::operator*=(Real x) const {
 return Factor();
}

Factor Factor::operator*=(const Factor &f) const {
  return Factor();
}

Factor Factor::operator*(const Factor &f) const {
  return *this;
}

Factor Factor::operator/(const Factor &f) const {
  return Factor();
}

Factor Factor::marginal(const NodeSet &nodes, bool normed) const {
  return Factor();
}

Factor Factor::maxMarginal(const NodeSet &nodes, bool normed) const {
  return Factor();
} 

Factor Factor::normalized() const {
  return Factor();
}

Real Factor::normalize() {
  return 0.0;
}

void Factor::InitPot() {}

} // end of namespace dai
