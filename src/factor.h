/// \file
/// \brief Defines Factor class which represent factors in probability distributions.

#ifndef PDBNTK_FACTOR_H_ 
#define PDBNTK_FACTOR_H_ 

#include "node.h"
#include "utils.h"

#include <iostream>
#include <functional>
#include <cmath>
#include <memory>

namespace pdbntk {

/// Represents a potential factor.
/** Mathematically, a \e factor is a function mapping joint states of some
 *  variables to the nonnegative real numbers.
 *  More formally, denoting a variable with label \f$l\f$ by
 *  \f$x_l\f$,
 *  a factor depending on the variables \f$\{x_l\}_{l\in L}\f$ is
 *  a function \f$f_L : \prod_{l\in L} X_l \to [0,\infty)\f$.
 *
 *  In libDAI, a factor is represented by a Factor object, which has two
 *  components:
 *  \arg a VarSet, corresponding with the set of variables \f$\{x_l\}_{l\in L}\f$
 *  that the factor depends on;
 *  \arg a TProb, a vector containing the value of the factor for each possible
 *  joint state of the variables.
 *
 *  The factor values are stored in the entries of the TProb in a particular
 *  ordering, which is defined by the one-to-one correspondence of a joint state
 *  in \f$\prod_{l\in L} X_l\f$ with a linear index in
 *  \f$\{0,1,\dots,\prod_{l\in L} S_l-1\}\f$ according to the mapping \f$\sigma\f$
 *  induced by dai::calcLinearState().
 *
 *  \tparam T Should be a scalar that is castable from and to double and should support elementary arithmetic operations.
 *  \todo Define a better fileformat for .fg files (maybe using XML)?
 *  \todo Add support for sparse factors.
 */
class Factor {
  private:
    /// Stores the nodes on which the factor depends.
    NodeSet ns_;
    /// Number of states of this potential factor.
    size_t states_;

  public:
    Factor() : ns_() {}
    /// Constructs factor depending on the node \a n. 
    Factor(Node* n) : ns_(n) {}

    /// Constructs factor depending on variables in \a vars with uniform distribution
    Factor(const NodeSet& nodes) : ns_(nodes) {}

    /// Returns constant reference to variable set (i.e., the variables on which the factor depends)
    const NodeSet& nodes() const { return ns_; }

    /// Returns reference to variable set (i.e., the variables on which the factor depends)
    NodeSet& nodes() { return ns_; }

    size_t states() { return states_; }

    /// Comparison
    bool operator==( const Factor& y ) const {
      return (ns_ == y.ns_);
    }

    /// Returns marginal on \a nodes, obtained by summing out all variables except those in \a nodes, and normalizing the result if \a normed == \c true
    Factor marginal(const NodeSet &nodes, bool normed=true) const;

    /// Returns max-marginal on \a vars, obtained by maximizing all variables except those in \a vars, and normalizing the result if \a normed == \c true
    Factor maxMarginal(const NodeSet &nodes, bool normed=true) const;
    //@}

    Factor normalized();
    Real normalize();
    Factor operator* (Real x) const;
    Factor operator*= (Real x) const;
    Factor operator* (const Factor& f) const;
    Factor operator/ (const Factor& f) const;
    Factor operator*= (const Factor& f) const;
    Factor operator/= (const Factor& f) const;
};

std::ostream& operator<< (std::ostream& os, const Factor& f);

} 

#endif // PDBNTK_FACTOR_H_
