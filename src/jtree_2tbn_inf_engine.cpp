#include "jtree_2tbn_inf_engine.h"

#include "dai/properties.h"

#include <string>

namespace pdbntk {

JTree2TBNInfEngine::JTree2TBNInfEngine(const FactorGraph &fg15,
    const FactorGraph &fg1) {
  dai::PropertySet infprops;
  infprops.set("updates", std::string("HUGIN"));
  // Creates a jtree engine for 1.5 slice.
  jtree_engine_.reset(new JTree(fg15, infprops)); 
  jtree_engine1_.reset(new JTree(fg1, infprops)); 
}

double JTree2TBNInfEngine::EnterEvidence(const mocapy::Sequence &evidence) {
  return 0.0;
}

std::vector<mocapy::ESSBase*> JTree2TBNInfEngine::GetResetESS() const {
  return std::vector<mocapy::ESSBase*>();
}

void JTree2TBNInfEngine::Fwd(const Evidence::Observation &o, int t) {
  /// Computes prior

  /// Clamps the observed nodes.
  for (Evidence::Observation::const_iterator i = o.begin(); i != o.end(); i++)
    jtree_engine_->clamp(jtree_engine_->fg().findNode(i->first), i->second);

  /// Collect evidence to root
}


}
