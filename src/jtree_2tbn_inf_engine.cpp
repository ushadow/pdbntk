#include "jtree_2tbn_inf_engine.h"

#include "dai/properties.h"

#include <string>

namespace pdbntk {

JTree2TBNInfEngine::JTree2TBNInfEngine(const FactorGraph &fg1,
    const FactorGraph &fg15, const NodeSet &interface1,
    const NodeSet &interface2) {
  dai::PropertySet infprops;
  infprops.set("updates", std::string("HUGIN"));
  infprops.set("root", interface1);
  jtree_engine1_.reset(new JTree(fg1, infprops)); 

  // Creates a jtree engine for 1.5 slice.
  infprops.set("root", interface2);
  jtree_engine_.reset(new JTree(fg15, infprops)); 
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

void JTree2TBNInfEngine::Fwd1(const Evidence::Observation &o, int t) {
  for (Evidence::Observation::const_iterator i = o.begin(); i != o.end(); i++)
    jtree_engine1_->clamp(jtree_engine1_->fg().findNode(i->first), i->second);
  jtree_engine1_->init();
  jtree_engine1_->run();
}


}
