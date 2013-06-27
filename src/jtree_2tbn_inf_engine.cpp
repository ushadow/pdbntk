#include "jtree_2tbn_inf_engine.h"

#include "dai/properties.h"

#include <string>

namespace pdbntk {

JTree2TBNInfEngine::JTree2TBNInfEngine(const dai::FactorGraph &fg15,
    const dai::FactorGraph &fg1) {
  dai::PropertySet infprops;
  infprops.set("updates", std::string("HUGIN"));
  jtree_engine_.reset(new dai::JTree(fg15, infprops)); 
  jtree_engine1_.reset(new dai::JTree(fg1, infprops)); 
}

double JTree2TBNInfEngine::EnterEvidence(const mocapy::Sequence &evidence) {
  return 0.0;
}

std::vector<mocapy::ESSBase*> JTree2TBNInfEngine::GetResetESS() const {
  return std::vector<mocapy::ESSBase*>();
}


}
