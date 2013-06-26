#include "jtree_2tbn_inf_engine.h"

namespace pdbntk {

JTree2TBNInfEngine::JTree2TBNInfEngine(mocapy::DBN *dbn,
    const std::vector<mocapy::NodeID> &interface) {
}

double JTree2TBNInfEngine::EnterEvidence(const mocapy::Sequence &evidence) {
  return 0.0;
}

std::vector<mocapy::ESSBase*> JTree2TBNInfEngine::GetResetESS() const {
  return std::vector<mocapy::ESSBase*>();
}


}
