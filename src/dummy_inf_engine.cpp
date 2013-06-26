#include "dummy_inf_engine.h"

#include <math.h>

namespace pdbntk {

double DummyInfEngine::EnterEvidence(const mocapy::Sequence &evidence) {
  for (int i = 0; i < 10000; i++) {
    for (int j = 0; j < 10000; j++) {
      double a = sqrt(sqrt(i) + sin(cos(j)));  
    }
  }
  return 1;
}

std::vector<mocapy::ESSBase*> DummyInfEngine::GetResetESS() const {
  std::vector<mocapy::ESSBase*> ess;
  return ess;
}
}
