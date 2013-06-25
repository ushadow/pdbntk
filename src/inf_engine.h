#ifndef INF_ENGINE_H_
#define INF_ENGINE_H_

#include "mocapy.h"

#include <vector>

namespace pdbntk {

class InfEngine {
public:
  double EnterEvidence(const mocapy::Sequence &evidence);
  std::vector<mocapy::ESSBase*> GetResetESS();
};
}
#endif // INF_ENGINE_H_
