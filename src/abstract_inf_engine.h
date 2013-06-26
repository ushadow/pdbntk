#ifndef ABSTRACT_INF_ENGINE_H_
#define ABSTRACT_INF_ENGINE_H_

#include "mocapy.h"

#include <vector>

namespace pdbntk {

class AbstractInfEngine {
public:
  AbstractInfEngine() {};
  virtual ~AbstractInfEngine() {};
  virtual double EnterEvidence(const mocapy::Sequence &evidence) = 0;
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const = 0;
};
}
#endif // ABSTRACT_INF_ENGINE_H_
