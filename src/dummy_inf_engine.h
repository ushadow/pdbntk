#ifndef DUMMY_INF_ENGINE_H_
#define DUMMY_INF_ENGINE_H_

#include "abstract_inf_engine.h"
#include "mocapy.h"

#include <vector>

namespace pdbntk {

class DummyInfEngine : public AbstractInfEngine {
public:
  virtual double EnterEvidence(const mocapy::Sequence &evidence);
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const;
};
}
#endif // DUMMY_INF_ENGINE_H_
