#ifndef JTREE_2TBN_INF_ENGINE_
#define JTREE_2TBN_INF_ENGINE_

#include "abstract_inf_engine.h"
#include "jtree_inf_engine.h"

#include "mocapy.h"

#include <boost/scoped_ptr.hpp>
#include <vector>

namespace pdbntk {
class JTree2TBNInfEngine : public AbstractInfEngine {
public:
  JTree2TBNInfEngine(mocapy::DBN *dbn,
      const std::vector<mocapy::NodeID> &interface);
  virtual double EnterEvidence(const mocapy::Sequence &evidence);
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const;

private:
  boost::scoped_ptr<JTreeInfEngine> jtree_engine_, jtree_engine1_;
};
}
#endif // JTREE_2TBN_INF_ENGINE_
