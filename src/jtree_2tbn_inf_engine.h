#ifndef JTREE_2TBN_INF_ENGINE_
#define JTREE_2TBN_INF_ENGINE_

#include <dai/factorgraph.h>
#include <dai/jtree.h>
#include "abstract_inf_engine.h"

#include <boost/scoped_ptr.hpp>
#include <vector>

namespace pdbntk {
class JTree2TBNInfEngine : public AbstractInfEngine {
public:
  JTree2TBNInfEngine(const dai::FactorGraph &fg15, const dai::FactorGraph &fg1);
  virtual double EnterEvidence(const mocapy::Sequence &evidence);
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const;
  void Fwd();

private:
  boost::scoped_ptr<dai::JTree> jtree_engine_, jtree_engine1_;
  
  void CollectEvidence();
};
}
#endif // JTREE_2TBN_INF_ENGINE_
