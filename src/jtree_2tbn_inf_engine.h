#ifndef JTREE_2TBN_INF_ENGINE_
#define JTREE_2TBN_INF_ENGINE_

#include "jtree.h"
#include "abstract_inf_engine.h"
#include "graph/factor_graph.h"
#include "evidence.h"

#include <memory> 
#include <vector>

namespace pdbntk {
class JTree2TBNInfEngine : public AbstractInfEngine {
public:
  /// \param fg15 1.5 slice in the DBN. Slice 0.5 is the interface nodes
  /// from slice 1. The interface nodes in slice 0.5 should be in a factor
  /// (cluster) and the interface nodes in slice 2 should also be in a factor.
  /// The nodes in the interface for both slice 1 and 2 must be in a factor.
  JTree2TBNInfEngine(const FactorGraph &fg1, const FactorGraph &fg15,
      const NodeSet &interface1, const NodeSet &interface2);
  virtual double EnterEvidence(const mocapy::Sequence &evidence);
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const;

  /// Forward step
  /// \param o observation at time t.
  void Fwd(const Evidence::Observation &o, int t);
  /// Forward pass for slice 1
  void Fwd1(const Evidence::Observation &o, int t);

private:
  std::unique_ptr<JTree> jtree_engine_, jtree_engine1_;
  
  void CollectEvidence();
};
}
#endif // JTREE_2TBN_INF_ENGINE_
