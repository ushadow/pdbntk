#ifndef JTREE_INF_ENGINE_H_ 
#define JTREE_INF_ENGINE_H_ 

#include "mocapy.h"
#include "abstract_inf_engine.h"

#include <vector>

namespace pdbntk {

class Node;

/// Junction tree inference engine for a Bayesian network.
class JTreeInfEngine : public AbstractInfEngine {
public:
  /** \param bnet a directed graph.
   */
  JTreeInfEngine(std::vector<Node*> *bnet);
  virtual double EnterEvidence(const mocapy::Sequence &evidence);
  virtual std::vector<mocapy::ESSBase*> GetResetESS() const;
};
}
#endif // JTREE_INF_ENGINE_H_
