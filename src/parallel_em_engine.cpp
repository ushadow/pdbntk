#include "parallel_em_engine.h"
#include "inf_engine.h"

#include "mocapy.h" 

namespace pdbntk {
  
ParallelEmEngine::ParallelEmEngine() : max_iter_(kMaxIter) {
} 

ParallelEmEngine& ParallelEmEngine::set_max_iter(int max_iter) {
  max_iter_ = max_iter;
  return *this;
}

void ParallelEmEngine::Learn(InfEngine *inf_engine,
    const std::vector<mocapy::Sequence> &evidence) {
  bool converged = false;
  int num_iter = 1;
  while (!converged && num_iter <= max_iter_) {
    DoEStep(inf_engine, evidence); 
  }
}

void ParallelEmEngine::DoEStep(InfEngine *inf_engine, 
    std::vector<mocapy::Sequence> evidence) {
  using std::vector;
  using mocapy::ESSBase;
  vector<ESSBase*> ess = inf_engine->GetResetESS();
  #pragma omp parallel
  {
    vector<ESSBase*> ess_private = inf_engine->GetResetESS();
    double loglik = 0;
    #pragma omp for nowait 
    for (int i = 0; i < evidence.size(); i++) {
      loglik += inf_engine->EnterEvidence(evidence[i]);
      for (int j = 0; j < ess_private.size(); j++) {
        ess_private[j]->add_ptv(evidence[i].get_values());
      }
    }
    
    #pragma omp critical (ess_combine)  // Single thread
    for (int k = 0; k < ess.size(); k++) {
      ess[k]->combine(*ess_private[k]);
    }
    
  }
}

void ParallelEmEngine::DoMStep() {

}
}
