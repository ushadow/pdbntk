#include "parallel_em_engine.h"
#include "abstract_inf_engine.h"

#include "mocapy.h" 

#include <omp.h>
#include <iostream>

namespace pdbntk {
  
ParallelEmEngine::ParallelEmEngine(int max_iter) : max_iter_(max_iter) {
} 

ParallelEmEngine& ParallelEmEngine::set_max_iter(int max_iter) {
  max_iter_ = max_iter;
  return *this;
}

double ParallelEmEngine::Iterate(AbstractInfEngine *inf_engine,
    std::vector<mocapy::Sequence> *evidence) {
  bool converged = false;
  int num_iter = 1;
  double loglik = 0;
  while (!converged && num_iter <= max_iter_) {
    loglik = DoEStep(inf_engine, evidence); 
    num_iter++;
  }
  return loglik;
}

double ParallelEmEngine::DoEStep(AbstractInfEngine *inf_engine, 
    std::vector<mocapy::Sequence> *evidence) {
  using std::vector;
  using mocapy::ESSBase;
  vector<ESSBase*> ess = inf_engine->GetResetESS();
  double loglik = 0;
  #pragma omp parallel
  {
    vector<ESSBase*> ess_private = inf_engine->GetResetESS();
    double loglik_private = 0;
    std::cout << "thread " << omp_get_thread_num() << std::endl;
    #pragma omp for nowait 
    for (int i = 0; i < evidence->size(); i++) {
      loglik_private += inf_engine->EnterEvidence((*evidence)[i]);
      for (int j = 0; j < ess_private.size(); j++) {
        ess_private[j]->add_ptv((*evidence)[i].get_values());
      }
    }
    
    #pragma omp critical (ess_combine)  // Single thread
    {
      for (int k = 0; k < ess.size(); k++) {
        ess[k]->combine(*ess_private[k]);
      }
      loglik += loglik_private;
    }
      
  }
  return loglik; 
}

void ParallelEmEngine::DoMStep() {

}
}
