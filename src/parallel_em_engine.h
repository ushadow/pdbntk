#ifndef PARALLEL_EM_ENGINE_H_
#define PARALLEL_EM_ENGEIN_H_

#include "mocapy.h"

#include<vector>

namespace pdbntk {

class InfEngine;

class ParallelEmEngine {
public:
  ParallelEmEngine();
  
  ParallelEmEngine& set_max_iter(int max_inter); 
  
  void Learn(InfEngine *inf_engine,
      const std::vector<mocapy::Sequence> &evidence);
 
private: 
  static const int kMaxIter = 100;
  int max_iter_;

  void DoEStep(InfEngine *inf_engine, std::vector<mocapy::Sequence> evidence);
  void DoMStep();
};  
}

#endif // PARALLEL_EM_ENGINE_H_
