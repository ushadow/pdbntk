#ifndef PARALLEL_EM_ENGINE_H_
#define PARALLEL_EM_ENGEIN_H_

#include "mocapy.h"

#include <vector>

namespace pdbntk {

class AbstractInfEngine;

class ParallelEmEngine {
public:
  ParallelEmEngine(int max_iter);
  
  ParallelEmEngine& set_max_iter(int max_iter); 
  
  double Iterate(AbstractInfEngine *inf_engine, std::vector<mocapy::Sequence> *evidence);
 
private: 
  static const int kMaxIter = 100;
  int max_iter_;

  double DoEStep(AbstractInfEngine *inf_engine, std::vector<mocapy::Sequence> *evidence);
  void DoMStep();
};  
}

#endif // PARALLEL_EM_ENGINE_H_
