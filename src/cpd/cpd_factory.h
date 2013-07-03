// The CPDFactory creates the desired objects for you. You can also do this yourself, but it should be easier using the NodeFactory.

#ifndef PDBNTK_CPD_FACTORY_H_
#define PDBNTK_CPD_FACTORY_H_ 

#include "cond_prob_dist.h"

namespace pdbntk {

/// Factory for creating ConProbDist objects.
class CPDFactory {
public:
	virtual ~CPDFactory() {};

	static CondProbDist* NewDiscreteCPD(uint node_size);
	static CondProbDist* NewGaussianCPD(uint dimension);

private:
	CPDFactory() {};
  CPDFactory(const CPDFactory&);
  CPDFactory& operator=(const CPDFactory&);
};
}

#endif // PDBNTK_CPD_FACTORY_H_
