#include "cpd_factory.h"

#include "discrete/discretedensities.h"
#include "discrete/discreteess.h"
#include "gaussian/gaussiandensities.h"
#include "gaussian/gaussianess.h"

namespace pdbntk {

CondProbDist* CPDFactory::NewDiscreteCPD(uint node_size) {
  return new CondProbDist(DISCRETE, new mocapy::DiscreteESS(),
                         new mocapy::DiscreteDensities(node_size));
}

CondProbDist* CPDFactory::NewGaussianCPD(uint dimension) {
  return new CondProbDist(GAUSSIAN, new mocapy::GaussianESS(),
                         new mocapy::GaussianDensities(dimension));
}

}
