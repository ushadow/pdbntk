#ifndef PDBNTK_UTILS_H_
#define PDBNTK_UTILS_H_

#include "mocapy.h"
#include <dai/smallset.h>
#include <vector>
#include <glog/logging.h>

namespace pdbntk {

class Node;

typedef float Real;
typedef dai::SmallSet<Node*> NodeSet;

std::vector<Node*> MkSliceAndHalfBNet(const mocapy::DBN &dbn, 
    const std::vector<mocapy::NodeID> &interface);
}

#endif // PDBNTK_UTILS_H_
