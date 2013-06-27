#ifndef UTILS_H_
#define UtILS_H_

#include "mocapy.h"

#include <vector>
#include <glog/logging.h>

namespace pdbntk {

class Node;

std::vector<Node*> MkSliceAndHalfBNet(const mocapy::DBN &dbn, 
    const std::vector<mocapy::NodeID> &interface);
}

#endif // UTILS_H_
