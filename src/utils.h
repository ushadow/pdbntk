#ifndef UTILS_H_
#define UtILS_H_

#include "mocapy.h"

#include <vector>
#include <glog/logging.h>

namespace pdbntk {

class Node;

template<class ESS, class Densities>
mocapy::ChildNode<ESS, Densities>* Clone(const mocapy::ChildNode<ESS, Densities> &node) {
  mocapy::eNodeType type = node.get_node_type();
  switch (type) {
    case mocapy::DISCRETE:
      break;
    case mocapy::GAUSSIAN:
      break;
    default:
      LOG(FATAL) << "Node type not supported: " << type ;
  }
}
}

#endif // UTILS_H_
