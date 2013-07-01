#include "jtree.h"
#include "gtest/gtest.h"

#include <memory>
#include <vector>
#include <ostream>

TEST(JTreeTest, MaximalCliques) {
  using pdbntk::Node;
  using std::shared_ptr;
  using std::vector;

  pdbntk::FactorGraph ahmm;
  vector<shared_ptr<Node> > nodes; 
  int n = 8;
  for (int i = 0; i < n; i++) {
    std::shared_ptr<Node> node(new Node(i + 1, NULL));
    nodes.push_back(node);
  }
  
  pdbntk::NodeSet ns;
  for (vector<shared_ptr<Node> >::const_iterator it = nodes.begin();
       it != nodes.end(); it++)
    ns.insert(it->get());
}
