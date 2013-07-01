#include "jtree.h"
#include "gtest/gtest.h"
#include "glog/logging.h"

#include <memory>
#include <vector>
#include <ostream>

TEST(JTreeTest, MaximalCliques) {
  using pdbntk::Node;
  using pdbntk::NodeSet;
  using pdbntk::Factor;
  using std::shared_ptr;
  using std::vector;

  vector<Factor> factors;

  vector<shared_ptr<Node> > nodes; 
  int n = 8;
  for (int i = 0; i < n; i++) {
    std::shared_ptr<Node> node(new Node(i, NULL));
    nodes.push_back(node);
  }
  
  NodeSet ns1(nodes[0].get(), nodes[1].get());
  factors.push_back(Factor(ns1 | nodes[2].get()));

  NodeSet ns2(nodes[1].get(), nodes[4].get());
  factors.push_back(Factor(ns2 | nodes[5].get()));

  factors.push_back(Factor(NodeSet(nodes[5].get(), nodes[7].get())));

  NodeSet ns4(nodes[4].get(), nodes[5].get());
  factors.push_back(Factor(ns4 | nodes[6].get()));

  pdbntk::FactorGraph fg(factors);
  pdbntk::JTree jtree(fg);

  for (vector<Factor>::const_iterator it = jtree.Qa.begin(); it != jtree.Qa.end();
       it++)
    LOG(INFO) << (*it); 
}
