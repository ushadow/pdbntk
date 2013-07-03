#include "jtree.h"
#include "cpd/cpd_factory.h"
#include "gtest/gtest.h"
#include "glog/logging.h"

#include <memory>
#include <vector>
#include <ostream>
#include <memory>

TEST(JTreeTest, MaximalCliques) {
  using pdbntk::Node;
  using pdbntk::NodeSet;
  using pdbntk::Factor;
  using pdbntk::CondProbDist;
  using pdbntk::CPDFactory;
  using std::shared_ptr;
  using std::vector;
  using std::unique_ptr;

  vector<Factor> factors;

  unique_ptr<CondProbDist> cpd1(CPDFactory::NewDiscreteCPD(13));
  unique_ptr<CondProbDist> cpd2(CPDFactory::NewDiscreteCPD(44));
  unique_ptr<CondProbDist> cpd3(CPDFactory::NewDiscreteCPD(2));
  unique_ptr<CondProbDist> cpd4(CPDFactory::NewGaussianCPD(9));
  unique_ptr<CondProbDist> cpd5(CPDFactory::NewDiscreteCPD(13));
  unique_ptr<CondProbDist> cpd6(CPDFactory::NewDiscreteCPD(44));

  Node n1(1, cpd1.get());
  Node n2(2, cpd2.get());
  Node n3(3, cpd3.get());
  Node n4(4, cpd4.get());
  Node n5(5, cpd5.get());
  Node n6(6, cpd6.get());
  Node n7(7, cpd3.get());
  Node n8(8, cpd4.get());

  NodeSet ns1(&n1, &n2);
  factors.push_back(Factor(ns1 | &n3));

  NodeSet ns2(&n2, &n5);
  factors.push_back(Factor(ns2 | &n6));

  factors.push_back(Factor(NodeSet(&n6, &n8)));

  NodeSet ns4(&n5, &n6);
  factors.push_back(Factor(ns4 | &n7));

  NodeSet ns5(&n1, &n3);
  factors.push_back(Factor(ns5 | &n5));

  pdbntk::FactorGraph fg(factors);
  pdbntk::JTree jtree(fg);

  for (vector<Factor>::const_iterator it = jtree.Qa.begin(); it != jtree.Qa.end();
       it++)
    LOG(INFO) << (*it); 
}
