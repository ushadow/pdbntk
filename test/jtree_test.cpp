#include "jtree.h"
#include "cpd/cpd_factory.h"
#include "gtest/gtest.h"
#include "glog/logging.h"

#include <memory>
#include <vector>
#include <ostream>

using pdbntk::Node;
using pdbntk::NodeSet;
using pdbntk::Factor;
using pdbntk::CondProbDist;
using pdbntk::CPDFactory;
using pdbntk::JTree;
using pdbntk::Real;
using dai::PropertySet;
using std::shared_ptr;
using std::vector;
using std::unique_ptr;

class JTreeTestF : public testing::Test {
 protected:
  virtual void SetUp() {
    vector<Factor> factors;
    cpd1_.reset(CPDFactory::NewDiscreteCPD(13));
    cpd2_.reset(CPDFactory::NewDiscreteCPD(44));
    cpd3_.reset(CPDFactory::NewDiscreteCPD(2));
    cpd4_.reset(CPDFactory::NewGaussianCPD(9));
    cpd5_.reset(CPDFactory::NewDiscreteCPD(13));
    cpd6_.reset(CPDFactory::NewDiscreteCPD(44));

    n1_.reset(new Node(1, cpd1_.get(), true));
    n2_.reset(new Node(2, cpd2_.get()));
    n3_.reset(new Node(3, cpd3_.get(), true));
    n4_.reset(new Node(4, cpd4_.get(), true));
    n5_.reset(new Node(5, cpd5_.get(), true));
    n6_.reset(new Node(6, cpd6_.get()));
    n7_.reset(new Node(7, cpd3_.get(), true));
    n8_.reset(new Node(8, cpd4_.get(), true));

    NodeSet ns1(n1_.get(), n2_.get());
    factors.push_back(Factor(ns1 | n3_.get()));

    NodeSet ns2(n2_.get(), n5_.get());
    factors.push_back(Factor(ns2 | n6_.get()));

    factors.push_back(Factor(NodeSet(n6_.get(), n8_.get())));

    NodeSet ns4(n5_.get(), n6_.get());
    factors.push_back(Factor(ns4 | n7_.get()));

    NodeSet ns5(n1_.get(), n3_.get());
    factors.push_back(Factor(ns5 | n5_.get()));

    pdbntk::FactorGraph fg(factors);

    PropertySet ps;
    ps.set("root", ns4 | n7_.get());

    jtree_.reset(new JTree(fg, ps));
  } 

  unique_ptr<JTree> jtree_;
  unique_ptr<CondProbDist> cpd1_, cpd2_, cpd3_, cpd4_, cpd5_, cpd6_;
  unique_ptr<Node> n1_, n2_, n3_, n4_, n5_, n6_, n7_, n8_;
};


bool IsEqFactor(const uint expected[], Factor f) {
  bool ret = true;
  const vector<Node*> &ns = f.nodes().elements();
  for (uint i = 0; i < ns.size(); i++) {
    ret = ret && (expected[i] == ns[i]->index()); 
  }
  return ret;
}

TEST_F(JTreeTestF, MaximalCliques) {
  uint qa0[4] = {1, 2, 3, 5};
  uint qa1[3] = {2, 5, 6};
  uint qa2[3] = {5, 6, 7};
  uint qa3[2] = {6, 8};

  EXPECT_PRED2(IsEqFactor, qa0, jtree_->Qa[0]);
  EXPECT_PRED2(IsEqFactor, qa1, jtree_->Qa[1]);
  EXPECT_PRED2(IsEqFactor, qa2, jtree_->Qa[2]);
  EXPECT_PRED2(IsEqFactor, qa3, jtree_->Qa[3]);

  uint qb0[2] = {5, 6};
  uint qb1[1] = {6};
  uint qb2[2] = {2, 5};

  EXPECT_PRED2(IsEqFactor, qb0, jtree_->Qb[0]);
  EXPECT_PRED2(IsEqFactor, qb1, jtree_->Qb[1]);
  EXPECT_PRED2(IsEqFactor, qb2, jtree_->Qb[2]);
 
  size_t root = jtree_->RTree[0].first; 
  EXPECT_EQ(2, root);
}

TEST_F(JTreeTestF, Properties) {
  PropertySet ps = jtree_->getProperties();
  NodeSet ns = ps.getAs<NodeSet>("root");
  NodeSet expected(n5_.get(), n6_.get());
  EXPECT_EQ(expected | n7_.get(), ns);
}

TEST_F(JTreeTestF, Clamp) {
  vector<Real> x(1, 1);
  jtree_->clamp(1, x); 
}

TEST(JTreeTest, DifferentFG) {
  vector<Factor> factors;

  unique_ptr<CondProbDist> cpd1(CPDFactory::NewDiscreteCPD(13));
  unique_ptr<CondProbDist> cpd2(CPDFactory::NewDiscreteCPD(44));
  unique_ptr<CondProbDist> cpd3(CPDFactory::NewDiscreteCPD(2));
  unique_ptr<CondProbDist> cpd4(CPDFactory::NewGaussianCPD(9));

  Node n1(1, cpd1.get(), true);
  Node n2(2, cpd2.get());
  Node n3(3, cpd3.get(), true);
  Node n4(4, cpd4.get(), true);

  NodeSet ns1(&n1, &n2);
  factors.push_back(Factor(ns1 | &n3));

  factors.push_back(Factor(NodeSet(&n2, &n4)));

  pdbntk::FactorGraph fg(factors);
  
  PropertySet ps;
  ps.set("root", ns1 | &n3);

  JTree jtree(fg, ps);
  uint qa0[3] = {1, 2, 3};
  uint qa1[2] = {2, 4};

  EXPECT_PRED2(IsEqFactor, qa0, jtree.Qa[0]);
  EXPECT_PRED2(IsEqFactor, qa1, jtree.Qa[1]);

  uint qb0[1] = {2};
  EXPECT_PRED2(IsEqFactor, qb0, jtree.Qb[0]);

  EXPECT_EQ(0, jtree.RTree[0].first);
}
