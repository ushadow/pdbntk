#include "pdbntk.h"
#include "mocapy.h" 

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <utility>

int main(void) {
  using mocapy::Node;
  using mocapy::NodeFactory;
  using mocapy::DBN;
  using mocapy::MDArray;
  using mocapy::eMISMASK;
  using mocapy::vec;
  using std::vector;
  using std::pair;
  
  uint hsize = 2;
  uint osize = 2;
  bool init_random = false;
  bool use_shrinkage = true;

  // Number of training sequences.
  int n = 40000;

  // Sequence lengths
  int t = 100;

  mocapy::CPD th0_cpd;
  th0_cpd.set_shape(2); 
  th0_cpd.set_values(mocapy::vec(0.1, 0.9));

  CPD th1_cpd;
  th1_cpd.set_shape(2, 2); 
  th1_cpd.set_values(mocapy::vec(0.95, 0.05, 0.1, 0.9));

  MDArray<double> means;
  means.set_shape(2, 2);
  means.set_values(vec(0.0, 0.0, 10.0, 10.0));

  Node* th0 = NodeFactory::new_discrete_node(hsize, "th0", init_random, 
                                             th0_cpd);
  Node* th1 = NodeFactory::new_discrete_node(hsize, "th1", init_random,
                                             th1_cpd);
  Node* to0 = NodeFactory::new_gaussian_node(osize, "to0", false, use_shrinkage,
                                             means); 

  DBN tdbn;
  tdbn.set_slices(vec(th0, to0), vec(th1, to0));

  tdbn.add_intra("th0", "to0");
  tdbn.add_inter("th0", "th1");
  tdbn.construct();
  
  Node* mh0 = NodeFactory::new_discrete_node(hsize, "mh0", init_random, CPD(),
                                             th0, true);
  Node* mh1 = NodeFactory::new_discrete_node(hsize, "mh1", init_random);
  Node* mo0 = NodeFactory::new_gaussian_node(osize, "mo0", true, use_shrinkage);

  DBN mdbn;  
  mdbn.set_slices(vec(mh0, mo0), vec(mh1, mo0));

  mdbn.add_intra("mh0", "mo0");
  mdbn.add_inter("mh0", "mh1");
  mdbn.construct();

  vector<mocapy::Sequence> seq_list;
  vector<MDArray<eMISMASK> > mismask_list;
  
  MDArray<eMISMASK> mismask;
  mismask.repeat(t, vec(mocapy::MOCAPY_HIDDEN, mocapy::MOCAPY_OBSERVED));

  for (int i = 0; i < n; i++) {
    mocapy::Sequence seq(vec(2, t));
    seq_list.push_back(seq);
    mismask_list.push_back(mismask);
  }

  pdbntk::DummyInfEngine inf_engine;
  pdbntk::ParallelEmEngine em(10);
  double loglik = em.Iterate(&inf_engine, &seq_list); 
  std::cout << "loglik = " << loglik << std::endl;
  return 0;     
}
