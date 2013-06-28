#include "region_graph.h"
#include "factor.h"
#include "cluster_graph.h"

#include <algorithm>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
#include <glog/logging.h>

namespace pdbntk {

void RegionGraph::construct( const FactorGraph &fg, const std::vector<NodeSet> &ors, const std::vector<Region> &irs, const std::vector<std::pair<size_t,size_t> > &edges ) {
    // Copy factor graph structure
    FactorGraph::operator=( fg );

    // Copy inner regions
    _IRs = irs;

    // Construct outer regions (giving them counting number 1.0)
    _ORs.clear();
    _ORs.reserve( ors.size() );
    bforeach( const NodeSet &alpha, ors )
        _ORs.push_back( FRegion(Factor(alpha), 1.0) );

    // For each factor, find an outer region that subsumes that factor.
    // Then, multiply the outer region with that factor.
    _fac2OR.clear();
    _fac2OR.reserve( nrFactors() );
    for( size_t I = 0; I < nrFactors(); I++ ) {
        size_t alpha;
        for( alpha = 0; alpha < nrORs(); alpha++ )
            if( OR(alpha).nodes() >> factor(I).nodes() ) {
                _fac2OR.push_back( alpha );
                break;
            }
        DAI_ASSERT( alpha != nrORs() );
    }
    recomputeORs();

    // Create bipartite graph
    _G.construct( nrORs(), nrIRs(), edges.begin(), edges.end() );
}


void RegionGraph::constructCVM( const FactorGraph &fg, const std::vector<NodeSet> &cl) {
  using std::pair;

  DLOG(INFO) << "constructCVM called (" << fg.nrVars() << " vars, " << fg.nrFactors() << " facs, " << cl.size() << " clusters)";

  // Retain only maximal clusters
  DLOG(INFO) << "  Constructing ClusterGraph";
  ClusterGraph cg( cl );
  DLOG(INFO) << "  Erasing non-maximal clusters";
  cg.eraseNonMaximal();

  // Create inner regions - first pass
  DLOG(INFO) << "  Creating inner regions (first pass)";
  std::set<NodeSet> betas;
  for( size_t alpha = 0; alpha < cg.nrClusters(); alpha++ )
    for( size_t alpha2 = alpha; (++alpha2) != cg.nrClusters(); ) {
      NodeSet intersection = cg.cluster(alpha) & cg.cluster(alpha2);
      if( intersection.size() > 0 )
        betas.insert( intersection );
    }

  // Create inner regions - subsequent passes
  DLOG(INFO) << "  Creating inner regions (next passes)";
  std::set<NodeSet> new_betas;
  do {
    new_betas.clear();
    for (std::set<NodeSet>::const_iterator gamma = betas.begin(); gamma != betas.end(); gamma++ )
      for (std::set<NodeSet>::const_iterator gamma2 = gamma; (++gamma2) != betas.end(); ) {
        NodeSet intersection = (*gamma) & (*gamma2);
        if( (intersection.size() > 0) && (betas.count(intersection) == 0) )
          new_betas.insert( intersection );
      }
    betas.insert(new_betas.begin(), new_betas.end());
  } while( new_betas.size() );

  // Create inner regions - final phase
  DLOG(INFO) << "  Creating inner regions (final phase)";
  std::vector<Region> irs;
  irs.reserve( betas.size() );
  for (std::set<NodeSet>::const_iterator beta = betas.begin(); beta != betas.end(); beta++ )
    irs.push_back( Region(*beta,0.0) );

  // Create edges
  DLOG(INFO) << "  Creating edges";
  std::vector<std::pair<size_t,size_t> > edges;
  for( size_t beta = 0; beta < irs.size(); beta++ )
    for( size_t alpha = 0; alpha < cg.nrClusters(); alpha++ )
      if( cg.cluster(alpha) >> irs[beta] )
        edges.push_back( pair<size_t,size_t>(alpha,beta) );

  // Construct region graph
  DLOG(INFO) << "  Constructing region graph";
  construct( fg, cg.clusters(), irs, edges );

  // Calculate counting numbers
  DLOG(INFO) << "  Calculating counting numbers";
  calcCVMCountingNumbers();

  DLOG(INFO) << "Done.";
}


void RegionGraph::calcCVMCountingNumbers() {
  using std::vector;
  // Calculates counting numbers of inner regions based upon counting numbers of outer regions

  vector<vector<size_t> > ancestors(nrIRs());
  boost::dynamic_bitset<> assigned(nrIRs());
  for( size_t beta = 0; beta < nrIRs(); beta++ ) {
    IR(beta).c() = 0.0;
    for( size_t beta2 = 0; beta2 < nrIRs(); beta2++ )
      if( (beta2 != beta) && IR(beta2) >> IR(beta) )
        ancestors[beta].push_back(beta2);
  }

  bool new_counting;
  do {
    new_counting = false;
    for( size_t beta = 0; beta < nrIRs(); beta++ ) {
      if( !assigned[beta] ) {
        bool has_unassigned_ancestor = false;
        for( vector<size_t>::const_iterator beta2 = ancestors[beta].begin(); (beta2 != ancestors[beta].end()) && !has_unassigned_ancestor; beta2++ )
          if( !assigned[*beta2] )
            has_unassigned_ancestor = true;
        if( !has_unassigned_ancestor ) {
          Real c = 1.0;
          bforeach( const dai::Neighbor &alpha, nbIR(beta) )
            c -= OR(alpha).c();
          for( vector<size_t>::const_iterator beta2 = ancestors[beta].begin(); beta2 != ancestors[beta].end(); beta2++ )
            c -= IR(*beta2).c();
          IR(beta).c() = c;
          assigned.set(beta, true);
          new_counting = true;
        }
      }
    }
  } while( new_counting );
}


bool RegionGraph::checkCountingNumbers() const {
  // Checks whether the counting numbers satisfy the fundamental relation

  bool all_valid = true;
  for (std::vector<Node*>::const_iterator n = nodes().begin(); n != nodes().end(); n++ ) {
    Real c_n = 0.0;
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
      if( OR(alpha).nodes().contains( *n ) )
        c_n += OR(alpha).c();
    for( size_t beta = 0; beta < nrIRs(); beta++ )
      if( IR(beta).contains( *n ) )
        c_n += IR(beta).c();
    if( fabs(c_n - 1.0) > 1e-15 ) {
      all_valid = false;
      LOG(FATAL) << "WARNING: counting numbers do not satisfy relation for " << *n << "(c_n = " << c_n << ").";
    }
  }

  return all_valid;
}


void RegionGraph::recomputeORs() {
}


void RegionGraph::recomputeORs( const NodeSet &ns ) {
}


void RegionGraph::recomputeOR( size_t I ) {
}


/// Send RegionGraph to output stream
std::ostream & operator << (std::ostream & os, const RegionGraph & rg) {
  using std::endl;
  os << "digraph RegionGraph {" << endl;
  os << "node[shape=box];" << endl;
  for( size_t alpha = 0; alpha < rg.nrORs(); alpha++ )
    os << "\ta" << alpha << " [label=\"a" << alpha << ": " << rg.OR(alpha).nodes() << ", c=" << rg.OR(alpha).c() << "\"];" << endl;
  os << "node[shape=ellipse];" << endl;
  for( size_t beta = 0; beta < rg.nrIRs(); beta++ )
    os << "\tb" << beta << " [label=\"b" << beta << ": " << (NodeSet)rg.IR(beta) << ", c=" << rg.IR(beta).c() << "\"];" << endl;
  for( size_t alpha = 0; alpha < rg.nrORs(); alpha++ )
    bforeach( const dai::Neighbor &beta, rg.nbOR(alpha) )
      os << "\ta" << alpha << " -> b" << beta << ";" << endl;
  os << "}" << endl;
  return os;
}


} // end of namespace dai
