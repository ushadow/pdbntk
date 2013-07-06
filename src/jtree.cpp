#include "jtree.h"
#include "graph/cluster_graph.h"
#include "graph/node.h"
#include "graph/factor.h"
#include "graph/factor_graph.h"

#include <iostream>
#include <stack>
#include <glog/logging.h>

namespace pdbntk {

using std::vector;
using dai::DEdge;
using dai::UEdge;
using dai::WeightedGraph;
using dai::RootedTree;
using dai::GraphEL;

void JTree::setProperties(const dai::PropertySet &opts) {
  using dai::PropertySet;

  // Default update type is HUGIN.
  if (opts.hasKey("updates"))
    props_.updates = opts.getStringAs<Properties::UpdateType>("updates");
  else
    props_.updates = Properties::UpdateType::HUGIN;
      
  if (opts.hasKey("inference"))
    props_.inference = opts.getStringAs<Properties::InfType>("inference");
  else
    props_.inference = Properties::InfType::SUMPROD;
  
  if (opts.hasKey("heuristic"))
    props_.heuristic = opts.getStringAs<Properties::HeuristicType>("heuristic");
  else
    props_.heuristic = Properties::HeuristicType::WEIGHTEDMINFILL;

  if (opts.hasKey("maxmem"))
    props_.maxmem = opts.getStringAs<size_t>("maxmem");
  else
    props_.maxmem = 0;

  if (opts.hasKey("root"))
    props_.root = opts.getAs<NodeSet>("root");
}

dai::PropertySet JTree::getProperties() const {
  dai::PropertySet opts;
  opts.set("updates", props_.updates);
  opts.set("inference", props_.inference);
  opts.set("heuristic", props_.heuristic);
  opts.set("maxmem", props_.maxmem);
  opts.set("root", props_.root);
  return opts;
}


std::string JTree::printProperties() const {
  std::stringstream s(std::stringstream::out );
  s << "[";
  s << "updates=" << props_.updates << ",";
  s << "heuristic=" << props_.heuristic << ",";
  s << "inference=" << props_.inference << ",";
  s << "maxmem=" << props_.maxmem << "]";
  return s.str();
}


JTree::JTree(const FactorGraph &fg, const dai::PropertySet &opts, bool automatic)
    : DAIAlgRG(), mes_(), logz_(), RTree(), Qa(), Qb(), props_() {
  setProperties(opts);

  if(automatic) {
    // Create ClusterGraph which contains maximal factors as clusters
    ClusterGraph _cg(fg, true);
    DLOG(INFO) << "Initial clusters: " << _cg;

    // Use heuristic to guess optimal elimination sequence
    greedyVariableElimination::eliminationCostFunction ec(NULL);
    switch( (size_t)props_.heuristic ) {
      case Properties::HeuristicType::MINNEIGHBORS:
        ec = eliminationCost_MinNeighbors;
        break;
      case Properties::HeuristicType::MINWEIGHT:
        ec = eliminationCost_MinWeight;
        break;
      case Properties::HeuristicType::MINFILL:
        ec = eliminationCost_MinFill;
        break;
      case Properties::HeuristicType::WEIGHTEDMINFILL:
        ec = eliminationCost_WeightedMinFill;
        break;
      default:
        DAI_THROW(UNKNOWN_ENUM_VALUE);
    }
    size_t fudge = 6; // this yields a rough estimate of the memory needed (for some reason not yet clearly understood)
    std::vector<NodeSet> ElimVec = _cg.VarElim(greedyVariableElimination(ec), props_.maxmem / (sizeof(Real) * fudge)).eraseNonMaximal().clusters();
    DLOG(INFO) << "VarElim result: " << ElimVec; 
    // Generate the junction tree corresponding to the elimination sequence
    GenerateJT(fg, ElimVec);
  }
}


void JTree::construct(const FactorGraph &fg, const vector<NodeSet> &cl,
                      bool verify) {
  // Copy the factor graph
  FactorGraph::operator=(fg);

  // Construct a weighted graph (each edge is weighted with the cardinality
  // of the intersection of the nodes, where the nodes are the elements of cl).
  WeightedGraph<int> JuncGraph;
  // Start by connecting all clusters with cluster zero, and weight zero,
  // in order to get a connected weighted graph
  for (size_t i = 1; i < cl.size(); i++)
    JuncGraph[UEdge(i, 0)] = 0;
  for (size_t i = 0; i < cl.size(); i++) {
    for (size_t j = i + 1; j < cl.size(); j++) {
      size_t w = (cl[i] & cl[j]).size();
      if (w)
        JuncGraph[UEdge(i, j)] = w;
    }
  }
  DLOG(INFO) << "Weightedgraph: " << JuncGraph;

  // Construct maximal spanning tree using Prim's algorithm
  RTree = MaxSpanningTree(JuncGraph, true);
  size_t root = FindRoot(cl);
  if (root < cl.size());
    RTree = RootedTree(GraphEL(RTree.begin(), RTree.end()), root);
  DLOG(INFO) << "Spanning tree: " << RTree;
  DAI_DEBASSERT( RTree.size() == cl.size() - 1 );

  // Construct corresponding region graph

  // Create outer regions
  _ORs.clear();
  _ORs.reserve( cl.size() );
  for( size_t i = 0; i < cl.size(); i++ )
    _ORs.push_back(FRegion(Factor(cl[i]), 1.0 ));

  // For each factor, find an outer region that subsumes that factor.
  // Then, multiply the outer region with that factor.
  _fac2OR.clear();
  _fac2OR.resize(nrFactors(), -1U);
  for (size_t I = 0; I < nrFactors(); I++) {
    size_t alpha;
    for (alpha = 0; alpha < nrORs(); alpha++)
      if (OR(alpha).nodes() >> factor(I).nodes()) {
        _fac2OR[I] = alpha;
        break;
      }
    if (verify)
      DAI_ASSERT(alpha != nrORs());
  }
  RecomputeORs();

  // Create inner regions and edges
  // Inner regions are analogous to separators.
  _IRs.clear();
  _IRs.reserve(RTree.size());
  vector<dai::Edge> edges;
  edges.reserve(2 * RTree.size());
  for( size_t i = 0; i < RTree.size(); i++ ) {
    edges.push_back(dai::Edge(RTree[i].first, nrIRs()));
    edges.push_back(dai::Edge(RTree[i].second, nrIRs()));
    // inner clusters have counting number -1, except if they are empty
    NodeSet intersection = cl[RTree[i].first] & cl[RTree[i].second];
    _IRs.push_back(Region(intersection, intersection.size() ? -1.0 : 0.0));
  }

  // create bipartite graph
  _G.construct(nrORs(), nrIRs(), edges.begin(), edges.end());

  // Check counting numbers
#ifdef DAI_DEBUG
  checkCountingNumbers();
#endif

  // Create beliefs
  Qa.clear();
  Qa.reserve(nrORs());
  for (size_t alpha = 0; alpha < nrORs(); alpha++)
    Qa.push_back(OR(alpha));

  Qb.clear();
  Qb.reserve(nrIRs());
  for (size_t beta = 0; beta < nrIRs(); beta++)
    Qb.push_back(Factor(IR(beta)));
}


void JTree::GenerateJT(const FactorGraph &fg, const std::vector<NodeSet> &cl) {
  construct(fg, cl, true);

  // Create messages
  mes_.clear();
  mes_.reserve(nrORs());
  for (size_t alpha = 0; alpha < nrORs(); alpha++) {
    mes_.push_back(std::vector<Factor>());
    mes_[alpha].reserve(nbOR(alpha).size());
    bforeach(const dai::Neighbor &beta, nbOR(alpha))
      mes_[alpha].push_back(Factor(IR(beta)));
  }

  DLOG(INFO) << "Regiongraph generated by JTree::GenerateJT: " << *this;
}


Factor JTree::belief(const NodeSet &vs) const {
  using std::vector;
  vector<Factor>::const_iterator beta;
  for( beta = Qb.begin(); beta != Qb.end(); beta++ )
    if( beta->nodes() >> vs )
      break;
  if( beta != Qb.end() ) {
    if( props_.inference == Properties::InfType::SUMPROD )
      return( beta->marginal(vs) );
    else
      return( beta->maxMarginal(vs) );
  } else {
    vector<Factor>::const_iterator alpha;
    for( alpha = Qa.begin(); alpha != Qa.end(); alpha++ )
      if( alpha->nodes() >> vs )
        break;
    if( alpha == Qa.end() ) {
      DAI_THROW(BELIEF_NOT_AVAILABLE);
      return Factor();
    } else {
      if( props_.inference == Properties::InfType::SUMPROD )
        return( alpha->marginal(vs) );
      else
        return( alpha->maxMarginal(vs) );
    }
  }
}


std::vector<Factor> JTree::beliefs() const {
  std::vector<Factor> result;
  for( size_t beta = 0; beta < nrIRs(); beta++ )
    result.push_back( Qb[beta] );
  for( size_t alpha = 0; alpha < nrORs(); alpha++ )
    result.push_back( Qa[alpha] );
  return result;
}


void JTree::runHUGIN() {
  for( size_t alpha = 0; alpha < nrORs(); alpha++ )
    Qa[alpha] = OR(alpha);

  // CollectEvidence
  logz_ = 0.0;
  for( size_t i = RTree.size(); (i--) != 0; ) {
    //      Make outer region RTree[i].first consistent with outer region RTree[i].second
    //      IR(i) = seperator OR(RTree[i].first) && OR(RTree[i].second)
    Factor new_Qb;
    if( props_.inference == Properties::InfType::SUMPROD )
      new_Qb = Qa[RTree[i].second].marginal( IR( i ), false );
    else
      new_Qb = Qa[RTree[i].second].maxMarginal( IR( i ), false );

    logz_ += log(new_Qb.normalize());
    Qa[RTree[i].first] *= new_Qb / Qb[i];
    Qb[i] = new_Qb;
  }
  if( RTree.empty() )
    logz_ += log(Qa[0].normalize() );
  else
    logz_ += log(Qa[RTree[0].first].normalize());

  // DistributeEvidence
  for( size_t i = 0; i < RTree.size(); i++ ) {
    //      Make outer region RTree[i].second consistent with outer region RTree[i].first
    //      IR(i) = seperator OR(RTree[i].first) && OR(RTree[i].second)
    Factor new_Qb;
    if( props_.inference == Properties::InfType::SUMPROD )
      new_Qb = Qa[RTree[i].first].marginal( IR( i ) );
    else
      new_Qb = Qa[RTree[i].first].maxMarginal( IR( i ) );

    Qa[RTree[i].second] *= new_Qb / Qb[i];
    Qb[i] = new_Qb;
  }

  // Normalize
  for( size_t alpha = 0; alpha < nrORs(); alpha++ )
    Qa[alpha].normalize();
}


void JTree::runShaferShenoy() {
  // First pass
  logz_ = 0.0;
  for( size_t e = nrIRs(); (e--) != 0; ) {
    // send a message from RTree[e].second to RTree[e].first
    // or, actually, from the seperator IR(e) to RTree[e].first

    size_t i = nbIR(e)[1].node; // = RTree[e].second
    size_t j = nbIR(e)[0].node; // = RTree[e].first
    size_t _e = nbIR(e)[0].dual;

    Factor msg = OR(i);
    bforeach( const dai::Neighbor &k, nbOR(i) )
      if( k != e )
        msg *= message(i, k.iter);
    if( props_.inference == Properties::InfType::SUMPROD )
      message( j, _e ) = msg.marginal( IR(e), false );
    else
      message( j, _e ) = msg.maxMarginal( IR(e), false );
    logz_ += log( message(j,_e).normalize() );
  }

  // Second pass
  for( size_t e = 0; e < nrIRs(); e++ ) {
    size_t i = nbIR(e)[0].node; // = RTree[e].first
    size_t j = nbIR(e)[1].node; // = RTree[e].second
    size_t _e = nbIR(e)[1].dual;

    Factor msg = OR(i);
    bforeach( const dai::Neighbor &k, nbOR(i) )
      if( k != e )
        msg *= message(i, k.iter);
    if( props_.inference == Properties::InfType::SUMPROD )
      message( j, _e ) = msg.marginal( IR(e) );
    else
      message( j, _e ) = msg.maxMarginal( IR(e) );
  }

  // Calculate beliefs
  for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
    Factor piet = OR(alpha);
    bforeach( const dai::Neighbor &k, nbOR(alpha) )
      piet *= message( alpha, k.iter );
    if( nrIRs() == 0 ) {
      logz_ += log( piet.normalize() );
      Qa[alpha] = piet;
    } else if( alpha == nbIR(0)[0].node /*RTree[0].first*/ ) {
      logz_ += log( piet.normalize() );
      Qa[alpha] = piet;
    } else
      Qa[alpha] = piet.normalized();
  }

  // Only for logZ (and for belief)...
  for( size_t beta = 0; beta < nrIRs(); beta++ ) {
    if( props_.inference == Properties::InfType::SUMPROD )
      Qb[beta] = Qa[nbIR(beta)[0].node].marginal( IR(beta) );
    else
      Qb[beta] = Qa[nbIR(beta)[0].node].maxMarginal( IR(beta) );
  }
}


Real JTree::run() {
  if( props_.updates == Properties::UpdateType::HUGIN )
    runHUGIN();
  else if( props_.updates == Properties::UpdateType::SHSH )
    runShaferShenoy();
  return 0.0;
}


Real JTree::logZ() const {
  /*    Real s = 0.0;
        for( size_t beta = 0; beta < nrIRs(); beta++ )
        s += IR(beta).c() * Qb[beta].entropy();
        for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        s += OR(alpha).c() * Qa[alpha].entropy();
        s += (OR(alpha).log(true) * Qa[alpha]).sum();
        }
        DAI_ASSERT( abs( logz_ - s ) < 1e-8 );
        return s;*/
  return logz_;
}


size_t JTree::findEfficientTree(const NodeSet& vs, RootedTree &Tree,
    size_t PreviousRoot) const {
  using dai::BigInt;
  using std::set;
  using std::vector;

  // find new root clique (the one with maximal statespace overlap with vs)
  BigInt maxval = 0;
  size_t maxalpha = 0;

  // reorder the tree edges such that maxalpha becomes the new root
  RootedTree newTree(GraphEL(RTree.begin(), RTree.end()), maxalpha);

  // identify subtree that contains all variables of vs which are not in the new root
  set<DEdge> subTree;
  // for each variable in vs
  for( NodeSet::const_iterator n = vs.begin(); n != vs.end(); n++ ) {
    for( size_t e = 0; e < newTree.size(); e++ ) {
      if( OR(newTree[e].second).nodes().contains( *n ) ) {
        size_t f = e;
        subTree.insert( newTree[f] );
        size_t pos = newTree[f].first;
        for( ; f > 0; f-- )
          if( newTree[f-1].second == pos ) {
            subTree.insert( newTree[f-1] );
            pos = newTree[f-1].first;
          }
      }
    }
  }
  if( PreviousRoot != (size_t)-1 && PreviousRoot != maxalpha) {
    // find first occurence of PreviousRoot in the tree, which is closest to the new root
    size_t e = 0;
    for( ; e != newTree.size(); e++ ) {
      if( newTree[e].second == PreviousRoot )
        break;
    }
    DAI_ASSERT( e != newTree.size() );

    // track-back path to root and add edges to subTree
    subTree.insert( newTree[e] );
    size_t pos = newTree[e].first;
    for( ; e > 0; e-- )
      if( newTree[e-1].second == pos ) {
        subTree.insert( newTree[e-1] );
        pos = newTree[e-1].first;
      }
  }

  // Resulting Tree is a reordered copy of newTree
  // First add edges in subTree to Tree
  Tree.clear();
  vector<DEdge> remTree;
  for (RootedTree::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
    if( subTree.count( *e ) )
      Tree.push_back( *e );
    else
      remTree.push_back( *e );
  size_t subTreeSize = Tree.size();
  // Then add remaining edges
  copy( remTree.begin(), remTree.end(), back_inserter( Tree ) );

  return subTreeSize;
}

Factor JTree::calcMarginal(const NodeSet& vs) {
  using std::vector;
  using std::map;

  vector<Factor>::const_iterator beta;
  for( beta = Qb.begin(); beta != Qb.end(); beta++ )
    if( beta->nodes() >> vs )
      break;
  if( beta != Qb.end() ) {
    if( props_.inference == Properties::InfType::SUMPROD )
      return( beta->marginal(vs) );
    else
      return( beta->maxMarginal(vs) );
  } else {
    vector<Factor>::const_iterator alpha;
    for( alpha = Qa.begin(); alpha != Qa.end(); alpha++ )
      if( alpha->nodes() >> vs )
        break;
    if( alpha != Qa.end() ) {
      if( props_.inference == Properties::InfType::SUMPROD )
        return( alpha->marginal(vs) );
      else
        return( alpha->maxMarginal(vs) );
    } else {
      // Find subtree to do efficient inference
      RootedTree T;
      size_t Tsize = findEfficientTree( vs, T );

      // Find remaining variables (which are not in the new root)
      NodeSet vsrem = vs / OR(T.front().first).nodes();
      Factor Pvs (vs);

      // Save Qa and Qb on the subtree
      map<size_t,Factor> Qa_old;
      map<size_t,Factor> Qb_old;
      vector<size_t> b(Tsize, 0);
      for( size_t i = Tsize; (i--) != 0; ) {
        size_t alpha1 = T[i].first;
        size_t alpha2 = T[i].second;
        size_t beta;
        for( beta = 0; beta < nrIRs(); beta++ )
          if (UEdge( RTree[beta].first, RTree[beta].second ) == UEdge( alpha1, alpha2 ) )
            break;
        DAI_ASSERT( beta != nrIRs() );
        b[i] = beta;

        if( !Qa_old.count( alpha1 ) )
          Qa_old[alpha1] = Qa[alpha1];
        if( !Qa_old.count( alpha2 ) )
          Qa_old[alpha2] = Qa[alpha2];
        if( !Qb_old.count( beta ) )
          Qb_old[beta] = Qb[beta];
      }

      return( Pvs.normalized() );
    }
  }
}


std::pair<size_t, dai::BigInt> boundTreewidth( const FactorGraph &fg, greedyVariableElimination::eliminationCostFunction fn, size_t maxStates ) {
  using dai::BigInt;

  // Create cluster graph from factor graph
  ClusterGraph _cg(fg, true);

  // Obtain elimination sequence
  std::vector<NodeSet> ElimVec = _cg.VarElim(greedyVariableElimination(fn), maxStates ).eraseNonMaximal().clusters();

  // Calculate treewidth
  size_t treewidth = 0;
  BigInt nrstates = 0.0;
  for( size_t i = 0; i < ElimVec.size(); i++ ) {
    if( ElimVec[i].size() > treewidth )
      treewidth = ElimVec[i].size();
    // BigInt s = ElimVec[i].nrStates();
    // if( s > nrstates )
    //   nrstates = s;
  }

  return std::make_pair(treewidth, nrstates);
}

std::vector<Node*> JTree::findMaximum() const {
  using std::vector;
  using std::stack;

  vector<Node*> maximum(nrNodes());
  vector<bool> visitedVars(nrNodes(), false);
  vector<bool> visitedORs(nrORs(), false);
  stack<size_t> scheduledORs;
  scheduledORs.push(0);
  while( !scheduledORs.empty() ) {
    size_t alpha = scheduledORs.top();
    scheduledORs.pop();
    if( visitedORs[alpha] )
      continue;
    visitedORs[alpha] = true;

    // The allowed configuration is restrained according to the variables assigned so far:
    // pick the argmax amongst the allowed states
    Real maxProb = -std::numeric_limits<Real>::max();
    size_t maxcount = 0;
  }
  return maximum;
}

size_t JTree::FindRoot(const vector<NodeSet> &cl) const {
  size_t root = cl.size();
  if (props_.root.size() > 0) {
    size_t lightest = -1UL; 
    root = 0;
    for (size_t i = 0; i < cl.size(); i++) {
      if (cl[i] >> props_.root) {
        size_t s = cl[i].size();
        if (s < lightest) {
          lightest = s;
          root = i;
        }
      }
    }
  }
  return root;
}
} // end of namespace dai
