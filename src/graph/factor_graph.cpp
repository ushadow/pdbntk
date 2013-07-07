#include "factor_graph.h"
#include "glog/logging.h"

#include <iostream>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <dai/util.h>
#include <dai/exceptions.h>
#include <boost/lexical_cast.hpp>

namespace pdbntk {

FactorGraph::FactorGraph(const std::vector<Factor> &P) : _G(), _backup() {
  using std::set;
  // add factors, obtain nodes. 
  set<Node*, NodeComparator> varset;
  _factors.reserve(P.size());
  size_t nrEdges = 0;
  for (std::vector<Factor>::const_iterator p2 = P.begin(); p2 != P.end(); p2++ ) {
    _factors.push_back( *p2 );
    copy(p2->nodes().begin(), p2->nodes().end(), inserter(varset, varset.begin()));
    nrEdges += p2->nodes().size();
  }

  // add nodes
  nodes_.reserve( varset.size() );
  for (set<Node*>::const_iterator p1 = varset.begin(); p1 != varset.end(); p1++ ) {
    nodes_.push_back(*p1);
  }

  // create graph structure
  constructGraph( nrEdges );
}

void FactorGraph::constructGraph( size_t nrEdges ) {
  using dai::Edge;

  // create a mapping for indices
  dai::hash_map<size_t, size_t> hashmap;

  for( size_t i = 0; i < nodes().size(); i++ )
    hashmap[node(i)->index()] = i;

  // create edge list
  std::vector<Edge> edges;
  edges.reserve( nrEdges );
  for( size_t i2 = 0; i2 < nrFactors(); i2++ ) {
    const NodeSet& ns = factor(i2).nodes();
    for( NodeSet::const_iterator q = ns.begin(); q != ns.end(); q++ )
      edges.push_back(Edge(hashmap[(*q)->index()], i2) );
  }

  // create bipartite graph
  _G.construct( nrNodes(), nrFactors(), edges.begin(), edges.end() );
}


/// Writes a FactorGraph to an output stream
std::ostream& operator<< ( std::ostream &os, const FactorGraph &fg ) {
  using std::endl;

  os << fg.nrFactors() << endl;

  for( size_t I = 0; I < fg.nrFactors(); I++ ) {
    os << endl;
    os << fg.factor(I).nodes().size() << endl;
    for( NodeSet::const_iterator i = fg.factor(I).nodes().begin(); i != fg.factor(I).nodes().end(); i++ )
      os << (*i)->index() << " ";
    os << endl;
    for( NodeSet::const_iterator i = fg.factor(I).nodes().begin(); i != fg.factor(I).nodes().end(); i++ )
    os << endl;
    size_t nr_nonzeros = 0;
    os << nr_nonzeros << endl;
  }

  return(os);
}


/// Reads a FactorGraph from an input stream
std::istream& operator>> ( std::istream& is, FactorGraph &fg ) {
  using std::vector;
  using std::string;

  vector<Factor> facs;
  size_t nr_Factors;
  string line;

  while( (is.peek()) == '#' )
    getline(is,line);
  is >> nr_Factors;
  if( is.fail() )
    DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of factors");
  DLOG(INFO) << "Reading " << nr_Factors << " factors...";

  getline (is,line);
  if( is.fail() || line.size() > 0 )
    DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Expecting empty line");

  std::map<long,size_t> vardims;
  for( size_t I = 0; I < nr_Factors; I++ ) {
    DLOG(INFO) << "Reading factor " << I << "...";
    size_t nr_members;
    while( (is.peek()) == '#' )
      getline(is,line);
    is >> nr_members;
    DLOG(INFO) << "  nr_members: " << nr_members;

    vector<long> labels;
    for( size_t mi = 0; mi < nr_members; mi++ ) {
      long mi_label;
      while( (is.peek()) == '#' )
        getline(is,line);
      is >> mi_label;
      labels.push_back(mi_label);
    }

    vector<size_t> dims;
    for( size_t mi = 0; mi < nr_members; mi++ ) {
      size_t mi_dim;
      while( (is.peek()) == '#' )
        getline(is,line);
      is >> mi_dim;
      dims.push_back(mi_dim);
    }

    // read values
    size_t nr_nonzeros;
    while( (is.peek()) == '#' )
      getline(is,line);
    is >> nr_nonzeros;
    DLOG(INFO) << "  nonzeroes: " << nr_nonzeros;
    for( size_t k = 0; k < nr_nonzeros; k++ ) {
      size_t li;
      Real val;
      while( (is.peek()) == '#' )
        getline(is,line);
      is >> li;
      while( (is.peek()) == '#' )
        getline(is,line);
      is >> val;
    }
  }

  DLOG(INFO) << "factors:" << facs;

  fg = FactorGraph(facs);

  return is;
}


NodeSet FactorGraph::Delta( size_t i ) const {
  // calculate Markov Blanket
  NodeSet Del;
  bforeach(const dai::Neighbor &I, nbV(i)) // for all neighboring factors I of i
    bforeach(const dai::Neighbor &j, nbF(I)) // for all neighboring variables j of I
    Del |= node(j);

  return Del;
}


NodeSet FactorGraph::Delta( const NodeSet &ns ) const {
  NodeSet result;
  for (NodeSet::const_iterator n = ns.begin(); n != ns.end(); n++)
    result |= Delta(findNode(*n));
  return result;
}


void FactorGraph::makeCavity( size_t i, bool backup ) {
  // fills all Factors that include var(i) with ones
  std::map<size_t,Factor> newFacs;
  bforeach( const dai::Neighbor &I, nbV(i) ) // for all neighboring factors I of i
    newFacs[I] = Factor(factor(I).nodes());
  setFactors( newFacs, backup );
}


void FactorGraph::ReadFromFile( const char *filename ) {
  std::ifstream infile;
  infile.open( filename );
  if( infile.is_open() ) {
    infile >> *this;
    infile.close();
  } else
    DAI_THROWE(CANNOT_READ_FILE,"Cannot read from file " + std::string(filename));
}


void FactorGraph::WriteToFile( const char *filename, size_t precision ) const {
  std::ofstream outfile;
  outfile.open( filename );
  if( outfile.is_open() ) {
    outfile.precision( precision );
    outfile << *this;
    outfile.close();
  } else
    DAI_THROWE(CANNOT_WRITE_FILE,"Cannot write to file " + std::string(filename));
}


void FactorGraph::printDot( std::ostream &os ) const {
  using std::endl;
  os << "graph FactorGraph {" << endl;
  os << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
  for( size_t i = 0; i < nrNodes(); i++ )
    os << "\tv" << node(i)->index() << ";" << endl;
  os << "node[shape=box,width=0.3,height=0.3,fixedsize=true];" << endl;
  for( size_t I = 0; I < nrFactors(); I++ )
    os << "\tf" << I << ";" << endl;
  for( size_t i = 0; i < nrNodes(); i++ )
    bforeach( const dai::Neighbor &I, nbV(i) )  // for all neighboring factors I of i
      os << "\tv" << node(i)->index() << " -- f" << I << ";" << endl;
  os << "}" << endl;
}


dai::GraphAL FactorGraph::MarkovGraph() const {
  dai::GraphAL G( nrNodes() );
  for( size_t i = 0; i < nrNodes(); i++ )
    bforeach( const dai::Neighbor &I, nbV(i) )
      bforeach( const dai::Neighbor &j, nbF(I) )
      if( i < j )
        G.addEdge( i, j, true );
  return G;
}


bool FactorGraph::isMaximal( size_t I ) const {
  const NodeSet& I_vars = factor(I).nodes();
  size_t I_size = I_vars.size();

  if( I_size == 0 ) {
    for( size_t J = 0; J < nrFactors(); J++ ) 
      if( J != I )
        if( factor(J).nodes().size() > 0 )
          return false;
    return true;
  } else {
    bforeach( const dai::Neighbor& i, nbF(I) ) {
      bforeach( const dai::Neighbor& J, nbV(i) ) {
        if( J != I )
          if( (factor(J).nodes() >> I_vars) && (factor(J).nodes().size() != I_size) )
            return false;
      }
    }
    return true;
  }
}


size_t FactorGraph::maximalFactor( size_t I ) const {
  const NodeSet& I_vars = factor(I).nodes();
  size_t I_size = I_vars.size();

  if( I_size == 0 ) {
    for( size_t J = 0; J < nrFactors(); J++ )
      if( J != I )
        if( factor(J).nodes().size() > 0 )
          return maximalFactor( J );
    return I;
  } else {
    bforeach( const dai::Neighbor& i, nbF(I) ) {
      bforeach( const dai::Neighbor& J, nbV(i) ) {
        if( J != I )
          if( (factor(J).nodes() >> I_vars) && (factor(J).nodes().size() != I_size) )
            return maximalFactor( J );
      }
    }
    return I;
  }
}


std::vector<NodeSet> FactorGraph::maximalFactorDomains() const {
  std::vector<NodeSet> result;

  for( size_t I = 0; I < nrFactors(); I++ )
    if( isMaximal( I ) )
      result.push_back( factor(I).nodes() );

  if( result.size() == 0 )
    result.push_back( NodeSet() );
  return result;
}


Real FactorGraph::logScore(const std::vector<size_t>& statevec) const {
  // Construct a State object that represents statevec
  // This decouples the representation of the joint state in statevec from the factor graph
  std::map<const Node*, size_t> statemap;
  for( size_t i = 0; i < statevec.size(); i++ )
    statemap[node(i)] = statevec[i];

  // Evaluate the log probability of the joint configuration in statevec
  // by summing the log factor entries of the factors that correspond to this joint configuration
  Real lS = 0.0;
  return lS;
}

void FactorGraph::clamp(size_t i, const std::vector<Real> &x, bool backup) {
  Factor mask(node(i));

  std::map<size_t, Factor> newFacs;
  bforeach(const dai::Neighbor &I, nbV(i)) {
    DLOG(INFO) << factor(I);
    newFacs[I] = factor(I) * mask;
  }
  setFactors(newFacs, backup);

  return;
}

void FactorGraph::clampVar( size_t i, const std::vector<size_t> &is, bool backup ) {
  Node *n = node(i);
  Factor mask_n(n);

  std::map<size_t, Factor> newFacs;
  bforeach(const dai::Neighbor &I, nbV(i))
    newFacs[I] = factor(I) * mask_n;
  setFactors( newFacs, backup );
}


void FactorGraph::clampFactor( size_t I, const std::vector<size_t> &is, bool backup ) {
  Factor newF(factor(I).nodes());

  setFactor(I, newF, backup);
}


void FactorGraph::backupFactor( size_t I ) {
  std::map<size_t,Factor>::iterator it = _backup.find( I );
  if( it != _backup.end() )
    DAI_THROW(MULTIPLE_UNDO);
  _backup[I] = factor(I);
}


void FactorGraph::restoreFactor( size_t I ) {
  std::map<size_t,Factor>::iterator it = _backup.find( I );
  if( it != _backup.end() ) {
    setFactor(I, it->second);
    _backup.erase(it);
  } else
    DAI_THROW(OBJECT_NOT_FOUND);
}

void FactorGraph::backupFactors( const NodeSet &ns ) {
  for( size_t I = 0; I < nrFactors(); I++ )
    if( factor(I).nodes().intersects( ns ) )
      backupFactor( I );
}


void FactorGraph::restoreFactors( const NodeSet &ns ) {
  using std::map;
  map<size_t,Factor> facs;
  for( map<size_t,Factor>::iterator uI = _backup.begin(); uI != _backup.end(); ) {
    if( factor(uI->first).nodes().intersects( ns ) ) {
      facs.insert( *uI );
      _backup.erase(uI++);
    } else
      uI++;
  }
  setFactors( facs );
}


void FactorGraph::restoreFactors() {
  setFactors( _backup );
  _backup.clear();
}


void FactorGraph::backupFactors( const std::set<size_t> & facs ) {
  for( std::set<size_t>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ )
    backupFactor( *fac );
}


bool FactorGraph::isPairwise() const {
  bool pairwise = true;
  for( size_t I = 0; I < nrFactors() && pairwise; I++ )
    if( factor(I).nodes().size() > 2 )
      pairwise = false;
  return pairwise;
}


bool FactorGraph::isBinary() const {
  bool binary = true;
  return binary;
}


FactorGraph FactorGraph::clamped( size_t i, size_t state ) const {
  Node* v = node(i);
  Real zeroth_order = (Real)1;
  std::vector<Factor> clamped_facs;
  for( size_t I = 0; I < nrFactors(); I++ ) {
    NodeSet v_I = factor(I).nodes();
    Factor new_factor;

    if( new_factor.nodes().size() != 0 ) {
      size_t J = 0;
      // if it can be merged with a previous one, do that
      for( J = 0; J < clamped_facs.size(); J++ )
        if( clamped_facs[J].nodes() == new_factor.nodes() ) {
          clamped_facs[J] *= new_factor;
          break;
        }
      // otherwise, push it back
      if( J == clamped_facs.size() || clamped_facs.size() == 0 )
        clamped_facs.push_back( new_factor );
    } 
  }
  *(clamped_facs.begin()) *= zeroth_order;
  return FactorGraph( clamped_facs );
}


FactorGraph FactorGraph::maximalFactors() const {
  using std::vector;

  vector<size_t> maxfac( nrFactors() );
  std::map<size_t,size_t> newindex;
  size_t nrmax = 0;
  for( size_t I = 0; I < nrFactors(); I++ ) {
    maxfac[I] = I;
    NodeSet maxfacvars = factor(maxfac[I]).nodes();
    for( size_t J = 0; J < nrFactors(); J++ ) {
      NodeSet Jvars = factor(J).nodes();
      if( Jvars >> maxfacvars && (Jvars != maxfacvars) ) {
        maxfac[I] = J;
        maxfacvars = factor(maxfac[I]).nodes();
      }
    }
    if( maxfac[I] == I )
      newindex[I] = nrmax++;
  }

  vector<Factor> facs( nrmax );
  for( size_t I = 0; I < nrFactors(); I++ )
    facs[newindex[maxfac[I]]] *= factor(I);

  return FactorGraph( facs.begin(), facs.end(), nodes().begin(), nodes().end(), facs.size(), nrNodes() );
}


} // end of namespace dai
