#ifndef PDBNTK_CLUSTER_GRAPH_H_ 
#define PDBNTK_CLUSTER_GRAPH_H_ 

#include "dai/bipgraph.h"
#include "node.h"
#include "utils.h"
#include "factor_graph.h"

#include <set>
#include <vector>

namespace pdbntk {

/// A ClusterGraph is a hypergraph with variables as nodes, and "clusters" (sets of variables) as hyperedges.
/** It is implemented as a bipartite graph with variable (Var) nodes and cluster (NodeSet) nodes.
 *  One may think of a ClusterGraph as a FactorGraph without the actual factor values.
 *  \todo Remove the nodes_ and _clusters variables and use only the graph and a contextual factor graph.
 */
class ClusterGraph {
private:
  /// Stores the neighborhood structure
  dai::BipartiteGraph       _G;

  /// Stores the variables corresponding to the nodes
  std::vector<Node*>     nodes_;

  /// Stores the clusters corresponding to the hyperedges
  std::vector<NodeSet>  _clusters;

public:
  /// \name Constructors and destructors
  //@{
  /// Default constructor
  ClusterGraph() : _G(), nodes_(), _clusters() {}

  /// Construct from vector of NodeSet 's
  ClusterGraph( const std::vector<NodeSet>& cls );

  /// Construct from a factor graph
  /** Creates cluster graph which has factors in \a fg as clusters if \a onlyMaximal == \c false,
   *  and only the maximal factors in \a fg if \a onlyMaximal == \c true.
   */
  ClusterGraph(const FactorGraph& fg, bool onlyMaximal);
  //@}

  /// \name Queries
  //@{
  /// Returns a constant reference to the graph structure
  const dai::BipartiteGraph& bipGraph() const { return _G; }

  /// Returns number of variables
  size_t nrVars() const { return nodes_.size(); }

  /// Returns a constant reference to the variables
  const std::vector<Node*>& nodes() const { return nodes_; }

  /// Returns a constant reference to the \a i'th variable
  const Node* node(size_t i) const {
    DAI_DEBASSERT( i < nrVars() );
    return nodes_[i]; 
  }

  /// Returns number of clusters
  size_t nrClusters() const { return _clusters.size(); }

  /// Returns a constant reference to the clusters
  const std::vector<NodeSet>& clusters() const { return _clusters; }

  /// Returns a constant reference to the \a I'th cluster
  const NodeSet& cluster( size_t I ) const {
    DAI_DEBASSERT( I < nrClusters() );
    return _clusters[I]; 
  }

  /// Returns the index of variable \a n
  size_t findNode(const Node *n) const {
    return find(nodes_.begin(), nodes_.end(), n) - nodes_.begin();
  }

  /// Returns the index of a cluster \a cl
  size_t findCluster( const NodeSet& cl ) const {
    return find( _clusters.begin(), _clusters.end(), cl ) - _clusters.begin();
  }

  /// Returns union of clusters that contain the \a i 'th variable
  NodeSet Delta( size_t i ) const {
    NodeSet result;
    bforeach( const dai::Neighbor& I, _G.nb1(i) )
      result |= _clusters[I];
    return result;
  }

  /// Returns union of clusters that contain the \a i 'th (except this variable itself)
  NodeSet delta( size_t i ) const {
    return Delta( i ) / nodes_[i];
  }

  /// Returns \c true if variables with indices \a i1 and \a i2 are adjacent, i.e., both contained in the same cluster
  bool adj( size_t i1, size_t i2 ) const {
    if( i1 == i2 )
      return false;
    bool result = false;
    bforeach( const dai::Neighbor& I, _G.nb1(i1) )
      if( find( _G.nb2(I).begin(), _G.nb2(I).end(), i2 ) != _G.nb2(I).end() ) {
        result = true;
        break;
      }
    return result;
  }

  /// Returns \c true if cluster \a I is not contained in a larger cluster
  bool isMaximal( size_t I ) const {
    DAI_DEBASSERT( I < _G.nrNodes2() );
    const NodeSet & clI = _clusters[I];
    bool maximal = true;
    // The following may not be optimal, since it may repeatedly test the same cluster *J
    bforeach( const dai::Neighbor& i, _G.nb2(I) ) {
      bforeach( const dai::Neighbor& J, _G.nb1(i) )
        if( (J != I) && (clI << _clusters[J]) ) {
          maximal = false;
          break;
        }
      if( !maximal )
        break;
    }
    return maximal;
  }
  //@}

  /// \name Operations
  //@{
  /// Inserts a cluster (if it does not already exist) and creates new variables, if necessary
  /** \note This function could be better optimized if the index of one variable in \a cl would be known.
   *        If one could assume nodes_ to be ordered, a binary search could be used instead of a linear one.
   */
  size_t insert(const NodeSet& cl) {
    size_t index = findCluster(cl);  // OPTIMIZE ME
    if( index == _clusters.size() ) {
      _clusters.push_back( cl );
      // add variables (if necessary) and calculate neighborhood of new cluster
      std::vector<size_t> nbs;
      for( NodeSet::const_iterator n = cl.begin(); n != cl.end(); n++ ) {
        size_t iter = findNode( *n );  // OPTIMIZE ME
        nbs.push_back( iter );
        if( iter == nodes_.size() ) {
          _G.addNode1();
          nodes_.push_back( *n );
        }
      }
      _G.addNode2( nbs.begin(), nbs.end(), nbs.size() );
    }
    return index;
  }

  /// Erases all clusters that are not maximal
  ClusterGraph& eraseNonMaximal() {
    for( size_t I = 0; I < _G.nrNodes2(); ) {
      if( !isMaximal(I) ) {
        _clusters.erase( _clusters.begin() + I );
        _G.eraseNode2(I);
      } else
        I++;
    }
    return *this;
  }

  /// Erases all clusters that contain the \a i 'th variable
  ClusterGraph& eraseSubsuming( size_t i ) {
    DAI_ASSERT( i < nrVars() );
    while( _G.nb1(i).size() ) {
      _clusters.erase( _clusters.begin() + _G.nb1(i)[0] );
      _G.eraseNode2( _G.nb1(i)[0] );
    }
    return *this;
  }

  /// Eliminates variable with index \a i, without deleting the variable itself
  /** \note This function can be better optimized
  */
  NodeSet elimNode(size_t i) {
    DAI_ASSERT( i < nrVars() );
    NodeSet Di = Delta( i );
    insert( Di / const_cast<Node*>(node(i)));
    eraseSubsuming( i );
    eraseNonMaximal();
    return Di;
  }
  //@}

  /// \name Input/Ouput
  //@{
  /// Writes a ClusterGraph to an output stream
  friend std::ostream& operator << (std::ostream& os, const ClusterGraph& cl) {
    os << cl.clusters();
    return os;
  }
  //@}

  /// \name Variable elimination
  //@{
  /// Performs Variable Elimination, keeping track of the interactions that are created along the way.
  /** \tparam EliminationChoice should support "size_t operator()( const ClusterGraph &cl, const std::set<size_t> &remainingVars )"
   *  \param f function object which returns the next variable index to eliminate; for example, a dai::greedyVariableElimination object.
   *  \param maxStates maximum total number of states of all clusters in the output cluster graph (0 means no limit).
   *  \throws OUT_OF_MEMORY if total number of states becomes larger than maxStates
   *  \return A set of elimination "cliques".
   */
  template<class EliminationChoice>
  ClusterGraph VarElim(EliminationChoice f, size_t maxStates = 0) const {
    // Make a copy
    ClusterGraph cl(*this);
    cl.eraseNonMaximal();

    // Creates an empty cluster graph.
    ClusterGraph result;

    // Construct set of variable indices
    std::set<size_t> varindices;
    for (size_t i = 0; i < nodes_.size(); ++i)
      varindices.insert( i );

    // Do variable elimination
    dai::BigInt totalStates = 0;
    while (!varindices.empty()) {
      size_t i = f(cl, varindices);
      NodeSet Di = cl.elimNode(i);
      result.insert( Di );
      varindices.erase( i );
    }

    return result;
  }
  //@}
};


/// Helper object for dai::ClusterGraph::VarElim()
/** Chooses the next variable to eliminate by picking them sequentially from a given vector of variables.
*/
class sequentialVariableElimination {
private:
  /// The variable elimination sequence
  std::vector<Node*> seq;
  /// Counter
  size_t i;

public:
  /// Construct from vector of variables
  sequentialVariableElimination( const std::vector<Node*> s ) : seq(s), i(0) {}

  /// Returns next variable in sequence
  size_t operator()( const ClusterGraph &cl, const std::set<size_t> &/*remainingVars*/ );
};


/// Helper object for dai::ClusterGraph::VarElim()
/** Chooses the next variable to eliminate greedily by taking the one that minimizes
*  a given heuristic cost function.
*/
class greedyVariableElimination {
public:
  /// Type of cost functions to be used for greedy variable elimination
  typedef size_t (*eliminationCostFunction)(const ClusterGraph &, size_t);

private:
  /// Pointer to the cost function used
  eliminationCostFunction heuristic;

public:
  /// Construct from cost function
  /** \note Examples of cost functions are eliminationCost_MinFill() and eliminationCost_WeightedMinFill().
  */
  greedyVariableElimination( eliminationCostFunction h ) : heuristic(h) {}

  /// Returns the best variable from \a remainingVars to eliminate in the cluster graph \a cl by greedily minimizing the cost function.
  /** This function calculates the cost for eliminating each variable in \a remaingVars and returns the variable which has lowest cost.
  */
  size_t operator()( const ClusterGraph &cl, const std::set<size_t>& remainingVars );
};


/// Calculates cost of eliminating the \a i 'th variable from cluster graph \a cl according to the "Mindai::Neighbors" criterion.
/** The cost is measured as "number of neigboring nodes in the current adjacency graph",
 *  where the adjacency graph has the variables as its nodes and connects
 *  nodes \a i1 and \a i2 iff \a i1 and \a i2 occur together in some common cluster.
 */
size_t eliminationCost_MinNeighbors( const ClusterGraph& cl, size_t i );


/// Calculates cost of eliminating the \a i 'th variable from cluster graph \a cl according to the "MinWeight" criterion.
/** The cost is measured as "product of weights of neighboring nodes in the current adjacency graph",
 *  where the adjacency graph has the variables as its nodes and connects
 *  nodes \a i1 and \a i2 iff \a i1 and \a i2 occur together in some common cluster.
 *  The weight of a node is the number of states of the corresponding variable.
 */
size_t eliminationCost_MinWeight( const ClusterGraph& cl, size_t i );


/// Calculates cost of eliminating the \a i 'th variable from cluster graph \a cl according to the "MinFill" criterion.
/** The cost is measured as "number of added edges in the adjacency graph",
 *  where the adjacency graph has the variables as its nodes and connects
 *  nodes \a i1 and \a i2 iff \a i1 and \a i2 occur together in some common cluster.
 */
size_t eliminationCost_MinFill( const ClusterGraph& cl, size_t i );


/// Calculates cost of eliminating the \a i 'th variable from cluster graph \a cl according to the "WeightedMinFill" criterion.
/** The cost is measured as "total weight of added edges in the adjacency graph",
 *  where the adjacency graph has the variables as its nodes and connects
 *  nodes \a i1 and \a i2 iff \a i1 and \a i2 occur together in some common cluster.
 *  The weight of an edge is the product of the number of states of the variables corresponding with its nodes.
 */
size_t eliminationCost_WeightedMinFill( const ClusterGraph& cl, size_t i );
} 

#endif // PDBNTK_CLUSTER_GRAPH_H_
