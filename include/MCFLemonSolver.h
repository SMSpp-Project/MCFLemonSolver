/*--------------------------------------------------------------------------*/
/*------------------------ File MCFLemonSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MCFLemonSolver class, implementing the Solver
 * interface, in particular in its CDASolver version, for Min-Cost Flow
 * problems as set by MCFBlock.
 *
 * This is based on interfacing algorithms implemented in the LEMON
 * (Library for Efficient Modeling and Optimization in Networks) project,
 * as currently found at
 *
 *     https://lemon.cs.elte.hu/trac/lemon
 *
 * \author Daniele Caliandro \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Daniele Caliandro, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFLemonSolver
 #define __MCFLemonSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "MCFBlock.h"

#include "lemon/capacity_scaling.h"

#include "lemon/cost_scaling.h"

#include "lemon/cycle_canceling.h"

#include "lemon/network_simplex.h"

#include "lemon/k_best_extreme_flows.h"

#include <lemon/quad_heap.h>

#include <lemon/list_graph.h>

#include <lemon/smart_graph.h>

#include <type_traits>

#include <array>

#include <chrono>

#include <memory>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/
/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace lemon;

 /*--------------------------------------------------------------------------*/
 /*------------------------------- CLASSES ----------------------------------*/
 /*--------------------------------------------------------------------------*/
 /** @defgroup MCFLemonSolver_CLASSES Classes in MCFLemonSolver.h
  *  @{ */

 /*--------------------------------------------------------------------------*/
 /*------------------------ CLASS MCFListDigraph ----------------------------*/
 /*--------------------------------------------------------------------------*/
 /*--------------------------- GENERAL NOTES --------------------------------*/
 /*--------------------------------------------------------------------------*/
 /*--------------------------------------------------------------------------*/
 /** MCFListDigraph extends ListDigraph customizing functions that provides
  * access to the structure of the graph used by f_algo. This class comes
  * from two needs:
  *
  * - In MCFBlock, the when a new arc is created its "name" is that of the
  *   arc with smallest name that has been previously deleted and not
  *   re-added yet (the name of the last arc if there are no such arcs).
  *   This is not so in the original ListDigraph, so the implementation
  *   of addArc() has to be changed accordingly.
  *
  * - The need of implementing the openArc() and closeArc() methods, that
  *   realise the arc closing / opening operations of MCFBlock. A "closed" arc
  *   is basically a deleted arc in the original ListDigraph, except that its
  *   "name" is not available when new arcs are constructed, so that it can
  *   be re-opened keeping its original name, capacity and cost. The
  *   implementation is based on making the startID of the arc negative but
  *   encoding the original value ( - ( startID + 1 ) ), so that when opening
  *   it back it is possible to restore the original value of the startID. An
  *   implementation of eraseClosed() is also provided since close() modifies
  *   the structure of the arc object without "marking" it deleted, and since
  *   an arc can be removed while it is closed, we need to "mark" it with the
  *   standard LEMON operations to recognise deleted arcs (that is, set the
  *   prev_in field to -2 and add it to the list of deleted arcs).
  *
  * The MCFArc class is also provided that derives from ListDigraph::Arc
  * since the constructor of that class must be used in the addArc()
  * method, because it sets the value of the arcID in the arc object
  * required for iterate over the graph structure of LEMON classes.
  * This is originally done by declaring ListDigraph as friend of Arc
  * (questionable decision ...), but friendship is not inherited, so it
  * has to be done again. */

class MCFListDigraph : public ListDigraph
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/

  using ListDigraphBase::Arc;
  using ListDigraphBase::first_free_arc;
  using ListDigraphBase::Node;

  /// basically ListDigraph::Arc, just made friend of MCFListDigraph
  class MCFArc : public ListDigraph::Arc
  {
   friend class MCFListDigraph;

   protected:
    MCFArc( int pid ) { id = pid; }
   };

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

 /// constructor, has nothing special to do
 MCFListDigraph() {}

 /// destructor = default
 ~MCFListDigraph() = default;

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// return true if the arc is closed

 bool isClosed( int n ) { return( arc_store()[ n ].source < 0 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// add a new arc with given endpoints

 MCFArc addArc( int u , int v )
 {
  auto & arcs = arc_store();
  auto & nodes = node_store();

  int n;

  if( first_free_arc == -1 ) {  // no previously deleted arcs
   n = arcs.size();            // just add at the end
   arcs.push_back( ArcT() );
   }
  else {                        // there are previously deleted arcs
   // look for the min-index arc in the deleted arc list
   int minidx = first_free_arc;
   int idx = first_free_arc;
   int prevmin = -1;
   for( auto nextidx = arcs[ idx ].next_in ; nextidx != -1 ;
	       idx = nextidx , nextidx = arcs[ idx ].next_in )
    if( nextidx < minidx ) {
     minidx = nextidx;
     prevmin = idx;
     }

   if( prevmin == -1 )
    first_free_arc = arcs[ first_free_arc ].next_in;
   else
    arcs[ prevmin ].next_in = arcs[ minidx ].next_in;

   arcs[ minidx ].source = u;
   arcs[ minidx ].target = v;

   arcs[ minidx ].next_out = nodes[ u ].first_out;
   if( nodes[ u ].first_out != -1 )
    arcs[ nodes[ u ].first_out ].prev_out = minidx;

   arcs[ minidx ].next_in = nodes[ v ].first_in;
   if( nodes[ v ].first_in != -1 )
    arcs[ nodes[ v ].first_in ].prev_in = minidx;

   arcs[ minidx ].prev_in = arcs[ minidx ].prev_out = -1;

   nodes[ u ].first_out = nodes[ v ].first_in = minidx;

   return( notify_added( minidx ) );
   }

  // otherwise add to the last position
  arcs[ n ].source = u;
  arcs[ n ].target = v;

  arcs[ n ].next_out = nodes[ u ].first_out;
  if( nodes[ u ].first_out != -1 )
   arcs[ nodes[ u ].first_out ].prev_out = n;

  arcs[ n ].next_in = nodes[ v ].first_in;
  if( nodes[ v ].first_in != -1 )
   arcs[ nodes[ v ].first_in ].prev_in = n;

  arcs[ n ].prev_in = arcs[ n ].prev_out = -1;

  nodes[ u ].first_out = nodes[ v ].first_in = n;

  return( notify_added( n ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// tell the ArcMap of the graph that arc n exists
 /** ListDigraph::addArc() notifies the observers of the arcs, i.e., every
  * ArcMap built on the graph, which is how a map grows with the graph;
  * addArc() above replaces it, so it has to do the same, or the next set()
  * on a map writes past the end of its storage. */

 MCFArc notify_added( int n )
 {
  MCFArc arc( n );
  notifier( Arc() ).add( arc );
  return( arc );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// close the given arc

 void closeArc( int n )
 {
  auto & arcs = arc_store();
  auto & nodes = node_store();

  // like in delete, remove from FS and BS lists
  if( arcs[ n ].next_in != -1 )
   arcs[ arcs[ n ].next_in ].prev_in = arcs[ n ].prev_in;

  if( arcs[ n ].prev_in != -1 )
   arcs[ arcs[ n ].prev_in ].next_in = arcs[ n ].next_in;
  else
   nodes[ arcs[ n ].target ].first_in = arcs[ n ].next_in;

  if( arcs[ n ].next_out != -1 )
   arcs[ arcs[ n ].next_out ].prev_out = arcs[ n ].prev_out;

  if( arcs[ n ].prev_out != -1 )
   arcs[ arcs[ n ].prev_out ].next_out = arcs[ n ].next_out;
  else
   nodes[ arcs[ n ].source ].first_out = arcs[ n ].next_out;

  // now just make the source negative
  arcs[ n ].source = -( arcs[ n ].source + 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// re-open the given arc

 void openArc( int n )
 {
  auto & arcs = arc_store();
  auto & nodes = node_store();

  // restore the original source
  arcs[ n ].source = -( arcs[ n ].source + 1 );
  auto uid = arcs[ n ].source;
  auto vid = arcs[ n ].target;

  // put it back into the FS and BS lists
  arcs[ n ].next_out = nodes[ uid ].first_out;
  if( nodes[ uid ].first_out != -1 )
   arcs[ nodes[ uid ].first_out ].prev_out = n;

  arcs[ n ].next_in = nodes[ vid ].first_in;
  if( nodes[ vid ].first_in != -1 )
   arcs[ nodes[ vid ].first_in ].prev_in = n;

  arcs[ n ].prev_in = arcs[ n ].prev_out = -1;

  nodes[ uid ].first_out = nodes[ vid ].first_in = n;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// erase a closed arc

 void eraseClosed( int n )
 {
  auto & arcs = arc_store();

  // just do the little that closeArc() did not
  arcs[ n ].next_in = first_free_arc;
  first_free_arc = n;
  arcs[ n ].prev_in = -2;
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
 /// the storage of the arcs of ListDigraphBase
 /** The vector is arcs in LEMON 1.3.1 and _arcs in the later versions (and
  * in lemon-ng), where arcs() is the range of a loop on the arcs: the one
  * that the base has is taken. */

 template< class G = MCFListDigraph >
 std::vector< ArcT > & arc_store( void ) {
  G & g = *this;
  if constexpr( requires { g._arcs; } )
   return( g._arcs );
  else
   return( g.arcs );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the storage of the nodes of ListDigraphBase [see arc_store()]

 template< class G = MCFListDigraph >
 std::vector< NodeT > & node_store( void ) {
  G & g = *this;
  if constexpr( requires { g._nodes; } )
   return( g._nodes );
  else
   return( g.nodes );
  }

 };  // end( class( MCFListDigraph ) )

 /*--------------------------------------------------------------------------*/
 /*-------------------------- TEMPLATE TYPES --------------------------------*/
 /*--------------------------------------------------------------------------*/
 /** Algorithms in the LEMON projects are template over at least three types:
  *
  * - GR, which is the type of graph. The possibilities are:
  *
  *   = SmartDigraph is a simple and fast directed graph implementation. It
  *     is quite memory efficient but at the price that it does not support
  *     node and arc deletion.
  *
  *   = MCFListDigraph, an ad-hoc simple extension of ListDigraph, which is
  *     a standard versatile and fast directed graph implementation based on
  *     linked lists that are stored in std::vector structures. This class
  *     provides only linear time counting for nodes and arcs. The original
  *     ListDigraph supports node and arc deletion, and MCFListDigraph extends
  *     this to arc closure and re-opening.
  *
  * - V, which is the type of flows / deficits; typically, double can be used
  *   for maximum compatibility, but int (or even smaller) would yield better
  *   performances;
  *
  * - C, which is the type of arc costs; typically, double can be used for
  *   maximum compatibility, but int (or even smaller) would yield better
  *   performances.
  *
  * Furthermore, scaling-type algorithms may behave in different ways
  * according to which combination of V and C is used, and there are different
  * "traits" for this which are another scaling parameter. However, we prefer
  * that the MCFLemonSolver class is always template over the first three
  * parameters only, which is why we define SMSppCapacityScaling and
  * SMSppCostScaling as template over < GR , V , C >, each with the trait
  * it is given below.
  *
  * Thus, the concept LEMONGraph is defined that only allows all possible
  * types of LEMON graphs, i.e., SmartDigraph and MCFListDigraph. */

 /// concept for "one of the LEMON graphs"
 template< typename Type >
 concept LEMONGraph = std::is_same< Type , SmartDigraph >::value ||
                      std::is_same< Type , MCFListDigraph >::value;

 /// the trait of CapacityScaling with a 4-ary heap
 /** The default trait of CapacityScaling with the heap of the shortest path
  * searches being a 4-ary heap rather than a binary one: a search may make
  * a decrease-key for every arc it relaxes and a delete-min for every node,
  * and the 4-ary heap, half as tall as the binary one, makes the former
  * cheaper at the price of comparing four children rather than two in the
  * latter. */
 template< LEMONGraph GR, typename V , typename C >
 struct SMSppCapacityScalingTraits
  : public CapacityScalingDefaultTraits< GR , V , C >
 {
  using Heap = QuadHeap< C , RangeMap< int > >;  ///< the heap
  };

 /// CapacityScaling algorithm using the trait with a 4-ary heap
 template< LEMONGraph GR, typename V , typename C >
 class SMSppCapacityScaling : public
  CapacityScaling< GR , V , C , SMSppCapacityScalingTraits< GR , V , C > >
 {
  public:

  SMSppCapacityScaling( const GR& dgp )
   : CapacityScaling< GR , V , C ,
		      SMSppCapacityScalingTraits< GR , V , C > >( dgp ) {}

  ~SMSppCapacityScaling() = default;
  };

 /// CostScaling algorithm using the trait of integer costs
 /** The trait of integer costs is the default one when C is integer; when C
  * is floating-point it makes the large cost type of the algorithm integer
  * as well, which is what keeps it exact, MCFLemonSolver giving it costs
  * that are integer [see MCFLemonSolver::pass_costs()]: with the default
  * trait the large costs would be double, and the rounding of the potentials
  * can make the algorithm read out of its own vectors. */
 template< LEMONGraph GR, typename V , typename C >
 class SMSppCostScaling : public CostScaling< GR , V , C ,
				  CostScalingDefaultTraits< GR , V , C , true > >
 {
  public:
  SMSppCostScaling( const GR& dgp )
   : CostScaling< GR , V , C ,
		  CostScalingDefaultTraits< GR , V , C , true > >( dgp ) {}

  ~SMSppCostScaling() = default;
  };

 /*--------------------------------------------------------------------------*/
 /*------------------------ CLASS MCFLemonSolver ----------------------------*/
 /*--------------------------------------------------------------------------*/
 /*--------------------------- GENERAL NOTES --------------------------------*/
 /*--------------------------------------------------------------------------*/
 /// CDASolver for MCFBlock based on the LEMON project
 /** The MCFLemonSolver implements Solver interface for MCFBlock that represent
  * (Linear) Min-Cost Flow (MCF) problems, using algorithms of LEMON library.
  * Because MCF is a Linear Program it has a(n exact) dual, and therefore
  * MCFLemonSolver implements the CDASolver interface for also giving out dual
  * information.
  *
  * MCFLemonSolver is template over four different types:
  *
  * - GR, which is the type of graph. The possibilities are:
  *
  *   = SmartDigraph is a simple and fast directed graph implementation. It
  *     is quite memory efficient but at the price that it does not support
  *     node and arc deletion.
  *
  *   = MCFListDigraph, an ad-hoc simple extension of ListDigraph, which is
  *     a standard versatile and fast directed graph implementation based on
  *     linked lists that are stored in std::vector structures. This class
  *     provides only linear time counting for nodes and arcs. The original
  *     ListDigraph supports node and arc deletion, and MCFListDigraph extends
  *     this to arc closure and re-opening.
  *
  * - V, which is the type of flows / deficits; typically, double can be used
  *   for maximum compatibility, but int (or even smaller) would yield better
  *   performances;
  *
  * - C, which is the type of arc costs; typically, double can be used for
  *   maximum compatibility, but int (or even smaller) would yield better
  *   performances;
  *
  * - Algo, which is the specific algorithm (itself, template over GR, V, and
  *   C) implemented in the LEMON package. The possibilities are:
  *
  *   = NetworkSimplex implements the primal Network Simplex algorithm for
  *     finding a minimum cost flow. This algorithm is a highly efficient
  *     specialized version of the linear programming simplex method directly
  *     for the minimum cost flow problem.
  *
  *   = CycleCanceling implements three different cycle-canceling algorithms
  *     for finding a minimum cost flow. The most efficient one is the
  *     Cancel-and-tighten algorithm, thus it is the default method. It runs
  *     in strongly polynomial time, but in practice, it is typically orders of
  *     magnitude slower than the scaling algorithms and NetworkSimplex.
  *
  *   = CostScaling implements a cost scaling algorithm that performs
  *     push/augment and relabel operations for finding a minimum cost flow.
  *     It is a highly efficient primal-dual solution method, which can be
  *     viewed as the generalization of the preflow push-relabel algorithm for
  *     the maximum flow problem. It is a polynomial algorithm.
  *
  *   = CapacityScaling implements the capacity scaling version of the
  *     successive shortest path algorithm for finding a minimum cost flow. It
  *     is an efficient dual solution method, which runs in polynomial time.
  *     In special cases it can be more efficient than CostScaling and
  *     NetworkSimplex algorithms.
  *
  *   In general, NetworkSimplex and CostScaling are the fastest
  *   implementations available in LEMON for solving this problem.
  *
  *   Note that scaling-type algorithms may behave in different ways according
  *   to which combination of V and C is used, and there are different
  *   "traits" for this which are another scaling parameter. However, in order
  *   to make MCFLemonSolver class template over always the same number of
  *   template parameters we fix the trait, which is why
  *   SMSppCapacityScaling and SMSppCostScaling are defined (as template over
  *   < GR , V , C >, each with its own trait) that are meant to be used
  *   instead of the original CapacityScaling and CostScaling.
  *
  * abLEMON takes fractional data in the four algorithms, comparing the
  * reduced costs and the flows against tolerances; yet CostScaling and
  * CycleCanceling are given the costs scaled by a power of 2 and rounded,
  * which are integer [see pass_costs()], on which their scaling ends at an
  * exact optimum, and CostScaling, CycleCanceling and CapacityScaling the
  * capacities and the supplies scaled and rounded in the same way; the
  * value and the potentials are those of the true costs. Besides,
  * CycleCanceling and CapacityScaling do not support negative costs on arcs
  * with infinite capacity. The deficits of a MCFBlock that sum to zero only
  * up to the rounding of their computation are balanced before each run,
  * since LEMON would otherwise read the problem as infeasible, whereas
  * deficits that really do not sum to zero make the problem infeasible
  * rather than being taken as demands that may be left unmet. */

template< template< typename , typename , typename > class Algo ,
           LEMONGraph GR, typename V , typename C >
class MCFLemonSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using ThisAlgo = Algo< GR , V , C >;
 using Index = unsigned int;
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;
 typedef typename ThisAlgo::ProblemType ProblemType;

 enum int_par_type_LEMON {
  kKBest = intLastParCDAS ,  ///< how the flows of new_var_solution() come
  intLastParLEMON
  ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

 enum str_par_type_LEMON {
  strDMXFile = strLastParCDAS ,  ///< DMX filename to output the instance
  strLastParLEMON
  ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

 enum LEMON_sol_type { INFEASIBLE , OPTIMAL , UNBOUNDED };

 static constexpr int kErrorStatus = -1;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 MCFLemonSolver( void ) :  um( nullptr ) , cm( nullptr ) , bm( nullptr ) {}

 ~MCFLemonSolver( void ) override { guts_of_destructor(); }

/*--------------------------------------------------------------------------*/
/*-------------------------- PUBLIC METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /// register the Solver to a MCFBlock
 /** Provides the [MCF]Block encoding the problem instance that the
  * Solver will solve. It is entirely implemented in this class because this
  * is done uniformly for all LEMON algorithms. */

 void set_Block( Block* block ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
 /// (try to) solve the MCF encoded in the MCFBlock
 /** Basically invokes the run() method of the underlying LEMON Algo. A lot
  * of the preparatory steps (locking the Block and the Solver, printing the
  * DMX file if required ...) are common to all the Algo and therefore are
  * implemented in this class; a guts_of_compute() method, which the
  * specialised classes implement, is invoked at the right time. */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 double get_elapsed_time( void ) const override { return( this->ticks ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// the optimal value, whichever the current solution is

 OFValue get_lb( void ) override { return( opt_value() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the value of the current solution [see new_var_solution()]

 OFValue get_ub( void ) override {
  if( f_cur >= 0 )
   return( OFValue( f_kbest->totalCost( f_cur ) ) );
  return( opt_value() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int get_status( void ) const { return( this->status ); }

/*--------------------------------------------------------------------------*/

 /// the flows are those of an optimal solution, none is given otherwise

 bool has_var_solution( void ) override {
  return( this->get_status() == ThisAlgo::ProblemType::OPTIMAL );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves to the next extreme flow, in the order of their costs
 /** After the optimal flow that compute() gives, up to intMaxSol - 1 more
  * flows are given, one per call, each of them a vertex of the polyhedron
  * of the flows: first the other optimal ones, if any, then the others in
  * the order of their costs, as long as their cost exceeds the optimal one
  * by at most dblAAccSol and by at most dblRAccSol times the largest
  * between 1 and the absolute value of the optimal one. They are found, the
  * first time the method is called after compute(), by the ranking of the
  * vertices of KBestExtremeFlows [see k_best_extreme_flows.h of abLEMON],
  * from the optimal flow of the algorithm when it is given the data as they
  * are, and solving the problem again with NetworkSimplex when it is given
  * them scaled; a vertex equal to the optimal flow already given is
  * skipped. A flow inside a face of the polyhedron is not given, since its
  * (sub)gradient in a Lagrangian or Benders' decomposition is a convex
  * combination of those of the vertices of the face.
  *
  * The flows come from the cycles of least reduced cost [see CYCLES in
  * KBestExtremeFlows], not from the ranking of the bases, since a
  * degenerate vertex may have very many, as a path in a shortest path
  * problem has as many as the trees containing it: the other optimal
  * vertices are those that the cycles of zero reduced cost lead to, and
  * the others are the best ones one cycle of least reduced cost away from
  * an optimal one, in the order of their costs, but they may not be the
  * best vertices. Each optimal vertex costs O( n m log n ) time; hence 2
  * flows are looked for first, and then twice as many as before whenever
  * they are all given, so that a caller asking for few of intMaxSol flows,
  * as a LagBFunction that sets intMaxSol to the size of its pool does,
  * pays for few. The potentials [see get_dual_solution()] remain the
  * optimal ones. */

 bool new_var_solution( void ) override {
  if( ( f_max_sol <= 1 ) ||
      ( this->get_status() != ThisAlgo::ProblemType::OPTIMAL ) )
   return( false );

  // the flows are found a few at a time, as many more each time as were
  // found before, so that a caller asking for few of many pays for few
  if( f_asked == 0 )
   run_kbest( std::min( f_max_sol , 2 ) );

  for( ;; ) {
   while( ( f_given < f_max_sol ) && ( f_next < f_kbest->flowNum() ) ) {
    const int j = f_next++;
    if( same_as_optimal( j ) )
     continue;
    f_cur = j;
    ++f_given;
    return( true );
    }
   if( ( f_given >= f_max_sol ) || ( f_asked >= f_max_sol ) ||
       ( f_kbest->flowNum() < f_asked ) )
    return( false );
   run_kbest( std::min( 2 * f_asked , f_max_sol ) );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// the potentials are those of an optimal solution, the certificate of
 /// infeasibility being a direction [see has_dual_direction()]

 bool has_dual_solution( void ) override {
  return( this->get_status() == ThisAlgo::ProblemType::OPTIMAL );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves to the next optimal potentials
 /** After the optimal potentials that compute() gives, up to intMaxSol - 1
  * more are given, one per call, each optimal and a vertex of the dual
  * polyhedron: those of the trees of the shortest paths of the reduced
  * costs, from and to each node in turn, of the optimal flow [see
  * runPotentials() of KBestExtremeFlows]; potentials that differ from those
  * already given by a constant are skipped. They are looked for 2 first,
  * and then twice as many as before whenever they are all given, as the
  * flows of new_var_solution() are. The flows [see get_var_solution()]
  * remain the current ones. */

 bool new_dual_solution( void ) override {
  if( ( f_max_sol <= 1 ) ||
      ( this->get_status() != ThisAlgo::ProblemType::OPTIMAL ) )
   return( false );

  if( f_pasked == 0 )
   run_potentials( std::min( f_max_sol , 2 ) );

  for( ;; ) {
   while( ( f_pgiven < f_max_sol ) &&
	  ( f_pnext < f_kbest->potentialNum() ) ) {
    const int j = f_pnext++;
    if( same_as_optimal_potentials( j ) )
     continue;
    f_pcur = j;
    ++f_pgiven;
    return( true );
    }
   if( ( f_pgiven >= f_max_sol ) || ( f_pasked >= f_max_sol ) ||
       ( f_kbest->potentialNum() < f_pasked ) )
    return( false );
   run_potentials( std::min( 2 * f_pasked , f_max_sol ) );
   }
  }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration* solc = nullptr ) override {
  auto MCFB = static_cast< MCFBlock* >( f_Block );

  /* The arc is identified by its LEMON id, which is the position it was
   * created in and therefore the index it has in the MCFBlock: ArcIt does
   * *not* enumerate the arcs in that order, and using the order of the
   * enumeration would silently give the right flows on the wrong arcs. The
   * flows are written all at once, since the dynamic arcs of the MCFBlock
   * are a list that setting them one by one would scan once per arc; the
   * closed arcs, which ArcIt does not enumerate, have flow 0. */
  MCFBlock::Vec_FNumber X( MCFB->get_NArcs() , 0 );
  for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a ) {
   auto i = Index( dgp->id( a ) );
   if( i < X.size() )
    X[ i ] = current_flow( a );
   }
  MCFB->set_x( X.cbegin() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_dual_solution( Configuration* solc = nullptr ) override {
  auto MCFB = static_cast< MCFBlock* >( f_Block );

  // the node is identified by its LEMON id, see get_var_solution(); the
  // digraph may also have more nodes than the MCFBlock has
  MCFBlock::Vec_CNumber Pi( MCFB->get_NNodes() , 0 );
  for( typename GR::NodeIt n( *dgp ) ; n != INVALID ; ++n ) {
   auto i = Index( dgp->id( n ) );
   if( i < Pi.size() )
    Pi[ i ] = current_potential( n );
   }
  MCFB->set_pi( Pi.cbegin() );

  // the reduced costs, i.e., the dual values of the bound constraints, which
  // the potentials alone do not give and which a dual solution of the
  // abstract representation needs, as the :MCFSolver give them [see
  // MCFSolution::write()]; a deleted arc has none
  MCFBlock::Vec_CNumber RC( MCFB->get_NArcs() , 0 );
  for( Index i = 0 ; i < RC.size() ; ++i )
   if( ! MCFB->is_deleted( i ) )
    RC[ i ] = ( MCFB->get_C().empty() ? 0 : MCFB->get_C( i ) )
              + Pi[ MCFB->get_SN( i ) - 1 ]
                               - Pi[ MCFB->get_EN( i ) - 1 ];
  MCFB->set_rc( RC.cbegin() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the current solution as a MCFSolution, out of the Solver's data
 /** Returns the current solution as a MCFSolution [see MCFBlock.h], built out
  * of the data structures of the LEMON algorithm rather than by writing it in
  * the Variable and the Constraint of the MCFBlock and having it read back
  * from there: no abstract representation is therefore required to exist, and
  * the MCFBlock is not written into at all, hence it is not lock()-ed and any
  * number of Solver attached to it can produce their own Solution at the same
  * time.
  *
  * Which parts of the solution are saved, i.e., the flows, the potentials or
  * both, is not decided here: the MCFBlock is asked for an empty MCFSolution
  * with the very same Configuration it would be asked for the full one, and
  * the shape of what it returns is what says it [see
  * MCFBlock::get_Solution()]. A part that is asked for but is not available,
  * as the flows are when the problem is infeasible, is left out rather than
  * being filled with junk.
  *
  * nullptr is returned if neither a solution nor a dual one is available. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override {
  if( ! f_Block )
   return( nullptr );

  // an unbounded instance has a direction to give rather than a solution:
  // the Solution says that what it holds is one, and has no potentials
  if( has_var_direction() ) {
   auto MCFB = static_cast< MCFBlock * >( f_Block );
   auto sol = static_cast< MCFSolution * >( MCFB->get_Solution( solc ,
								 true ) );
   if( ! sol->get_x().empty() ) {
    sol->is_direction( true );
    sol->set_x( cycle_flow() );
    }
   if( ! sol->get_pi().empty() )
    sol->set_pi( MCFBlock::Vec_CNumber() );
   return( sol );
   }

  if( ( ! has_var_solution() ) && ( ! has_dual_solution() ) )
   return( nullptr );

  auto MCFB = static_cast< MCFBlock * >( f_Block );
  auto sol = static_cast< MCFSolution * >( MCFB->get_Solution( solc , true ) );

  if( ! sol->get_x().empty() ) {
   if( has_var_solution() ) {
    MCFBlock::Vec_FNumber X( MCFB->get_NArcs() );
    for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a ) {
     auto i = Index( dgp->id( a ) );
     if( i < X.size() )
      X[ i ] = current_flow( a );
     }
    sol->set_x( std::move( X ) );
    }
   else
    sol->set_x( MCFBlock::Vec_FNumber() );
   }

  if( ! sol->get_pi().empty() ) {
   if( has_dual_solution() ) {
    MCFBlock::Vec_CNumber Pi( MCFB->get_NNodes() );
    for( typename GR::NodeIt n( *dgp ) ; n != INVALID ; ++n ) {
     auto i = Index( dgp->id( n ) );
     if( i < Pi.size() )
      Pi[ i ] = current_potential( n );
     }
    sol->set_pi( std::move( Pi ) );
    }
   else
    sol->set_pi( MCFBlock::Vec_CNumber() );
   }

  return( sol );
  }

/*--------------------------------------------------------------------------*/
 /// the Solution is filled from the data of the LEMON algorithm, not from
 /// the Variable

 [[nodiscard]] bool is_get_Solution_physical( void ) const override {
  return( true );
  }

 /// true if the instance is unbounded: the certificate is a direction
 /** On an unbounded instance the certificate is a cycle of negative cost
  * whose arcs have infinite capacity [see unboundedCycle() of abLEMON],
  * which the algorithm gives, or NetworkSimplex on the true data when the
  * algorithm was given an infinite capacity as a finite bound [see
  * pass_flows()]. With intMaxSol larger than 1 there are up to intMaxSol
  * such cycles, pairwise arc-disjoint and in the order of their mean costs
  * [see unboundedCycles() of abLEMON], given one per call of
  * new_var_direction(), the first call of which looks for them: the cycle
  * of the algorithm, which is the first one, costs nothing more than the
  * run, the others a cycle of minimum mean cost each. */

 bool has_var_direction( void ) override {
  if( this->get_status() != ThisAlgo::ProblemType::UNBOUNDED )
   return( false );
  if( f_dir < 0 ) {
   find_cycles( 1 );
   f_dir = 0;
   }
  return( ! f_dirs.empty() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves to the next cycle [see has_var_direction()]

 bool new_var_direction( void ) override {
  if( ! has_var_direction() )
   return( false );
  if( ( ! f_dirs_all ) && ( f_max_sol > 1 ) ) {
   f_dirs_all = true;
   find_cycles( f_max_sol );
   }
  if( f_dir + 1 >= int( f_dirs.size() ) )
   return( false );
  ++f_dir;
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes one unit of flow along the current cycle in the flow Variable
 /** Writes in the flow Variable of the MCFBlock one unit of flow along the
  * current cycle [see has_var_direction()], and tells the MCFBlock that
  * what they hold is a direction [see MCFBlock::is_direction()]. The
  * Configuration says which part to write, as in MCFBlock::get_Solution():
  * a value of 2 asks for the dual part alone, and then nothing is done.
  * Throws if the last run did not end declaring the instance unbounded. */

 void get_var_direction( Configuration * dirc = nullptr ) override {
  if( ! f_Block )  // no [MCF]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( dirc );
  if( tsolc && ( tsolc->f_value == 2 ) )
   return;

  if( ! has_var_direction() )
   throw( std::logic_error( "MCFLemonSolver::get_var_direction: no "
			    "certificate of unboundedness" ) );

  auto MCFB = static_cast< MCFBlock * >( f_Block );
  MCFB->is_direction( true );  // what the Variable hold is a direction
  MCFB->set_x( cycle_flow().begin() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the instance is infeasible: the certificate is a direction
 /** On an infeasible instance the certificate is a set S of nodes whose
  * supply exceeds the capacity of the arcs leaving it [see infeasibleCut()
  * of abLEMON], or, when the supplies do not sum to zero, the set of all
  * the nodes, whose supply is then positive or negative. It is the dual
  * direction whose potentials are -1 on S and 0 elsewhere, or +1 on all the
  * nodes when their supply is negative, along which the dual objective
  * grows by the violation of the cut [see get_dual_direction()]. The set
  * is that of the algorithm, or of NetworkSimplex on the true data when the
  * algorithm was given the flows scaled [see pass_flows()], and it is
  * checked on the true data. With intMaxSol larger than 1 there are up to
  * intMaxSol such sets, given one per call of new_dual_direction(): first
  * the weakly connected components of the set of the most violated cut
  * whose violation is positive, which are disjoint cuts with no arc between
  * two of them [see disjointInfeasibleCuts() of abLEMON], then the most
  * violated cuts in order [see infeasibleCuts() of abLEMON], looked for at
  * the first call of new_dual_direction(), since they cost up to n maximum
  * flows each while the first set is that of the algorithm. */

 bool has_dual_direction( void ) override {
  if( this->get_status() != ThisAlgo::ProblemType::INFEASIBLE )
   return( false );
  if( f_cut < 0 ) {
   find_cuts( 1 );
   f_cut = 0;
   }
  return( ! f_cuts.empty() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves to the next set [see has_dual_direction()]

 bool new_dual_direction( void ) override {
  if( ! has_dual_direction() )
   return( false );
  if( ( ! f_cuts_all ) && ( f_max_sol > 1 ) ) {
   f_cuts_all = true;
   find_cuts( f_max_sol );
   }
  if( f_cut + 1 >= int( f_cuts.size() ) )
   return( false );
  ++f_cut;
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes the dual direction of the current set in the MCFBlock
 /** Writes in the dual values of the flow conservation constraints of the
  * MCFBlock the potentials of the current set [see has_dual_direction()],
  * \f$ \pi \f$, and in those of the bound constraints the reduced costs
  * of the homogeneous system, \f$ \pi_{SN(a)} - \pi_{EN(a)} \f$, as the
  * ray of the dual is written by MILPSolver.
  * The Configuration says which part to write, as in
  * MCFBlock::get_Solution(): a value of 1 asks for the primal part alone,
  * and then nothing is done. Throws if the last run did not end declaring
  * the instance infeasible. */

 void get_dual_direction( Configuration * dirc = nullptr ) override {
  if( ! f_Block )  // no [MCF]Block to write to
   return;         // cowardly and silently return

  auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( dirc );
  if( tsolc && ( tsolc->f_value == 1 ) )
   return;

  if( ! has_dual_direction() )
   throw( std::logic_error( "MCFLemonSolver::get_dual_direction: no "
			    "certificate of infeasibility" ) );

  auto MCFB = static_cast< MCFBlock * >( f_Block );
  MCFBlock::Vec_CNumber Pi( MCFB->get_NNodes() , 0 );
  for( auto nd : f_cuts[ f_cut ] ) {
   auto i = Index( dgp->id( nd ) );
   if( i < Pi.size() )
    Pi[ i ] = f_cut_sign[ f_cut ];
   }
  MCFB->set_pi( Pi.cbegin() );

  MCFBlock::Vec_CNumber RC( MCFB->get_NArcs() , 0 );
  for( Index i = 0 ; i < RC.size() ; ++i )
   if( ! MCFB->is_deleted( i ) )
    RC[ i ] = Pi[ MCFB->get_SN( i ) - 1 ] - Pi[ MCFB->get_EN( i ) - 1 ];
  MCFB->set_rc( RC.cbegin() );
  }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
 /// gets the @return default value of the string parameter @param par

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override {
  if( par >= strLastParLEMON )
   throw( std::invalid_argument( "MCFLemonSolver::get_dflt_str_par: invalid "
				 "parameter " +
				 std::to_string( par ) ) );

  static const std::string _empty;
  if( par == strDMXFile )
   return( _empty );

  return( CDASolver::get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the @return default value of the double parameter @param par
 /** dblRelAcc, dblFAccSol, dblAAccDSol and dblAbsAcc are 0, i.e., the
  * tolerances are those that abLEMON sets on the data [see
  * pass_tolerances()]. */

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  if( ( par == dblRelAcc ) || ( par == dblFAccSol ) ||
      ( par == dblAAccDSol ) || ( par == dblAbsAcc ) )
   return( 0 );
  return( CDASolver::get_dflt_dbl_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the @return value of the double parameter @param par

 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  if( par == dblRelAcc )
   return( f_rel_cst );
  if( par == dblFAccSol )
   return( f_rel_flw );
  if( par == dblAAccDSol )
   return( f_abs_cst );
  if( par == dblAbsAcc )
   return( f_abs_flw );
  if( par == dblRAccSol )
   return( f_rel_sol );
  if( par == dblAAccSol )
   return( f_abs_sol );
  return( CDASolver::get_dbl_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameter @param par to @param value
 /** The tolerances of the algorithm on the reduced costs and on the flows
  * [see tolerances() of abLEMON] come from dblRelAcc and dblFAccSol, if
  * positive, as relative ones: the former times the largest absolute value
  * of a cost, the latter times the largest one among the finite capacities
  * and the deficits; otherwise from dblAAccDSol and dblAbsAcc, if positive
  * and finite, as absolute ones; otherwise abLEMON sets them on the data,
  * as it does by default [see pass_tolerances()]. */

 void set_par( idx_type par , double value ) override {
  if( par == dblRelAcc )
   f_rel_cst = value;
  else if( par == dblFAccSol )
   f_rel_flw = value;
  else if( par == dblAAccDSol )
   f_abs_cst = value;
  else if( par == dblAbsAcc )
   f_abs_flw = value;
  else if( par == dblRAccSol )
   f_rel_sol = value;
  else if( par == dblAAccSol )
   f_abs_sol = value;
  else
   CDASolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the @return value of the int parameter @param par

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == intMaxSol )
   return( f_max_sol );
  if( par == kKBest )
   return( f_kbest_method );
  return( CDASolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the int parameter @param par to @param value
 /** intMaxSol is the largest number of the flows given [see
  * new_var_solution()], and kKBest says how they are generated by
  * KBestExtremeFlows of abLEMON: 0, the default, by the cycles of least
  * reduced cost from the optimal vertex (CYCLES), 1 as the k best vertices
  * in order (RANKING), exactly when the capacities and the deficits are
  * integer. */

 void set_par( idx_type par , int value ) override {
  if( par == intMaxSol )
   f_max_sol = value;
  else if( par == kKBest ) {
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "MCFLemonSolver::set_par: invalid kKBest " +
				  std::to_string( value ) ) );
   f_kbest_method = value;
   forget_flows();
   }
  else
   CDASolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par == kKBest )
   return( 0 );
  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type
 int_par_str2idx( const std::string & name ) const override {
  if( name == "kKBest" )
   return( kKBest );
  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string &
 int_par_idx2str( idx_type idx ) const override {
  static const std::string name = "kKBest";
  if( idx == kKBest )
   return( name );
  return( CDASolver::int_par_idx2str( idx ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the @return value of string parameter @param par

 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override
 {
  if( par == strDMXFile )
   return( this->f_dmx_file );

  return( CDASolver::get_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the string parameter @param par to @param value

 using CDASolver::set_par;

 void set_par( idx_type par , std::string && value ) override {
  if( par == strDMXFile )
   f_dmx_file = std::move( value );
  else
   CDASolver::set_par( par , std::move( value ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @return number of string algorithmic parameters

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( strLastParLEMON );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  if( name == "strDMXFile" )
   return( strDMXFile );

  return( CDASolver::str_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  if( idx >= strLastParLEMON )
   throw( std::invalid_argument( "MCFLemonSolver::str_par_idx2str: invalid "
				 "index " +
				 std::to_string( idx ) ) );

  static const std::string par = "strDMXFile";
  if( idx == strDMXFile )
   return( par );

  return( CDASolver::str_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 /** The only reason why MCFSolver::add_Modification() needs be defined is to
  * properly react to NBModification. Indeed, the correct reaction is to
  * *immediately* reload the MCFBlock, besides clearing the list of
  * Modification as Solver::add_Modification() already does. The issue is
  * that if arcs/nodes are added/deleted after the NBModification is issued
  * but before it is processed, then the number of nodes/arcs at the moment
  * in which the NBModification is processed is different from that at the
  * moment in which is issued, which may break the "naming convention"
  * (because the name of, say, a newly created arc depends on the current
  * state and/or number of the arcs).
  *
  * The MCFBlock need not be "freshly minted": the graph is rebuilt with the
  * arcs the MCFBlock has closed or deleted closed or deleted in the graph as
  * well, which is what happens when the MCFBlock is a R3Block of another
  * MCFBlock that is loaded and then further modified, the NBModification to
  * this MCFBlock being generated by a map_forward_Modification() of the one
  * to the original MCFBlock. */

 void add_Modification( sp_Mod& mod ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_constructor( void ) { f_algo = nullptr; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void guts_of_destructor( void ) {
  forget_flows();
  delete f_algo;
  f_algo = nullptr;
  delete um;
  um = nullptr;
  delete cm;
  cm = nullptr;
  delete icm;
  icm = nullptr;
  f_cost_scale = 1;
  delete ium;
  ium = nullptr;
  delete ibm;
  ibm = nullptr;
  f_flow_scale = 1;
  delete bm;
  bm = nullptr;
  delete dgp;
  dgp = nullptr;
  }

/*--------------------------------------------------------------------------*/

 /// the optimal value, taken with the true data

 OFValue opt_value( void ) {
  switch( this->get_status() ) {
   case( ThisAlgo::ProblemType::OPTIMAL ) :
    if( ( f_cost_scale != 1 ) || ( f_flow_scale != 1 ) ) {
     // the Algo saw the costs, or the flows, scaled and rounded [see
     // pass_costs() and pass_flows()]: the value is taken with the true
     // costs and flows
     OFValue v = 0;
     for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
      v += OFValue( f_algo->flow( a ) / f_flow_scale ) * (*cm)[ a ];
     return( v );
     }
    return( OFValue( f_algo->totalCost() ) );
   case( ThisAlgo::ProblemType::UNBOUNDED ) :
    return( -Inf< OFValue >() );
   default : return( Inf< OFValue >() );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the flow on the arc a in the current solution [see new_var_solution()]

 [[nodiscard]] double current_flow( typename GR::Arc a ) const {
  if( f_cur >= 0 )
   return( double( f_kbest->flow( f_cur , a ) ) );
  return( double( f_algo->flow( a ) / f_flow_scale ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the j-th extreme flow is the optimal one that compute() gave

 [[nodiscard]] bool same_as_optimal( int j ) const {
  for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a ) {
   const double x = f_algo->flow( a ) / f_flow_scale;
   if( std::abs( double( f_kbest->flow( j , a ) ) - x ) >
       1e-9 * std::max( 1.0 , std::abs( x ) ) )
    return( false );
   }
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void make_kbest( void );

 void run_kbest( int k );

 void run_potentials( int k );

 void find_cycles( int k );

 void find_cuts( int k );

 [[nodiscard]] bool violated( const std::vector< typename GR::Node > & S ,
			      double sign );

 void true_supplies( MCFNodeMapV & s );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// one unit of flow along the current cycle [see has_var_direction()]

 [[nodiscard]] MCFBlock::Vec_FNumber cycle_flow( void ) const {
  auto MCFB = static_cast< MCFBlock * >( f_Block );
  MCFBlock::Vec_FNumber X( MCFB->get_NArcs() , 0 );
  for( auto a : f_dirs[ f_dir ] ) {
   auto i = MCFBlock::Index( dgp->id( a ) );
   if( i < X.size() )
    X[ i ] = 1;
   }
  return( X );
  }

 /// the potential of the node n in the current potentials [see
 /// new_dual_solution()]

 [[nodiscard]] double current_potential( typename GR::Node n ) const {
  if( f_pcur >= 0 )
   return( double( f_kbest->potential( f_pcur , n ) ) );
  return( double( f_algo->potential( n ) / f_cost_scale ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the j-th potentials are those that compute() gave, up to a
 /// constant

 [[nodiscard]] bool same_as_optimal_potentials( int j ) const {
  double lo = Inf< double >() , hi = -Inf< double >() , big = 1;
  for( typename GR::NodeIt n( *dgp ) ; n != INVALID ; ++n ) {
   const double p = f_algo->potential( n ) / f_cost_scale;
   const double d = double( f_kbest->potential( j , n ) ) - p;
   lo = std::min( lo , d );
   hi = std::max( hi , d );
   big = std::max( big , std::abs( p ) );
   }
  return( hi - lo <= 1e-9 * big );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// drops the extreme flows, the cycles and the cuts of the last compute()
 // the method of KBestExtremeFlows that kKBest says
 typename KBestExtremeFlows< GR , V , C >::Method kbest_method( void ) const {
  return( f_kbest_method == 1 ? KBestExtremeFlows< GR , V , C >::RANKING
			       : KBestExtremeFlows< GR , V , C >::CYCLES );
  }

 void forget_flows( void ) {
  f_dirs.clear();
  f_dir = -1;
  f_dirs_all = false;
  f_cuts.clear();
  f_cut_sign.clear();
  f_cut = -1;
  f_cuts_all = false;
  f_kbest.reset();
  f_asked = 0;
  f_next = 0;
  f_cur = -1;
  f_given = 1;
  f_pasked = 0;
  f_pnext = 0;
  f_pcur = -1;
  f_pgiven = 1;
  }

 void close_supplies( MCFNodeMapV & s , typename GR::Node big , V sum );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual void guts_of_compute( void ) = 0;

 void guts_of_poM( c_p_Mod mod );

 void guts_of_set_Block( MCFBlock* MCFB );

 void pass_tolerances( void );

 void pass_costs( void );

 void pass_flows( typename GR::Node big );

 /// true if the Algo is given the capacities and supplies scaled
 [[nodiscard]] bool flows_scaled( void ) const {
  return( f_flow_scale != 1 );
  }

 /*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

  std::string f_dmx_file;  ///< string for DMX file output

  Algo< GR , V , C >* f_algo = nullptr;  ///< the actual LEMON algo

  ProblemType status = ProblemType::INFEASIBLE;  ///< return status

  GR * dgp = nullptr;
  ///< the (di)graph, i.e., either MCFListDigraph or SmartDigraph

  bool cost_changed = false;    ///< true if a cost Modification is done
  bool cap_changed = false;     ///< true if a capacity Modification is done
  bool supply_changed = false;  ///< true if a supply Modification is done

  int n_arcs = 0;
  ///< number of arcs in the graph, closed arcs count as deleted
  int n_arcs_added = 0;
                   ///< number of arcs added to the graph with addArc/openArc
  int n_arcs_deleted = 0;
              ///< number of arcs deleted from the graph with rmvArc/closeArc

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------- PRIVATE FIELDS -------------------------------*/

  double ticks;  ///< Elapsed time in ticks for compute() method

  MCFArcMapV * um;   ///< ArcMap that contains the capacity of each arc
  MCFArcMapV * cm;   ///< ArcMap that contains the cost of each arc
  MCFNodeMapV * bm;  /// NodeMap that contains the supply of each node

  MCFArcMapV * icm = nullptr;
  double f_rel_cst = 0;  ///< relative tolerance of the costs [dblRelAcc]
 double f_rel_flw = 0;  ///< relative tolerance of the flows [dblFAccSol]
 double f_abs_cst = 0;  ///< absolute tolerance of the costs [dblAAccDSol]
 double f_abs_flw = 0;  ///< absolute tolerance of the flows [dblAbsAcc]

 int f_max_sol = 1;  ///< largest number of the flows given [intMaxSol]
 int f_kbest_method = 0;  ///< how the flows are generated [kKBest]
 double f_rel_sol = Inf< double >();  ///< relative gap of them [dblRAccSol]
 double f_abs_sol = Inf< double >();  ///< absolute gap of them [dblAAccSol]

 std::unique_ptr< KBestExtremeFlows< GR , V , C > > f_kbest;
 ///< the extreme flows [see new_var_solution()], made at its first call
 int f_next = 0;   ///< the next of them to look at
 int f_cur = -1;   ///< the current one, -1 for the flow of the Algo
 int f_pasked = 0;  ///< the number of the potentials last looked for
 int f_pnext = 0;   ///< the next of them to look at
 int f_pcur = -1;   ///< the current one, -1 for those of the Algo
 int f_pgiven = 1;  ///< the number of the potentials given since compute()
 int f_given = 1;  ///< the number of the flows given since compute()
 int f_asked = 0;  ///< the number of the flows the last ranking looked for
 C f_gap = -1;     ///< the gap of the ranking, -1 for none

 std::vector< std::vector< typename GR::Arc > > f_dirs;
 ///< the cycles of an unbounded instance [see has_var_direction()]
 int f_dir = -1;   ///< the current one, -1 before they are looked for
 bool f_dirs_all = false;  ///< true if more than the first were looked for

 std::vector< std::vector< typename GR::Node > > f_cuts;
 ///< the sets of an infeasible instance [see has_dual_direction()]
 std::vector< double > f_cut_sign;  ///< their potentials, -1 or +1
 int f_cut = -1;   ///< the current one, -1 before they are looked for
 bool f_cuts_all = false;  ///< true if more than the first were looked for

 ///< the costs scaled by f_cost_scale and rounded [see pass_costs()]
  double f_cost_scale = 1;
  ///< the scale of icm, 1 when the Algo is given the costs as they are

  MCFArcMapV * ium = nullptr;
  ///< the capacities scaled by f_flow_scale and rounded [see pass_flows()]
  MCFNodeMapV * ibm = nullptr;
  ///< the supplies scaled by f_flow_scale and rounded [see pass_flows()]
  double f_flow_scale = 1;
  ///< the scale of ium and ibm, 1 when the Algo is given the flows as they are
  V f_inf_cap = 0;
  ///< the scaled bound given for an infinite capacity [see pass_flows()]

  /// true if the Algo is given the capacities and supplies scaled and rounded
  static constexpr bool f_scaled_flows = std::is_floating_point< V >::value &&
   ( std::is_same< Algo< GR , V , C > , SMSppCostScaling< GR , V , C > >::value
     || std::is_same< Algo< GR , V , C > ,
                      CycleCanceling< GR , V , C > >::value
     || std::is_same< Algo< GR , V , C > ,
                      SMSppCapacityScaling< GR , V , C > >::value );

  /// true if the Algo bounds an infinite capacity in the vector of the
  /// capacities it is given, which a later run would read [see compute()]
  static constexpr bool f_bounds_inf_cap =
   std::is_same< Algo< GR , V , C > , SMSppCostScaling< GR , V , C > >::value
   || std::is_same< Algo< GR , V , C > , CycleCanceling< GR , V , C > >::value;

/*--------------------------------------------------------------------------*/

 }; // end( class MCFLemonSolver< Algo , GR , C , V > )

/*--------------------------------------------------------------------------*/
/*------------------------- SPECIALIZED CLASSES ----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------- MCFLemonSolverNetworkSimplex ------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverNetworkSimplex< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing NetworkSimplex specific algorithmic parameters and methods
 * set/get_*_par for manage them.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to NetworkSimplex. */

template< typename GR , typename V , typename C >
class MCFLemonSolverNetworkSimplex : public
 MCFLemonSolver< NetworkSimplex , GR , V , C >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 // large batch of "using", unfortunately needed due to the fact that on
 // the first pass of compiling a template, where only non-dependent names
 // are looked up, the compiler does not "see" the base class MCFLemonSolver
 // and all its ancestors (CDASolver, Solver, ThinComputeInterface) and all
 // their names

 using BaseClass = MCFLemonSolver< NetworkSimplex , GR , V , C >;
 using idx_type = ThinComputeInterface::idx_type;

 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::intLastParLEMON;
 using BaseClass::status;
 using BaseClass::strLastParLEMON;
 using CDASolver::kBlockLocked;
 using CDASolver::kInfeasible;
 using CDASolver::kStopTime;
 using Solver::f_Block;
 using Solver::f_id;
 using Solver::lock;
 using Solver::OFValue;
 using Solver::unlock;
 using ThinComputeInterface::kUnEval;

 using typename BaseClass::ThisAlgo;
 using NSPivotRule = typename NetworkSimplex< GR , V , C >::PivotRule;

/*--------------------------------------------------------------------------*/
 // enums for handling the extra parameters

 enum LEMON_NS_int_par_type {
  kPivot = intLastParLEMON ,  ///< pivot algorithm for the simplex
  kReopt ,                   ///< whether to re-optimize from the last basis
  intLastParLEMON_NS  ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

/*------- CONSTRUCTING AND DESTRUCTING MCFLemonSolverNetworkSimplex --------*/

 /// constructor: initializes algorithm parameters
 /** Void constructor. Define f_pivot_rule to the default algorithm
  * parameters used by NetworkSimplex, then calls guts_of_constructor().  */

 MCFLemonSolverNetworkSimplex( void ) : BaseClass()
 {
  BaseClass::guts_of_constructor();
  f_pivot_rule = NSPivotRule::BLOCK_SEARCH;
  f_reopt = true;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverNetworkSimplex() override { BaseClass::guts_of_destructor(); }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 /// set the int parameter  @param par with @param value

 void set_par( idx_type par , int value ) override {
  if( par == kPivot ) {
   if( ( value < 0 ) || ( value > 4 ) )
    throw( std::invalid_argument( "MCFLemonSolverNetworkSimplex::"
				 "set_par: invalid kPivot " +
				  std::to_string( value ) ) );

   if( value == f_pivot_rule )
    return; // nothing is changed

   f_pivot_rule = NSPivotRule( value );
   return;
   }

  if( par == kReopt ) {
   f_reopt = ( value != 0 );
   return;
   }

  BaseClass::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_NS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_NS )
   throw( std::invalid_argument( "MCFLemonSolverNetworkSimplex::"
				 "get_dflt_int_par: invalid parameter " +
				 std::to_string( par ) ) );
  if( par == kPivot )
   return( NSPivotRule::BLOCK_SEARCH );
  if( par == kReopt )
   return( 1 );

  return( BaseClass::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kPivot )
   return( f_pivot_rule );
  if( par == kReopt )
   return( f_reopt );

  return( BaseClass::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string& name )
  const override {
  if( name == "kPivot" )
   return( kPivot );
  if( name == "kReopt" )
   return( kReopt );

  return( BaseClass::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_NS )
   throw( std::invalid_argument( "MCFLemonSolverNetworkSimplex::"
				 "int_par_idx2str: invalid index " +
				 std::to_string( idx ) ) );

  static const std::array< std::string , 2 > pars = { "kPivot" , "kReopt" };
  if( ( idx >= kPivot ) && ( idx < intLastParLEMON_NS ) )
   return( pars[ idx - kPivot ] );

  return( BaseClass::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 /* With kReopt the network simplex is asked to re-optimize: the basis of the
  * previous run survives a change of the costs as it is, and a change of the
  * capacities or of the supplies after the flow of its tree has been
  * recomputed, while the algorithm drops it by itself when the graph changes
  * [see runWarm() in the shim of the build]; without, each run starts from
  * the artificial basis. */

 void guts_of_compute( void ) override {
  status = f_reopt ? f_algo->runWarm( NSPivotRule( f_pivot_rule ) )
                   : f_algo->run( NSPivotRule( f_pivot_rule ) );
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

  NSPivotRule f_pivot_rule;
  bool f_reopt;  ///< whether to re-optimize from the last basis [see kReopt]

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class MCFLemonSolverNetworkSimplex< GR , V , C >  )

/*--------------------------------------------------------------------------*/
/*--------------------- MCFLemonSolverCycleCanceling -----------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCycleCanceling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing CycleCanceling specific algorithmic parameters and methods
 * set/get_*_par for manage them.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to CycleCanceling. */

template< typename GR , typename V , typename C >
class MCFLemonSolverCycleCanceling : public
 MCFLemonSolver< CycleCanceling , GR , V , C >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 // large batch of "using", unfortunately needed due to the fact that on
 // the first pass of compiling a template, where only non-dependent names
 // are looked up, the compiler does not "see" the base class MCFLemonSolver
 // and all its ancestors (CDASolver, Solver, ThinComputeInterface) and all
 // their names

 using BaseClass = MCFLemonSolver< CycleCanceling , GR , V , C >;
 using idx_type = ThinComputeInterface::idx_type;

 using BaseClass::intLastParLEMON;

 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::status;
 using BaseClass::strLastParLEMON;
 using CDASolver::dblLastParCDAS;
 using CDASolver::kBlockLocked;
 using CDASolver::kInfeasible;
 using CDASolver::kStopTime;
 using Solver::f_Block;
 using Solver::f_id;
 using Solver::lock;
 using Solver::OFValue;
 using Solver::unlock;
 using ThinComputeInterface::kUnEval;

 using typename BaseClass::ThisAlgo;
 using CCMethod = typename ThisAlgo::Method;

 enum LEMON_CC_int_par_type {
  kMethod = intLastParLEMON ,  ///< the method parameter of CycleCanceling
  kReopt ,                    ///< whether to re-optimize from the last run
  intLastParLEMON_CC  ///< first allowed parameter value for derived classes
   /**< convenience value for easily allow derived classes
    * to further extend the set of types of return codes */
  };

/*----- CONSTRUCTING AND DESTRUCTING MCFLemonSolverCycleCanceling ----------*/

 /// constructor: initializes algorithm parameters
 /** Void constructor. Define the f_method of CycleCanceling, then calls
  * guts_of_constructor().  */

 MCFLemonSolverCycleCanceling( void ) : BaseClass() {
  BaseClass::guts_of_constructor();
  f_method = CycleCanceling< GR , V , C >::Method::CANCEL_AND_TIGHTEN;
  f_reopt = true;
  }

 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverCycleCanceling( void ) override {
  BaseClass::guts_of_destructor();
  }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == kMethod ) {
   if( ( value < 0 ) || ( value > 2 ) )
    throw( std::invalid_argument( "MCFLemonSolverCycleCanceling::"
				 "set_par: invalid kMethod " +
				  std::to_string( value ) ) );
   if( value == f_method )
    return; // nothing is changed

   f_method = CCMethod( value );
   return;
   }

  if( par == kReopt ) {
   f_reopt = ( value != 0 );
   return;
   }

  BaseClass::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_CC );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_CC )
   throw( std::invalid_argument( "MCFLemonSolverCycleCanceling::"
				 "get_dflt_int_par: invalid parameter " +
				 std::to_string( par ) ) );
  if( par == kMethod )
   return( CycleCanceling< GR , V , C >::Method::CANCEL_AND_TIGHTEN );
  if( par == kReopt )
   return( 1 );

  return( BaseClass::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kMethod )
   return( f_method );
  if( par == kReopt )
   return( f_reopt );

  return( BaseClass::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "kMethod" )
   return( kMethod );
  if( name == "kReopt" )
   return( kReopt );

  return( BaseClass::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_CC )
   throw( std::invalid_argument( "MCFLemonSolverCycleCanceling::"
				 "int_par_idx2str: invalid parameter " +
				 std::to_string( idx ) ) );

  static const std::array< std::string , 2 > pars = { "kMethod" , "kReopt" };
  if( ( idx >= kMethod ) && ( idx < intLastParLEMON_CC ) )
   return( pars[ idx - kMethod ] );

  return( BaseClass::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 /* With kReopt the negative cycles are canceled from the flow of the last
  * run, which after a change of the capacities or of the supplies is put
  * within the new bounds and made to meet the new supplies [see runWarm()
  * in the shim of the build]; without, each run finds a feasible flow
  * anew. */

 void guts_of_compute( void ) override {
  status = f_reopt ? f_algo->runWarm( CCMethod( f_method ) )
                   : f_algo->run( CCMethod( f_method ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------- PRIVATE FIELDS -------------------------------*/

 CCMethod f_method;
 bool f_reopt;  ///< whether to re-optimize from the last run [see kReopt]

/*--------------------------------------------------------------------------*/

 };  // end( class MCFLemonSolverCycleCanceling< GR , V , C > )

/*--------------------------------------------------------------------------*/
/*--------------------- MCFLemonSolverCapacityScaling ----------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCapacityScaling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, the enum
 * indexing the algorithmic parameter of CapacityScaling and the set/get_*_par
 * methods that manage it.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to CapacityScaling.
 *
 * kFactor is the scaling factor of the successive approximations, i.e., the
 * base of the geometric sequence of the deltas: it has to be at least 2, the
 * default is 4, and 1 would ask for the plain successive shortest path
 * algorithm, which is why CapacityScaling refuses it. */

template< typename GR , typename V , typename C >
 class MCFLemonSolverCapacityScaling : public
 MCFLemonSolver< SMSppCapacityScaling , GR , V , C >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 // large batch of "using", unfortunately needed due to the fact that on
 // the first pass of compiling a template, where only non-dependent names
 // are looked up, the compiler does not "see" the base class MCFLemonSolver
 // and all its ancestors (CDASolver, Solver, ThinComputeInterface) and all
 // their names

 using BaseClass = MCFLemonSolver< SMSppCapacityScaling , GR , V , C >;
 using idx_type = ThinComputeInterface::idx_type;

 using BaseClass::intLastParLEMON;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::status;
 using BaseClass::strLastParLEMON;
 using CDASolver::kBlockLocked;
 using CDASolver::kInfeasible;
 using CDASolver::kStopTime;
 using Solver::f_Block;
 using Solver::f_id;
 using Solver::lock;
 using Solver::OFValue;
 using Solver::unlock;
 using ThinComputeInterface::kUnEval;

/*--------------------------------------------------------------------------*/
 // enums for handling the extra parameters

 enum LEMON_CapS_int_par_type {
  kFactor = intLastParLEMON ,  ///< the scaling factor of CapacityScaling
  kReopt ,                    ///< whether to re-optimize from the last run
  intLastParLEMON_CapS  ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

/*------ CONSTRUCTING AND DESTRUCTING MCFLemonSolverCapacityScaling --------*/

 /// constructor: initializes algorithm parameters
 /** Void constructor. Sets f_factor to the default scaling factor of
  * CapacityScaling, then calls guts_of_constructor(). */

 MCFLemonSolverCapacityScaling( void ) : BaseClass() {
  BaseClass::guts_of_constructor();
  f_factor = 4;
  f_reopt = true;
  }

 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverCapacityScaling( void ) override {
  BaseClass::guts_of_destructor();
  }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == kFactor ) {
   if( value < 2 )
    throw( std::invalid_argument( "MCFLemonSolverCapacityScaling::"
				  "set_par: invalid kFactor " +
				  std::to_string( value ) ) );
   f_factor = value;
   return;
   }

  if( par == kReopt ) {
   f_reopt = ( value != 0 );
   return;
   }

  BaseClass::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_CapS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_CapS )
   throw( std::invalid_argument( "MCFLemonSolverCapacityScaling::"
				 "get_dflt_int_par: invalid parameter " +
				 std::to_string( par ) ) );
  if( par == kFactor )
   return( 4 );

  if( par == kReopt )
   return( 1 );

  return( BaseClass::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kFactor )
   return( f_factor );

  if( par == kReopt )
   return( f_reopt );

  return( BaseClass::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "kFactor" )
   return( kFactor );

  if( name == "kReopt" )
   return( kReopt );

  return( BaseClass::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_CapS )
   throw( std::invalid_argument( "MCFLemonSolverCapacityScaling::"
				 "int_par_idx2str: invalid parameter " +
				 std::to_string( idx ) ) );

  static const std::array< std::string , 2 > pars = { "kFactor" , "kReopt" };
  if( ( idx >= kFactor ) && ( idx < intLastParLEMON_CapS ) )
   return( pars[ idx - kFactor ] );

  return( BaseClass::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 /* With kReopt the algorithm starts from the flow and the potentials of the
  * last run, the flow put within the new bounds and at the one the reduced
  * costs say [see runWarm() in the shim of the build]; without, each run
  * starts from the zero flow and the zero potentials. */

 void guts_of_compute( void ) override {
  status = f_reopt ? f_algo->runWarm( f_factor ) : f_algo->run( f_factor );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------- PRIVATE FIELDS -------------------------------*/

 int f_factor;  ///< the scaling factor, at least 2
 bool f_reopt;  ///< whether to re-optimize from the last run [see kReopt]

/*--------------------------------------------------------------------------*/

 };  // end( class MCFLemonSolverCapacityScaling< GR, V, C> )

/*--------------------------------------------------------------------------*/
/*----------------------------- CostScaling --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCostScaling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing CostScaling specific algorithmic parameters and methods
 * set/get_*_par for manage them.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to CostScaling.
 *
 * kMethod takes the values of CostScaling::Method, i.e., 0 = PUSH,
 * 1 = AUGMENT and 2 = PARTIAL_AUGMENT, plus kAutoMethod = 3, the default,
 * which is PARTIAL_AUGMENT when the capacities and the supplies are given as
 * they are and AUGMENT when they are scaled [see pass_flows()]: the number of
 * operations of PARTIAL_AUGMENT grows with the magnitude of the flows, and on
 * scaled ones it can be 2 to 3 orders of magnitude slower than the other
 * two. */

template< typename GR , typename V , typename C >
class MCFLemonSolverCostScaling : public
 MCFLemonSolver< SMSppCostScaling , GR , V , C >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 // large batch of "using", unfortunately needed due to the fact that on
 // the first pass of compiling a template, where only non-dependent names
 // are looked up, the compiler does not "see" the base class MCFLemonSolver
 // and all its ancestors (CDASolver, Solver, ThinComputeInterface) and all
 // their names

 using BaseClass = MCFLemonSolver< SMSppCostScaling , GR , V , C >;
 using idx_type = ThinComputeInterface::idx_type;

 using BaseClass::intLastParLEMON;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::flows_scaled;
 using BaseClass::status;
 using BaseClass::strLastParLEMON;
 using CDASolver::kBlockLocked;
 using CDASolver::kInfeasible;
 using CDASolver::kStopTime;
 using Solver::f_Block;
 using Solver::f_id;
 using Solver::lock;
 using Solver::OFValue;
 using Solver::unlock;
 using ThinComputeInterface::kUnEval;

 using typename BaseClass::ThisAlgo;
 using CSMethod = typename ThisAlgo::Method;

 /// the value of kMethod that chooses the method [see the class]
 static constexpr int kAutoMethod = 3;

 enum LEMON_CS_int_par_type {
  kMethod = intLastParLEMON ,  ///< the method of cost scaling
  kReopt ,                    ///< whether to re-optimize from the last run
  intLastParLEMON_CS
  //< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

/*-------- CONSTRUCTING AND DESTRUCTING MCFLemonSolverCostScaling ----------*/

 /// constructor: initializes algorithm parameters
 /** Void constructor. Define the f_method of CostScaling, then calls
  * guts_of_constructor().  */

 MCFLemonSolverCostScaling( void ) : BaseClass() {
  BaseClass::guts_of_constructor();
  f_method = kAutoMethod;
  f_reopt = true;
  }

 /// destructor, calls guts_of_destructor()
 ~MCFLemonSolverCostScaling( void ) override {
  BaseClass::guts_of_destructor();
  }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == kMethod ) {
   if( ( value < 0 ) || ( value > kAutoMethod ) )
    throw( std::invalid_argument( "MCFLemonSolverCostScaling::"
				 "set_par: invalid kMethod " +
				  std::to_string( value ) ) );
   f_method = value;
   return;
   }

  if( par == kReopt ) {
   f_reopt = ( value != 0 );
   return;
   }

  BaseClass::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_CS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_CS )
   throw( std::invalid_argument( "MCFLemonSolverCostScaling::"
				 "get_dflt_int_par: invalid parameter " +
				 std::to_string( par ) ) );
  if( par == kMethod )
   return( kAutoMethod );
  if( par == kReopt )
   return( 1 );

  return( BaseClass::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kMethod )
   return( f_method );
  if( par == kReopt )
   return( f_reopt );

  return( BaseClass::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "kReopt" )
   return( kReopt );

  return( name == "kMethod" ? kMethod : BaseClass::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_CS )
   throw( std::invalid_argument( "MCFLemonSolverCostScaling::"
				 "int_par_idx2str: invalid parameter " +
				 std::to_string( idx ) ) );

  static const std::array< std::string , 2 > pars = { "kMethod" , "kReopt" };
  if( ( idx >= kMethod ) && ( idx < intLastParLEMON_CS ) )
   return( pars[ idx - kMethod ] );

  return( BaseClass::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 /* With kReopt the scaling starts from the flow and the potentials of the
  * last run, the flow put within the new bounds and made to meet the new
  * supplies after a change of them, and from the largest violation of the
  * optimality conditions under them [see runWarm() in the shim of the
  * build]; without, each run finds a feasible flow anew and starts from the
  * largest cost. */

 void guts_of_compute( void ) override {
  const CSMethod m = f_method != kAutoMethod ? CSMethod( f_method ) :
                     flows_scaled() ? ThisAlgo::AUGMENT :
                                      ThisAlgo::PARTIAL_AUGMENT;
  status = f_reopt ? f_algo->runWarm( m ) : f_algo->run( m );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------- PRIVATE FIELDS -------------------------------*/

 int f_method;  ///< a CSMethod, or kAutoMethod
 bool f_reopt;  ///< whether to re-optimize from the last run [see kReopt]

 };  // end( class MCFLemonSolverCostScaling< GR , V , C> )

 /** @} end( group( MCFLemonSolver_CLASSES ) ) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 } // namespace SMSpp_di_unipi_it

/*--------------------------------------------------------------------------*/

#endif /* __MCFLemonSolver */

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFLemonSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
