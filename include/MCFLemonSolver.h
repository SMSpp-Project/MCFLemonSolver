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

#include <lemon/list_graph.h>

#include <lemon/smart_graph.h>

#include <lemon/dimacs.h>

#include <type_traits>

#include <chrono>

#include <lemon/lgf_writer.h>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/
/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace lemon;
 using namespace lemon::concepts;
 using namespace std;

 /** @} ---------------------------------------------------------------------*/
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
  *   re-added yet (the name of th last arc if there are no such arcs).
  *   This is not so in the original ListDigraph, so the implementation
  *   of addArc() has to be changed accordingly.
  *
  * - The need of implementing the openArc() and closeArc() methods, that
  *   realise the arc losing /opening operations of MCFBlock. A "closed" arc
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
  class MCFArc : ListDigraph::Arc
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

 bool isClosed( int n ) { return( arcs[ n ].source < 0 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// add a new arc with given endpoints

 MCFArc addArc( int u , int v )
 {
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

   return( MCFArc( minidx ) );
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

  return( MCFArc( n ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// close the given arc

 void closeArc( int n )
 {
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
  // just do the little that closeArc() did not
  arcs[ n ].next_in = first_free_arc;
  first_free_arc = n;
  arcs[ n ].prev_in = -2;
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
  *   for maximum compatibility, but int (or even smaller) would yeld better
  *   performances;
  *
  * - C, which is the type of ar costs; typically, double can be used for
  *   maximum compatibility, but int (or even smaller) would yeld better
  *   performances.
  *
  * Furthermore, scaling-type algorithms may behave in different ways
  * according to which combination of V and C is used, and there are different
  * "traits" for this which are another scaling parameter. However, we prefer
  * that the MCFLemonSolver class is always template over the three first
  * parameter only, which is why we define SMSppCapacityScaling and
  * SMSppCostScaling as template over < GR , V , C > and using the default
  * trait.
  *
  * Thus, the concept LEMONGraph is defined that only allows all possible
  * types of LEMON graphs, i.e., SmartDigraph and MCFListDigraph.
  *  @{ */

 /// concept for "one of the LEMON graphs"
 template< typename Type >
 concept LEMONGraph = std::is_same< Type , SmartDigraph >::value ||
                      std::is_same< Type , MCFListDigraph >::value;

 /// CapacityScaling algorithm using the default trait
 template< LEMONGraph GR, typename V , typename C >
 class SMSppCapacityScaling : public CapacityScaling< GR , V , C >
 {
  public:

  SMSppCapacityScaling( const GR& dgp ) : CapacityScaling< GR , V , C >( dgp ) {}

  ~SMSppCapacityScaling() = default;
  };

 /// CostScaling algorithm using the default trait
 template< LEMONGraph GR, typename V , typename C >
 class SMSppCostScaling : public CostScaling< GR , V , C >
 {
  public:
  SMSppCostScaling( const GR& dgp ) : CostScaling< GR , V , C >( dgp ) {}

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
  *   for maximum compatibility, but int (or even smaller) would yeld better
  *   performances;
  *
  * - C, which is the type of ar costs; typically, double can be used for
  *   maximum compatibility, but int (or even smaller) would yeld better
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
  *     for finding a minimum cost flow. The most efficent one is the
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
  *   template parameters we fix the use of the default trait, which is why
  *   SMSppCapacityScaling and SMSppCostScaling are defined (as template over
  *   < GR , V , C >) that are meant to be used instead of the original
  *   CapacityScaling and CostScaling. */

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

 ~MCFLemonSolver( void ) override { delete um; delete cm; delete bm; }

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
  * implemented in this class; a guts_of_compute() method is inkoked at the
  * right time that eed be implemented in specialised classe. */

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 double get_elapsed_time( void ) const override { return( this->ticks ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_lb( void ) override { return( OFValue( f_algo->totalCost() ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OFValue get_ub( void ) override { return( OFValue( f_algo->totalCost() ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int get_status( void ) const { return( this->status ); }

/*--------------------------------------------------------------------------*/

 bool has_var_solution( void ) override {
  switch( this->get_status() ) {
   case( ThisAlgo::ProblemType::OPTIMAL ) :
   case( ThisAlgo::ProblemType::UNBOUNDED ) :
    return( true );
   default : return( false );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_dual_solution( void ) override {
  switch( this->get_status() ) {
   case( ThisAlgo::ProblemType::OPTIMAL ) :
   case( ThisAlgo::ProblemType::INFEASIBLE ) :
    return( true );
   default : return( false );
   }
  }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration* solc = nullptr ) override {
  auto MCFB = static_cast< MCFBlock* >( f_Block );

  /* The arc is identified by its LEMON id, which is the position it was
   * created in and therefore the index it has in the MCFBlock: ArcIt does
   * *not* enumerate the arcs in that order, and using the order of the
   * enumeration would silently give the right flows on the wrong arcs. */
  for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a ) {
   auto i = Index( dgp->id( a ) );
   if( i < MCFB->get_NArcs() )
    MCFB->set_x( i , f_algo->flow( a ) );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_dual_solution( Configuration* solc = nullptr ) override {
  auto MCFB = static_cast< MCFBlock* >( f_Block );

  // the node is identified by its LEMON id, see get_var_solution(); the
  // digraph may also have more nodes than the MCFBlock has
  for( typename GR::NodeIt n( *dgp ) ; n != INVALID ; ++n ) {
   auto i = Index( dgp->id( n ) );
   if( i < MCFB->get_NNodes() )    // NB: set_pi() takes the value first
    MCFB->set_pi( f_algo->potential( n ) , i );
   }
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
      X[ i ] = f_algo->flow( a );
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
      Pi[ i ] = f_algo->potential( n );
     }
    sol->set_pi( std::move( Pi ) );
    }
   else
    sol->set_pi( MCFBlock::Vec_CNumber() );
   }

  return( sol );
  }

/*--------------------------------------------------------------------------*/
 /// returns false until we understand if and how LEMON does is

 bool has_var_direction( void ) override { return( false ); }

 void get_var_direction( Configuration* dirc = nullptr ) override {
  throw( std::logic_error( "LEMONSolver:get_dual_direction() called" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns false until we understand if and how LEMON does is

 bool has_dual_direction( void ) override { return( false ); }

 void get_dual_direction( Configuration* dirc = nullptr ) override {
  throw( std::logic_error( "LEMONSolver:get_dual_direction() called" ) );
  }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
 /// gets the @return default value of the string parameter @param par

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override {
  if( par > strLastParLEMON )
   throw( std::invalid_argument( "Invalid str parameter: out_of_range " +
				 std::to_string( par ) ) );

  static const std::string _empty;
  if( par == strLastParLEMON )
   return( _empty );

  return( CDASolver::get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the @return value of string parameter @param par

 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override
 {
  if( par == strDMXFile )
   return( this->f_dmx_file );

  return( get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @return number of string algorithimc parameters

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( strLastParLEMON );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  if( name == "strDMXFile" )
   return( strDMXFile );

  return( CDASolver::str_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  if( idx > strLastParLEMON )
   throw( std::invalid_argument( "str_par_idx2str: index out_of_range " +
				 std::to_string( idx ) ) );

  static const std::string par = "strDMXFile";
  if( idx == strDMXFile )
   return( par );

  return( CDASolver::dbl_par_idx2str( idx ) );
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
  * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
  * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
  *
  * This should ordinarily always happen, as whenever the MCFBlock is changed
  * the NBModification is immediately issued. The problem may come if the
  * MCFBlock is a R3Block of another MCFBlock which is loaded and then
  * further modified, and the NBModification to this MCFBlock is generated by
  * a map_forward_Modification() of the NBModification to the original
  * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
  * closed or deleted arcs and this method would not work. */

 void add_Modification( sp_Mod& mod ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 void guts_of_constructor( void ) { f_algo = nullptr; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void guts_of_destructor( void ) {
  delete f_algo;
  dgp->clear();
  delete dgp;
  }

/*--------------------------------------------------------------------------*/

 virtual void guts_of_compute( void ) = 0;

 void guts_of_poM( c_p_Mod mod );

 void guts_of_set_Block( MCFBlock* MCFB );

 /*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ---------------------------*/
/*--------------------------------------------------------------------------*/

  std::string f_dmx_file;  ///< string for DMX file output

  Algo< GR , V , C >* f_algo;  ///< the (pointer to) the actual LEMON algo

  ProblemType status = ProblemType::INFEASIBLE;  ///< return status

  GR * dgp;  ///< the (di)graph, i.e., either MCFListDigraph or SmartDigraph

  bool cost_changed = false;    ///< true if a cost Modification is done
  bool cap_changed = false;     ///< true if a capacity Modification is done
  bool supply_changed = false;  ///< true if a supply Modification is done

  int n_arcs;  ///< number of arcs in the graph, closed arcs count as deleted
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

/*--------------------------------------------------------------------------*/

 }; // end( class MCFLemonSolver< Algo , GR , C , V > )

/*--------------------------------------------------------------------------*/
/*------------------------- SPECIALIZED CLASSES ----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------- MCFLemonSolverNetworkSimplex ------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverNetworkSimplex< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing NetworkSimplex specific algorithimc parameters and methods
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

 using BaseClass::f_algo;
 using BaseClass::intLastParCDAS;
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
  kPivot = intLastParCDAS ,  ///< pivot algorithm for the simplex
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
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverNetworkSimplex() override { BaseClass::guts_of_destructor(); }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 /// set the int parameter  @param par with @param value

 void set_par( idx_type par , int value ) override {
  if( par == kPivot ) {
   if( ( value < 0 ) || ( value > 4 ) )
    throw( std::invalid_argument( "Error: invalid kPivot " +
				  std::to_string( value ) ) );

   if( value == f_pivot_rule )
    return; // nothing is changed

   f_pivot_rule = NSPivotRule( value );
   return;
   }

  CDASolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_NS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_NS )
   throw( std::invalid_argument( "Invalid int parameter: out_of_range " +
				 std::to_string( par ) ) );
  if( par == kPivot )
   return( NSPivotRule::BLOCK_SEARCH );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kPivot )
   return( f_pivot_rule );

  return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string& name )
  const override {
  if( name == "kPivot" )
   return( kPivot );

  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_NS )
   throw( std::invalid_argument( "int_par_idx2str: index out_of_range " +
				 std::to_string( idx ) ) );

  static const std::string par = "kPivot";
  if( idx == kPivot )
   return( par );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void guts_of_compute( void ) override {
  status = f_algo->run( NSPivotRule( f_pivot_rule ) );
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

  NSPivotRule f_pivot_rule;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 }; // end( class MCFLemonSolverNetworkSimplex< GR , V , C >  )

/*--------------------------------------------------------------------------*/
/*--------------------- MCFLemonSolverCycleCanceling -----------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCycleCanceling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing CycleCanceling specific algorithimc parameters and methods
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

 using BaseClass::intLastParCDAS;

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
  kMethod = intLastParCDAS ,  ///< the method parameter of CycleCanceling
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
  }

 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverCycleCanceling( void ) override { BaseClass::guts_of_destructor(); }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == kMethod ) {
   if( ( value < 0 ) || ( value > 4 ) )
    throw( std::invalid_argument( "Error: invalid kMethod " +
				  std::to_string( value ) ) );
   if( value == f_method )
    return; // nothing is changed

   f_method = CCMethod( value );
   return;
   }

  CDASolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_CC );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_CC )
   throw( std::invalid_argument( "Invalid int parameter: out_of_range " +
				 std::to_string( par ) ) );
  if( par == kMethod )
   return( CycleCanceling< GR , V , C >::Method::CANCEL_AND_TIGHTEN );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kMethod )
   return( f_method );

  return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "kMethod" )
   return( kMethod );

  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_CC )
   throw( std::invalid_argument( "Invalid int parameter: out_of_range " +
				 std::to_string( idx ) ) );

  static const std::string par = "kMethod";
  if( idx == kMethod )
   return( par );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void guts_of_compute( void ) override {
  status = f_algo->run( CCMethod( f_method ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------- PRIVATE FIELDS -------------------------------*/

 CCMethod f_method;

/*--------------------------------------------------------------------------*/

 };  // end( class MCFLemonSolverCycleCanceling< GR , V , C > )

/*--------------------------------------------------------------------------*/
/*--------------------- MCFLemonSolverCapacityScaling ----------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCapacityScaling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to CapacityScaling. */

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

 using BaseClass::intLastParCDAS;
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

/*------ CONSTRUCTING AND DESTRUCTING MCFLemonSolverCapacityScaling --------*/

 /// constructor, calls guts_of_constructor()

 MCFLemonSolverCapacityScaling( void ) : BaseClass() {
  BaseClass::guts_of_constructor();
  }

 /// destructor, calls guts_of_destructor()

 ~MCFLemonSolverCapacityScaling( void ) override {
  BaseClass::guts_of_destructor();
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void guts_of_compute( void ) override { status = f_algo->run(); }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class MCFLemonSolverCapacityScaling< GR, V, C> )

/*--------------------------------------------------------------------------*/
/*----------------------------- CostScaling --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolverCostScaling< GR , V , C > that derives from
 * MCFLemonSolver and contains the specialized compute() method, enums for
 * indexing CostScaling specific algorithimc parameters and methods
 * set/get_*_par for manage them.
 *
 * The template parameters are the same as those of MCFLemonSolver, except
 * of course the first that is fixed to CostScaling. */

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

 using BaseClass::intLastParCDAS;
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

 using typename BaseClass::ThisAlgo;
 using CSMethod = typename ThisAlgo::Method;

 enum LEMON_CS_int_par_type {
  kMethod = intLastParCDAS ,  ///< the method of cost scaling
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
  f_method = SMSppCostScaling< GR , V , C >::Method::PARTIAL_AUGMENT;
  }

 /// destructor, calls guts_of_destructor()
 ~MCFLemonSolverCostScaling( void ) override { BaseClass::guts_of_constructor(); }

/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/

 void set_par( idx_type par , int value ) override {
  if( par == kMethod ) {
   if( ( value < 0 ) || ( value > 4 ) )
    throw( std::invalid_argument( "Error: invalid kMethod " +
				  std::to_string( value ) ) );
   if( value == f_method )
    return; // nothing is changed

   f_method = CSMethod( value );
   return;
   }

  CDASolver::set_par( par , value );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( intLastParLEMON_CS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  if( par > intLastParLEMON_CS )
   throw( std::invalid_argument( "Invalid int parameter: out_of_range " +
				 std::to_string( par ) ) );
  if( par == kMethod )
   return( SMSppCostScaling< GR , V , C >::Method::PARTIAL_AUGMENT );

  return( CDASolver::get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == kMethod )
   return( f_method );

  return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  return( name == "kMethod" ? kMethod : CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  if( idx > intLastParLEMON_CS )
   throw( std::invalid_argument( "Invalid int parameter: out_of_range " +
				 std::to_string( idx ) ) );

  static const std::string par = "kMethod";
  if( idx == kMethod )
   return( par );

  return( CDASolver::int_par_idx2str( idx ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void guts_of_compute( void ) override {
  status = f_algo->run( CSMethod( f_method ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*-------------------------- PRIVATE METHODS -------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------- PRIVATE FIELDS -------------------------------*/

 CSMethod f_method;

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
