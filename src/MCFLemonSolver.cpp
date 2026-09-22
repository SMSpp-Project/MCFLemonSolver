/*--------------------------------------------------------------------------*/
/*------------------------- File MCFLemonSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MCFLemonSolver class.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFLemonSolver.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_define_force_load( MCFLemonSolver );

// insert into the Solver factory all he 8 variants of MCFLemonSolver having
// the types C and V and all possible choice of the 4 algorithms and 2 graphs

#define SMSpp_insert_LEMON( C , V ) \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverNetworkSimplex< SmartDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverNetworkSimplex< MCFListDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCycleCanceling< SmartDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCycleCanceling< MCFListDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCostScaling< SmartDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCostScaling< MCFListDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCapacityScaling< SmartDigraph , C , V > ); \
 \
SMSpp_insert_in_factory_cpp_0_t(    \
 MCFLemonSolverCapacityScaling< MCFListDigraph , C , V > );

/*--------------------------------------------------------------------------*/
// this macro, coded bit-wise, decides which variants of the algorithms are
// inserted in the factory in terms of types of allowed costs and flows:
//
// - bit 0 ( SMSpp_which_insert_LEMON & 1 ): double costs and/or flows
//
// - bit 1 ( SMSpp_which_insert_LEMON & 2 ): long costs and/or flows
//
// - bit 2 ( SMSpp_which_insert_LEMON & 4 ): int costs and/or flows
//
// The macro can be either defined from outside (say, the makefile) or,
// if not, set by changing the definition below
//
// Note that with only one bit set to 1, the 8 variants having only that type
// as flows and costs are inserted. With two bits set to 1, 8 variants are
// inserted for each of the 4 possible combinations of the two types as
// flows and costs (i.e., 32 variants). With all three bits set to 1, 8
// variants are inserted for each of the 9 possible combinations of the three
// types as flows and costs (i.e., 72 variants)

#ifndef SMSpp_which_insert_LEMON
 #define SMSpp_which_insert_LEMON 3
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/
// register the various MCFLemonSolver*< GR , V , C > to the Solver factory
// it's all combinations of:
// - four algorithms ( NetworkSimplex , CycleCanceling , CostScaling ,
//   CapacityScaling )
// - two graphs ( SmartDigraph , ListDigraph )
// - three data types ( int , long , double ) in two places ( costs , flows ),
//   i.e., 3^2 = 9 possible data configurations
// so it is 4 * 2 * 9 = 72 versions

#if SMSpp_which_insert_LEMON & 1     // double costs and flows are allowed

 SMSpp_insert_LEMON( double , double )

 #if SMSpp_which_insert_LEMON & 2    // long costs and flows are allowed

  SMSpp_insert_LEMON( long , double )

  SMSpp_insert_LEMON( double , long )

 #endif

 #if SMSpp_which_insert_LEMON & 4    // int costs and flows are allowed

  SMSpp_insert_LEMON( int , double )

  SMSpp_insert_LEMON( double , int )

 #endif
#endif

#if SMSpp_which_insert_LEMON & 2     // long costs and flows are allowed

 SMSpp_insert_LEMON( long , long )

 #if SMSpp_which_insert_LEMON & 4    // int costs and flows are allowed

  SMSpp_insert_LEMON( int , long )

  SMSpp_insert_LEMON( long , int )

 #endif
#endif

#if SMSpp_which_insert_LEMON & 4     // int costs and flows are allowed

  SMSpp_insert_LEMON( int , int )

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of MCFLemonSolver ------------------------*/
/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename > class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::guts_of_set_Block( MCFBlock * MCFB )
{
 // whatever was built for the previous MCFBlock (if any) goes
 guts_of_destructor();

 n_arcs_added = n_arcs_deleted = 0;
 cost_changed = cap_changed = supply_changed = false;

 // create the new Graph (MCFListDigraph or SmartDigraph)
 dgp = new GR;

 // actually reserve get_MaxNNodes() node space in dgp
 auto maxNNodes = MCFB->get_MaxNNodes();
 dgp->reserveNode( maxNNodes );
 MCFBlock::Index n = MCFB->get_NNodes();
 // add all the nodes, one by one (LEMON has no way of adding them all at
 // once)
 for( MCFBlock::Index i = 0 ; i < n ; ++i )
  dgp->addNode();

 // reserve get_MaxNArcs() arcs space
 auto maxNArcs = MCFB->get_MaxNArcs();
 dgp->reserveArc( maxNArcs );
 MCFBlock::Index m = MCFB->get_NArcs();
 n_arcs = m;
 // get starting node subset and ending node subset
 MCFBlock::c_Subset& sn = MCFB->get_SN();
 MCFBlock::c_Subset& en = MCFB->get_EN();

 // add all arcs in dgp
 // sn[ i ] - 1 and en[ i ] - 1 are used because sn, en nodes start from 1
 for( MCFBlock::Index i = 0 ; i < m ; ++i )
  if constexpr( std::is_same< GR , SmartDigraph >::value )
   dgp->addArc( dgp->nodeFromId( sn[ i ] - 1 ) ,
                dgp->nodeFromId( en[ i ] - 1 ) );
  else
   dgp->addArc( sn[ i ] - 1 , en[ i ] - 1 );

 // defining names for types for readability
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;

 /* Now we are going to fill up ArcMap and NodeMap. The three maps exist
  * even when the MCFBlock has no data for them, holding what the empty
  * vector means (+Inf capacities, zero costs and deficits): a Modification
  * can fill the data in later, and a change of the graph passes all three
  * to the Algo again. */
 um = new MCFArcMapV( *dgp , Inf< V >() );  // create the upper map
 if( ! MCFB->get_U().empty() ) {
  auto& u = MCFB->get_U();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   um->set( dgp->arcFromId( i ) , u[ i ] );
  }

 cm = new MCFArcMapV( *dgp , 0 );  // create the cost map
 if( ! MCFB->get_C().empty() ) {
  auto& c = MCFB->get_C();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )  // a deleted arc costs NaN
   cm->set( dgp->arcFromId( i ) , std::isnan( c[ i ] ) ? 0 : c[ i ] );
  }

 bm = new MCFNodeMapV( *dgp , 0 );  // create the supply map
 if( ! MCFB->get_B().empty() ) {
  auto& b = MCFB->get_B();
  for( MCFBlock::Index i = 0 ; i < n ; i++ )
   bm->set( dgp->nodeFromId( i ) ,
            b[ i ] == 0 ? 0 : -b[ i ] );  // note that supply = - deficit
  }

 /* The MCFBlock may come with arcs already closed or deleted: MCFListDigraph
  * closes and erases them as the Modification would have done, so that they
  * can be re-opened or re-created later, whereas SmartDigraph, which can do
  * neither, gives them zero capacity, which is the same for the flow. An
  * arc is closed by fixing its flow Variable, so there can be closed arcs
  * only if the Variable have been generated. */
 const bool hasvars = ( MCFB->get_number_static_variables() > 0 ) ||
                      ( MCFB->get_number_dynamic_variables() > 0 );
 for( MCFBlock::Index i = 0 ; i < m ; ++i ) {
  if( ( ! MCFB->is_deleted( i ) ) && ! ( hasvars && MCFB->is_closed( i ) ) )
   continue;

  if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
   if( MCFB->is_deleted( i ) )
    dgp->erase( dgp->arcFromId( i ) );
   else
    dgp->closeArc( i );
   }
  else
   um->set( dgp->arcFromId( i ) , 0 );

  --n_arcs;
  }

 // new instance of LEMON Algorithm with no arc mixing, on the final graph
 if constexpr( std::is_same< Algo< GR , V , C > ,
	                     NetworkSimplex< GR , V , C > >::value )
  f_algo = new Algo< GR , V , C >( *dgp , false );
 else
  f_algo = new Algo< GR , V , C >( *dgp );

 f_algo->upperMap( *um );   // pass the maps to the Algo
 pass_costs();
 f_algo->supplyMap( *bm );

 // the scaled capacities and supplies are passed by compute()
 if constexpr( f_scaled_flows )
  cap_changed = supply_changed = true;
 }  // end( guts_of_set_Block )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename > class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::pass_costs( void )
{
 /* CostScaling, CycleCanceling and CapacityScaling are exact only on
  * integer costs: on fractional ones the first may read out of its own
  * vectors, the second end away from the optimum and the third report as
  * infeasible a problem that is not. They are therefore given the costs
  * scaled by a power of 2 and rounded, which are integer in C as long as
  * they stay below 2^52. CostScaling multiplies them by ( n + 1 ) * 16, its
  * default factor, in an integer large cost type [see SMSppCostScaling],
  * and its potentials grow up to about n times that, hence the scale also
  * keeps max |c| * 16 * n^2 below 2^62; CycleCanceling sums them along
  * cycles, and CapacityScaling along the shortest paths of its potentials,
  * of at most n arcs, hence the scale keeps max |c| * 16 * n below 2^52.
  * The value and the potentials are brought back to the true costs [see
  * get_ub()]. */
 constexpr bool scaled = std::is_floating_point< C >::value &&
  ( std::is_same< Algo< GR , V , C > , SMSppCostScaling< GR , V , C > >::value
    || std::is_same< Algo< GR , V , C > ,
                     CycleCanceling< GR , V , C > >::value
    || std::is_same< Algo< GR , V , C > ,
                     SMSppCapacityScaling< GR , V , C > >::value );

 if constexpr( ! scaled )
  f_algo->costMap( *cm );
 else {
  double max = 1;
  for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
   max = std::max( max , std::abs( double( (*cm)[ a ] ) ) );

  const double n = countNodes( *dgp ) + 1;
  double bound;
  if constexpr( std::is_same< Algo< GR , V , C > ,
		              SMSppCostScaling< GR , V , C > >::value )
   bound = std::min( std::exp2( 52 ) / max ,
		     std::exp2( 62 ) / ( max * 16 * n * n ) );
  else
   bound = std::exp2( 52 ) / ( max * n * 16 );
  f_cost_scale = std::exp2( std::floor( std::log2( bound ) ) );
  if( ! icm )
   icm = new typename GR::template ArcMap< V >( *dgp , 0 );
  for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
   icm->set( a , std::round( (*cm)[ a ] * f_cost_scale ) );

  f_algo->costMap( *icm );
  }
 }  // end( pass_costs )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename > class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::pass_flows( typename GR::Node big )
{
 /* CostScaling, CycleCanceling and CapacityScaling are exact only on integer
  * capacities and supplies: on fractional ones the residual capacities and
  * the excesses are left with a few ulp that keep arcs and nodes active, and
  * they may report infeasibility or not terminate. They are therefore given
  * the capacities and the supplies scaled by a power of 2 and rounded, which
  * are integer in V as long as the sum of n of them stays below 2^52, the
  * excess of a node being at most the sum of the supplies. The rounded
  * supplies are then made to sum to zero exactly, in integer arithmetic, on
  * the node whose supply is largest in absolute value, the problem having
  * already been found balanced up to rounding. The flows are brought back to
  * the true scale [see get_ub()]. Integer capacities and supplies whose sum
  * of n stays below 2^52 are given as they are, with scale 1: a large scale
  * makes the default method of CostScaling slow [see
  * MCFLemonSolverCostScaling].
  *
  * The three algorithms also saturate every arc of negative cost, and answer
  * "unbounded" when its capacity is infinite, even if the problem has a
  * finite optimum. Every capacity is therefore given as at most a bound that
  * the flow of no arc exceeds in an optimal flow of least total flow, which
  * only has cycles of negative cost, each of them through an arc of negative
  * cost and one of finite capacity: the sum of the positive supplies, plus
  * that of the finite capacities of the arcs of negative cost, or of all the
  * finite capacities if an arc of negative cost has an infinite one, plus
  * one. If the flow of an arc of infinite capacity reaches it the problem is
  * unbounded [see compute()], and otherwise the optimum is that of the true
  * capacities. The bound also keeps a large capacity that no optimal flow
  * uses from setting the scale, and with it the rounding of the other data;
  * it depends on the signs of the costs, hence it is computed anew when they
  * change. */
 V supply = 0;
 for( typename GR::NodeIt nd( *dgp ) ; nd != INVALID ; ++nd )
  if( (*bm)[ nd ] > 0 )
   supply += (*bm)[ nd ];

 V negative = 0;
 V finite = 0;
 bool inf_negative = false;
 for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
  if( (*um)[ a ] < Inf< V >() ) {
   finite += (*um)[ a ];
   if( (*cm)[ a ] < 0 )
    negative += (*um)[ a ];
   }
  else
   if( (*cm)[ a ] < 0 )
    inf_negative = true;

 const V bound = supply + 1 + ( inf_negative ? finite : negative );

 const double n = countNodes( *dgp ) + 1;
 bool integer = double( bound ) * n <= std::exp2( 52 );
 for( typename GR::ArcIt a( *dgp ) ; integer && ( a != INVALID ) ; ++a )
  if( ( (*um)[ a ] < bound ) && ( (*um)[ a ] != std::floor( (*um)[ a ] ) ) )
   integer = false;
 for( typename GR::NodeIt nd( *dgp ) ; integer && ( nd != INVALID ) ; ++nd )
  if( (*bm)[ nd ] != std::floor( (*bm)[ nd ] ) )
   integer = false;

 f_flow_scale = integer ? 1 :
  std::exp2( std::floor( std::log2( std::exp2( 52 ) /
				    ( double( bound ) * n ) ) ) );
 if( ! ium )
  ium = new typename GR::template ArcMap< V >( *dgp , 0 );
 if( ! ibm )
  ibm = new typename GR::template NodeMap< V >( *dgp , 0 );

 f_inf_cap = std::ceil( bound * f_flow_scale );
 for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
  ium->set( a , (*um)[ a ] < bound ?
	        V( std::round( (*um)[ a ] * f_flow_scale ) ) : f_inf_cap );

 V sum = 0;
 for( typename GR::NodeIt nd( *dgp ) ; nd != INVALID ; ++nd ) {
  const V b = std::round( (*bm)[ nd ] * f_flow_scale );
  ibm->set( nd , b );
  sum += b;
  }
 ibm->set( big , (*ibm)[ big ] - sum );

 f_algo->upperMap( *ium );
 f_algo->supplyMap( *ibm );
 }  // end( pass_flows )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename> class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::set_Block( Block * block )
{
 if( block == f_Block ) // actually doing nothing
  return; // cowardly and silently return

 Solver::set_Block( block ); // attach to the new Block

 if( ! block ) // this is just: go sit down in a corner and wait
  return; // all done

 auto MCFB = dynamic_cast< MCFBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument(
                 "MCFLemonSolver::set_Block: block must be a MCFBlock" ) );

 bool owned = MCFB->is_owned_by( f_id );
 if( ( ! owned ) && ( ! MCFB->read_lock() ) )
  throw( std::logic_error( "MCFLemonSolver::set_Block: cannot acquire "
			   "read_lock on MCFBlock" ) );

 guts_of_set_Block( MCFB );  // fill up the graph and f_algo

 if( ! owned )
  MCFB->read_unlock();

 } // end( MCFLemonSolver< Algo , GR , V , C >set_Block )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename> class Algo ,
	  LEMONGraph GR , typename V , typename C >
int MCFLemonSolver< Algo , GR , V , C >::compute( bool changedvars )
{
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;

 const static std::array< int , 3 > LEMONstatus_2_sol_type = {
  Solver::kInfeasible , Solver::kOK , Solver::kUnbounded
  };

 lock(); // first of all, acquire self-lock

 if( ! f_Block ) { // there is no [MCFBlock] to solve
  unlock();
  return( kBlockLocked ); // return error
  }

 bool owned = f_Block->is_owned_by( f_id ); // check if already locked
 if( ( ! owned ) && ( ! f_Block->read_lock() ) ) { // if not, try to read_lock
  unlock();
  return( kBlockLocked ); // return error on failure
  }

 // while [read_]locked, process any outstanding Modification
 process_outstanding_Modification();

 if( ! f_dmx_file.empty() ) {  // if so required
  // output the current instance (after the changes) to a DMX file: the
  // MCFBlock writes it, in the complete DIMACS format and with the closed
  // and the deleted arcs where they are, i.e., the data the graph has
  std::ofstream ProbFile( f_dmx_file ,
                          std::ios_base::out | std::ios_base::trunc );
  if( ! ProbFile.is_open() )
   throw( std::logic_error( "MCFLemonSolver::compute: cannot open DMX file "
			    + f_dmx_file ) );
  f_Block->print( ProbFile , 'C' );
  ProbFile.close();
  }

 if( ! owned ) // if the [MCF]Block was actually read_locked
  f_Block->read_unlock(); // read_unlock it

 // if arcs were added\removed then f_algo must be reset() 
 if( ( n_arcs_added > 0 ) || ( n_arcs_deleted > 0 ) ) {
  n_arcs += n_arcs_added;
  n_arcs -= n_arcs_deleted;
  n_arcs_added = 0;
  n_arcs_deleted = 0;

  f_algo->reset();

  cost_changed = true;
  cap_changed = true;
  supply_changed = true;
  }

 if( cost_changed ) {
  pass_costs();
  cost_changed = false;
  // the bound of the capacities depends on the signs of the costs [see
  // pass_flows()]
  if constexpr( f_scaled_flows )
   cap_changed = true;
  }

 if( ( ! f_scaled_flows ) && cap_changed ) {
  f_algo->upperMap( *um );
  cap_changed = false;
  }

 /* The LEMON algorithms decide on the sum of the supplies taken exactly: a
  * positive sum, even by one ulp, makes the problem infeasible, and a
  * negative one is read as demands that may be left unmet, i.e., another
  * problem. The deficits of a MCFBlock sum to zero only up to the rounding
  * of the values they are computed from, so a residual within the rounding
  * is moved, for this run only, onto the node whose supply is largest in
  * absolute value, while a larger one makes the problem infeasible, as it is
  * for the MCFBlock. Since the moved residual is itself rounded, the supply
  * of that node is then lowered until the sum, taken in the order of the
  * nodes of the graph as LEMON takes it, is not positive. */
 auto sum_of = [ this ]() {
  V s = 0;
  for( typename GR::NodeIt nd( *dgp ) ; nd != INVALID ; ++nd )
   s += (*bm)[ nd ];
  return( s );
  };

 const V sum = sum_of();
 V size = 0;
 typename GR::Node big = INVALID;
 for( typename GR::NodeIt nd( *dgp ) ; nd != INVALID ; ++nd ) {
  size += std::abs( (*bm)[ nd ] );
  if( ( big == INVALID ) ||
      ( std::abs( (*bm)[ nd ] ) > std::abs( (*bm)[ big ] ) ) )
   big = nd;
  }

 const bool balanced = std::abs( sum ) <= 1e-9 * std::max( V( 1 ) , size );

 if constexpr( f_scaled_flows ) {
  // the capacities and the supplies go scaled and rounded, and the latter
  // sum to zero exactly [see pass_flows()]
  if( balanced && ( cap_changed || supply_changed ) ) {
   pass_flows( big );
   cap_changed = supply_changed = false;
   }
  }
 else if( balanced && ( sum != 0 ) ) {
  const V supply = (*bm)[ big ];
  bm->set( big , supply - sum );
  if constexpr( std::is_floating_point< V >::value )
   while( sum_of() > 0 )
    bm->set( big , std::nextafter( (*bm)[ big ] , -Inf< V >() ) );
  f_algo->supplyMap( *bm );
  bm->set( big , supply );
  supply_changed = true;  // the next run passes the map as it is again
  }
 else
  if( supply_changed ) {
   f_algo->supplyMap( *bm );
   supply_changed = false;
   }

 auto start = std::chrono::system_clock::now();
 if( balanced )
  guts_of_compute(); // here the actual magic is done by specialised classes
 else
  status = ThisAlgo::INFEASIBLE;

 // an infinite capacity was given as a finite bound [see pass_flows()]: an
 // arc whose flow reaches it makes the problem unbounded
 if constexpr( f_scaled_flows )
  if( status == ThisAlgo::OPTIMAL )
   for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a )
    if( ( (*um)[ a ] == Inf< V >() ) && ( f_algo->flow( a ) >= f_inf_cap ) ) {
     status = ThisAlgo::UNBOUNDED;
     break;
     }

 /* The rounding of scaled capacities and supplies perturbs them by up to
  * half the inverse of the scale, which is tiny with respect to the supplies
  * unless the bound of the capacities is many orders of magnitude larger
  * than them [see pass_flows()]: then the flow is not one of the true data,
  * and rather than giving it out the run ends in error. */
 if constexpr( f_scaled_flows )
  if( ( status == ThisAlgo::OPTIMAL ) && ( f_flow_scale != 1 ) ) {
   const V tol = 1e-9 * std::max( V( 1 ) , size );
   MCFNodeMapV excess( *dgp , 0 );
   for( typename GR::ArcIt a( *dgp ) ; a != INVALID ; ++a ) {
    const V x = f_algo->flow( a ) / f_flow_scale;
    if( x > (*um)[ a ] + tol )
     status = ProblemType( kErrorStatus );
    excess.set( dgp->source( a ) , excess[ dgp->source( a ) ] - x );
    excess.set( dgp->target( a ) , excess[ dgp->target( a ) ] + x );
    }
   for( typename GR::NodeIt nd( *dgp ) ; nd != INVALID ; ++nd )
    if( std::abs( excess[ nd ] + (*bm)[ nd ] ) > tol )
     status = ProblemType( kErrorStatus );
   }

 auto end = std::chrono::system_clock::now();

 std::chrono::duration< double > elapsed = end - start;
 ticks = elapsed.count();

 unlock(); // release self-lock

 // now give out the result
 if( this->get_status() == kErrorStatus )
  return( Solver::kError );
 return( LEMONstatus_2_sol_type[ this->get_status() ] );

 }  // end( MCFLemonSolver< Algo , GR , V , C >::compute )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename> class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::add_Modification( sp_Mod & mod )
{
 if( std::dynamic_pointer_cast< const NBModification >( mod ) ) {
  // this is the "nuclear option": the MCFBlock has been re-loaded, so
  // the LEMON solver also has to (immediately)
  auto MCFB = static_cast< MCFBlock * >( f_Block );
  guts_of_set_Block( MCFB );
  // besides, any outstanding modification makes no sense any longer
  mod_clear();
  }
 else
  push_back( mod );
 }

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename> class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V ,
                     C >::process_outstanding_Modification( void )
{
 // no-frills loop: do them in order, with no attempt at optimizing
 // note that NBModification have already been dealt with and therefore need
 // not be considered here

 for( ;; ) {
  auto mod = pop();
  if( ! mod )
   break;

  guts_of_poM( mod.get() );
  }
 }  // end( MCFLemonSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

template< template < typename , typename , typename> class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::guts_of_poM( c_p_Mod mod )
{
 auto MCFB = static_cast< MCFBlock * >( f_Block );

 // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* This requires patiently sifting through the possible Modification types
  * to find what this Modification exactly is, and call the appropriate
  * method of LEMON. */

 // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const GroupModification* >( mod ) ) {
  for( const auto& submod : tmod->sub_Modifications() )
   guts_of_poM( submod.get() );

  return;
  }

 // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Note: in the following, we can assume that C, B and U are nonempty. This
  * is because they can be empty only if they are so when the object is
  * loaded. But if a Modification has been issued, they are no longer empty (a
  * Modification changing nothing from the "empty" state is not issued). */

 if( auto tmod = dynamic_cast< const MCFBlockRngdMod* >( mod ) ) {
  auto rng = tmod->rng();

  switch( tmod->type() ) {
   // changing costs - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgCost ) : {
    cost_changed = true;
    for( MCFBlock::Index i = rng.first ; i < rng.second ; i++ )
     if( dgp->valid( dgp->arcFromId( i ) ) ) {
      auto cost = MCFB->get_C( i );
      cm->set( dgp->arcFromId( i ) , std::isnan( cost ) ? 0 : cost );
      }

    return;
    }

   // changing capacities- - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgCaps ) : {
    cap_changed = true;
    for( MCFBlock::Index i = rng.first ; i < rng.second ; i++ )
     if( dgp->valid( dgp->arcFromId( i ) ) )
      um->set( dgp->arcFromId( i ) , MCFB->get_U( i ) );

    return;
    }

   // changing deficits- - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgDfct ) : {
    supply_changed = true;
    for( MCFBlock::Index i = rng.first ; i < rng.second ; i++ )
     bm->set( dgp->nodeFromId( i ) , - MCFB->get_B( i ) );

    return;
    }

   // opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eOpenArc ) : {
    // only if the graph is ListDigraph
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     for( ; rng.first < rng.second ; ++rng.first )
      if( ( ! MCFB->is_deleted( rng.first ) ) &&
	  dgp->valid( dgp->arcFromId( rng.first ) ) &&
	  dgp->isClosed( rng.first ) ) {
       dgp->openArc( rng.first );
       n_arcs_added++;
       }
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   // closing arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eCloseArc ) : {
    // only if the graph is ListDigraph
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     for( ; rng.first < rng.second ; ++rng.first )
      if( ( ! MCFB->is_deleted( rng.first ) ) &&
	  ! dgp->isClosed( rng.first ) &&
	  dgp->valid( dgp->arcFromId( rng.first ) ) ) {
       dgp->closeArc( rng.first );
       n_arcs_deleted++;
       }
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   // adding arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eAddArc ) : {
    // only if the graph is ListDigraph
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     auto ca = MCFB->get_C( rng.first );
     auto startNode = MCFB->get_SN( rng.first ) - 1;
     auto endNode = MCFB->get_EN( rng.first ) - 1;
     auto capacity = MCFB->get_U( rng.first );

     // the graph must give the new arc the name the MCFBlock gave it
     if( MCFBlock::Index( dgp->id( dgp->addArc( startNode , endNode ) ) )
	 != rng.first )
      throw( std::logic_error( "MCFLemonSolver::guts_of_poM: new arc " +
			       std::to_string( rng.first ) +
			       " has a different name in the graph" ) );
     auto arc = dgp->arcFromId( rng.first );

     n_arcs_added++;
     cm->set( arc , std::isnan( ca ) ? 0 : ca );
     um->set( arc , capacity );

     cost_changed = true;
     cap_changed = true;
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   // removing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eRmvArc ) : {
    // only if ListDigraph is used
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     // if the arc was previously deleted, return
     if( ! dgp->valid( dgp->arcFromId( rng.second - 1 ) ) )
      return;

     auto arc = dgp->arcFromId( rng.second - 1 );

     // if the arc isClosed, only 'mark' it as eliminated, otherwise
     // normally erase it
     if( dgp->isClosed( rng.second - 1 ) )
      dgp->eraseClosed( rng.second - 1 );
     else
      dgp->erase( arc );

     n_arcs_deleted++;
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   default :
    throw( std::invalid_argument(
		  "MCFLemonSolver::guts_of_poM: unknown MCFBlockRngdMod type" ) );

   }  // end( case )
  }  // end( MCFBlockRngdMod )

 // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const MCFBlockSbstMod* >( mod ) ) {
  switch( tmod->type() ) {

   // changing costs - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgCost ) : {
    cost_changed = true;
    for( auto i : tmod->nms() )
     if( dgp->valid( dgp->arcFromId( i ) ) ) {
      auto arc = dgp->arcFromId( i );
      auto cost = MCFB->get_C( i );
      cm->set( arc , std::isnan( cost ) ? 0 : cost );
      }

    return;
    }

   // changing capacities- - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgCaps ) : {
    cap_changed = true;
    auto & CC = MCFB->get_C();
    auto & U = MCFB->get_U();
    for( auto i : tmod->nms() )
     if( ( CC.empty() || ! std::isnan( CC[ i ] ) ) &&
	 dgp->valid( dgp->arcFromId( i ) ) )
      um->set( dgp->arcFromId( i ) , U[ i ] );

    return;
    }

   // changing deficits- - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eChgDfct ) : {
    supply_changed = true;
    auto & B = MCFB->get_B();
    for( auto i : tmod->nms() )
     bm->set( dgp->nodeFromId( i ) , -B[ i ] );

    return;
    }

   // opening arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eOpenArc ) : {
    // only if ListDigraph is used
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     for( auto arc : tmod->nms() )
      // checking if the arc is closed
      if( ( ! MCFB->is_deleted( arc ) ) && dgp->isClosed( arc ) &&
	  dgp->valid( dgp->arcFromId( arc ) ) ) {
       dgp->openArc( arc );
       n_arcs_added++;
       }
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   // closing arcs - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eCloseArc ) : {
    // only if ListDigraph is used
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     for( auto arc : tmod->nms() )
      // checking if the arc is not closed
      if( ( ! MCFB->is_deleted( arc ) ) &&
	  dgp->valid( dgp->arcFromId( arc ) ) &&
	  ! dgp->isClosed( arc ) ) {
       dgp->closeArc( arc );
       n_arcs_deleted++;
       }
     }
    else
     throw( std::logic_error( "MCFLemonSolver::guts_of_poM: SmartDigraph "
			      "does not support operations on arcs" ) );
    return;
    }

   default :
    throw( std::invalid_argument(
		  "MCFLemonSolver::guts_of_poM: unknown MCFBlockSbstMod type" ) );

   }  // end( case )
  }  // end( MCFBlockSbstMod )

 // any remaining Modification is plainly ignored, since it must be an
 // "abstract" Modification, which this Solver does not need to look at

 }  // end( guts_of_poM )

/*--------------------------------------------------------------------------*/
/*------------------- End File MCFLemonSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
