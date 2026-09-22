/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of the MCFLemonSolver: the four LEMON algorithms, on the graph that
 * can change and on the static one, are attached to a small MCFBlock whose
 * optimal value is known at every step, and the MCFBlock is then changed in
 * all the ways a MCFBlock can change, i.e., costs, capacities, closing and
 * opening arcs, adding and removing arcs. After each change every Solver
 * that has followed the changes is compared with the known value, and so is
 * a Solver of each kind attached anew to the changed MCFBlock, which has to
 * build its graph with the arcs closed and deleted as they are.
 *
 * Besides, the cases that do not come up when an instance is simply solved
 * are checked on their own: a Solver destroyed without ever having been
 * attached to a MCFBlock, a MCFBlock that has no capacities and no costs to
 * begin with and is given them by a Modification, fractional costs and
 * flows, negative costs on arcs of infinite capacity, a capacity many orders
 * of magnitude larger than the other data, each method of CostScaling,
 * deficits that sum to zero only up to rounding and deficits that do not,
 * the parameters looked up by their name, the invalid values of the
 * parameters of the algorithms, and the DMX file of the instance.
 *
 * The exit code is the number of failed checks.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFLemonSolver.h"

#include <cmath>

#include <cstdio>

#include <fstream>

#include <iostream>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;  // number of failed checks

static const std::vector< std::string > ALGOS = {
 "NetworkSimplex" , "CycleCanceling" , "CostScaling" , "CapacityScaling" };

static const std::vector< std::string > GRAPHS = {
 "MCFListDigraph" , "SmartDigraph" };

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/

static std::string name( const std::string & algo , const std::string & gr )
{
 return( "MCFLemonSolver" + algo + "<" + gr + ",double,double>" );
 }

/*--------------------------------------------------------------------------*/
/// checks the value a Solver gives against the known one, +Inf = infeasible,
/// -Inf = unbounded

static void check_value( Solver * slv , double expected ,
			 const std::string & what )
{
 const auto who = what + ", " + slv->classname();
 const int status = slv->compute();

 if( expected == Inf< double >() ) {
  check( status == Solver::kInfeasible , who + ": infeasible expected" );
  check( slv->get_ub() == Inf< double >() , who + ": ub is not +Inf" );
  return;
  }

 if( expected == -Inf< double >() ) {
  check( status == Solver::kUnbounded , who + ": unbounded expected" );
  check( slv->get_ub() == -Inf< double >() , who + ": ub is not -Inf" );
  return;
  }

 check( status == Solver::kOK , who + ": status " + std::to_string( status ) );
 const double ub = slv->get_ub();
 check( std::abs( ub - expected ) <= 1e-9 * std::max( 1.0 ,
						   std::abs( expected ) ) ,
	who + ": value " + std::to_string( ub ) + " instead of " +
	std::to_string( expected ) );
 check( slv->get_lb() == ub , who + ": lb differs from ub" );
 }

/*--------------------------------------------------------------------------*/
/// every Solver attached to the MCFBlock, and one of each kind attached anew

static void check_all( MCFBlock & mcf , double expected ,
		       const std::string & what ,
		       const std::vector< std::string > & algos = ALGOS )
{
 for( auto slv : mcf.get_registered_solvers() )
  check_value( slv , expected , what + " (warm)" );

 for( const auto & gr : GRAPHS )
  for( const auto & algo : algos ) {
   auto slv = Solver::new_Solver( name( algo , gr ) );
   mcf.register_Solver( slv );
   check_value( slv , expected , what + " (new)" );
   mcf.unregister_Solver( slv , true );
   }
 }

/*--------------------------------------------------------------------------*/
/* The instance: 4 nodes, 4 units of flow from node 1 to node 4, and room
 * for two more arcs besides the five below (costs, capacities):
 *
 *   0: 1 -> 2  (1, 3)    2: 2 -> 3  (1, 2)    4: 3 -> 4  (1, 5)
 *   1: 1 -> 3  (4, 5)    3: 2 -> 4  (5, 3)
 *
 * The optimum sends 2 units along 1-2-3-4 and 2 along 1-3-4, i.e., 16.
 * load() generates the flow Variable, whose fixing closes an arc. */

static void load( MCFBlock & mcf , bool withdata ,
		  const MCFBlock::Vec_FNumber & b = { -4 , 0 , 0 , 4 } ,
		  const MCFBlock::Vec_CNumber & c = { 1 , 4 , 1 , 5 , 1 } ,
		  const MCFBlock::Vec_FNumber & u = { 3 , 5 , 2 , 3 , 5 } )
{
 const MCFBlock::Subset sn = { 1 , 1 , 2 , 2 , 3 };
 const MCFBlock::Subset en = { 2 , 3 , 3 , 4 , 4 };

 if( withdata )
  mcf.load( 4 , 5 , en , sn , u , c , b , 0 , 0 , 0 , 2 );
 else
  mcf.load( 4 , 5 , en , sn , {} , {} , b , 0 , 0 , 0 , 2 );

 }

/*--------------------------------------------------------------------------*/
/// the MCFBlock changed in every way, followed by the Solver on MCFListDigraph

static void test_changes( void )
{
 MCFBlock mcf;
 load( mcf , true );

 for( const auto & algo : ALGOS )
  mcf.register_Solver( Solver::new_Solver( name( algo , "MCFListDigraph" ) ) );

 check_all( mcf , 16 , "initial instance" );

 mcf.chg_cost( 10 , 2 );  // 2 -> 3 now costs 10: 4 units along 1-3-4
 check_all( mcf , 20 , "cost changed" );

 mcf.close_arc( 1 );      // 1 -> 3 closed: 1 -> 2 alone carries 3 < 4
 check_all( mcf , Inf< double >() , "arc closed" );

 // 1 -> 4 is arc 5: 3 units along 1-2-4 and 1 along 1-4
 check( mcf.add_arc( 1 , 4 , 7 , 10 ) == 5 , "name of the added arc" );
 check_all( mcf , 25 , "arc added" );

 mcf.open_arc( 1 );       // 1 -> 3 back: 4 units along 1-3-4
 check_all( mcf , 20 , "arc opened" );

 mcf.chg_ucap( 2 , 4 );   // 3 -> 4 carries 2: 2 along 1-3-4, 2 along 1-2-4
 check_all( mcf , 22 , "capacity changed" );

 mcf.remove_arc( 5 );     // 1 -> 4 was not used
 check_all( mcf , 22 , "arc removed" );

 // the new arc takes the name of the one removed, and all the flow
 check( mcf.add_arc( 1 , 4 , 1 , 10 ) == 5 , "name of the re-added arc" );
 check_all( mcf , 4 , "arc re-added" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/// no capacities and no costs to begin with, given by a Modification

static void test_empty_data( void )
{
 MCFBlock mcf;
 load( mcf , false );

 for( const auto & algo : ALGOS )
  mcf.register_Solver( Solver::new_Solver( name( algo , "MCFListDigraph" ) ) );

 check_all( mcf , 0 , "no costs, no capacities" );

 const MCFBlock::Vec_CNumber c = { 1 , 4 , 1 , 5 , 1 };
 mcf.chg_costs( c.begin() );  // uncapacitated: 4 units along 1-2-3-4
 check_all( mcf , 12 , "costs given" );

 mcf.chg_ucap( 2 , 2 );       // 2 along 1-2-3-4, 2 along 1-3-4
 check_all( mcf , 16 , "one capacity given" );

 check( mcf.add_arc( 1 , 4 , 1 ) == 5 , "name of the uncapacitated arc" );
 check_all( mcf , 4 , "uncapacitated arc added" );

 mcf.unregister_Solvers( true );
 }

/*--------------------------------------------------------------------------*/
/* The deficits sum to zero only up to rounding: in double, the supplies
 * 3.7, -0.1, -0.3, -3.3 sum to 3.6e-16 in any order, which LEMON alone
 * reads as infeasible. The optimum sends 0.1 to node 2 along 1-2, 0.3 to
 * node 3 along 1-2-3, and to node 4 1.7 along 1-2-3-4 and 1.6 along 1-3-4,
 * i.e., 0.1 + 0.6 + 5.1 + 8 = 13.8. Deficits that really do not sum to zero
 * make the problem infeasible, which LEMON alone does not say when the
 * demand exceeds the supply: it leaves part of the demand unmet. */

static void test_balance( void )
{
 MCFBlock mcf;
 load( mcf , true , { -3.7 , 0.1 , 0.3 , 3.3 } );
 check_all( mcf , 13.8 , "deficits summing to zero up to rounding" );

 MCFBlock mcf2;
 load( mcf2 , true , { -4 , 0 , 0 , 4.1 } );
 check_all( mcf2 , Inf< double >() , "demand exceeding supply" );
 }

/*--------------------------------------------------------------------------*/
/* The costs of the instance divided by 10, i.e., fractional and not exact in
 * binary: the optimum is the same flow, and its value 1.6. LEMON is exact
 * only on integer costs, which is why MCFLemonSolver compiles the network
 * simplex with tolerances and gives CostScaling and CycleCanceling the costs
 * scaled and rounded; every algorithm has to give the value of the true
 * costs. */

static void test_fractional_costs( void )
{
 MCFBlock mcf;
 load( mcf , true , { -4 , 0 , 0 , 4 } , { 0.1 , 0.4 , 0.1 , 0.5 , 0.1 } );
 check_all( mcf , 1.6 , "fractional costs" );

 mcf.chg_cost( 1.0 , 2 );  // 2 -> 3 now costs 1: 4 units along 1-3-4
 check_all( mcf , 2 , "fractional costs changed" );
 }

/*--------------------------------------------------------------------------*/
/* 8 units from node 1 to node 3 on three arcs of fractional costs (costs,
 * capacities):
 *
 *   0: 1 -> 2  (37.41, 3)    1: 2 -> 3  (22.43, 6)    2: 1 -> 3  (92.89, 10)
 *
 * The optimum sends 3 units along 1-2-3 and 5 along 1-3, i.e., 643.97.
 * CapacityScaling, run on these costs as they are, answers "infeasible". */

static void test_fractional_costs_paths( void )
{
 MCFBlock mcf;
 const MCFBlock::Subset sn = { 1 , 2 , 1 };
 const MCFBlock::Subset en = { 2 , 3 , 3 };
 const MCFBlock::Vec_CNumber c = { 37.41 , 22.43 , 92.89 };
 const MCFBlock::Vec_FNumber u = { 3 , 6 , 10 };
 const MCFBlock::Vec_FNumber d = { -8 , 0 , 8 };
 mcf.load( 3 , 3 , en , sn , u , c , d );
 check_all( mcf , 643.97 , "fractional costs along paths" );
 }

/*--------------------------------------------------------------------------*/
/* The capacities and the deficits of the instance divided by 10, i.e.,
 * fractional and not exact in binary: the optimum is the same flow divided
 * by 10, and its value 1.6. LEMON is exact only on integer capacities and
 * supplies, which is why MCFLemonSolver gives CostScaling, CycleCanceling and
 * CapacityScaling the flows scaled and rounded; every algorithm has to give
 * the value of the true flows. */

static void test_fractional_flows( void )
{
 MCFBlock mcf;
 load( mcf , true , { -0.4 , 0 , 0 , 0.4 } , { 1 , 4 , 1 , 5 , 1 } ,
       { 0.3 , 0.5 , 0.2 , 0.3 , 0.5 } );
 check_all( mcf , 1.6 , "fractional flows" );

 // 2 -> 4, which no optimal flow uses, gets a capacity many orders of
 // magnitude larger than the others, which must not set the rounding of
 // the others
 mcf.chg_ucap( 1e12 , 3 );
 check_all( mcf , 1.6 , "fractional flows, a large unused capacity" );

 // 2 -> 3 carries 0.05: 0.05 along 1-2-3-4 and 0.35 along 1-3-4
 mcf.chg_ucap( 0.05 , 2 );
 check_all( mcf , 1.9 , "fractional flows, capacity changed" );
 }

/*--------------------------------------------------------------------------*/
/* The fractional flows of test_fractional_flows(), with the arc of negative
 * cost 3 -> 2 and a capacity of 1e15 on it: its capacity then bounds the
 * flow of the optimal ones, and scaled so that it is integer the supplies
 * are rounded to multiples of 2, i.e., a flow of the scaled data is not one
 * of the true data. CostScaling, CycleCanceling and CapacityScaling have to
 * end in error rather than give it out, the network simplex, which is given
 * the data as they are, has to solve it. */

static void test_scale_out_of_range( void )
{
 MCFBlock mcf;
 load( mcf , true , { -0.4 , 0 , 0 , 0.4 } , { 1 , 4 , 1 , 5 , 1 } ,
       { 0.3 , 0.5 , 0.2 , 0.3 , 0.5 } );
 check( mcf.add_arc( 3 , 2 , -0.5 , 1e15 ) == 5 , "name of the large arc" );

 for( const auto & algo : ALGOS ) {
  auto slv = Solver::new_Solver( name( algo , "MCFListDigraph" ) );
  mcf.register_Solver( slv );
  if( algo == "NetworkSimplex" )
   check_value( slv , 1.6 , "scale out of range" );
  else {
   const int status = slv->compute();
   check( status == Solver::kError , "scale out of range, " +
	  slv->classname() + ": status " + std::to_string( status ) );
   check( ! slv->has_var_solution() , "scale out of range, " +
	  slv->classname() + ": a solution is given" );
   }
  mcf.unregister_Solver( slv , true );
  }
 }

/*--------------------------------------------------------------------------*/
/* No capacities at all, and an arc of negative cost: 3 -> 2 at -0.5 makes
 * with 2 -> 3 a cycle of positive cost, so that the optimum is finite and
 * the one without it, 4 units along 1-2-3-4, i.e., 12. CostScaling,
 * CycleCanceling and CapacityScaling saturate an arc of negative cost and
 * answer "unbounded" when its capacity is infinite, which is why
 * MCFLemonSolver gives them a finite bound in its place. At -2 the cycle is
 * of negative cost and infinite capacity, and the problem is unbounded. */

static void test_negative_costs( void )
{
 MCFBlock mcf;
 load( mcf , false );
 const MCFBlock::Vec_CNumber c = { 1 , 4 , 1 , 5 , 1 };
 mcf.chg_costs( c.begin() );
 check( mcf.add_arc( 3 , 2 , -0.5 ) == 5 , "name of the negative arc" );
 check_all( mcf , 12 , "negative cost, finite optimum" );

 mcf.chg_cost( -2 , 5 );
 check_all( mcf , -Inf< double >() , "negative cycle of infinite capacity" );
 }

/*--------------------------------------------------------------------------*/
/// a Solver that never saw a MCFBlock is destroyed, of every kind

static void test_unattached( void )
{
 for( const auto & gr : GRAPHS )
  for( const auto & algo : ALGOS ) {
   auto slv = Solver::new_Solver( name( algo , gr ) );
   check( slv != nullptr , name( algo , gr ) + " is not in the factory" );
   delete slv;
   }
 }

/*--------------------------------------------------------------------------*/
/// the parameters, by name and by value

static void test_parameters( void )
{
 for( const auto & algo : ALGOS ) {
  auto slv = Solver::new_Solver( name( algo , "MCFListDigraph" ) );
  const auto who = slv->classname();

  const auto s = slv->str_par_str2idx( "strDMXFile" );
  check( slv->str_par_idx2str( s ) == "strDMXFile" ,
	 who + ": strDMXFile by name" );

  const auto d = slv->dbl_par_str2idx( "dblMaxTime" );
  check( slv->dbl_par_idx2str( d ) == "dblMaxTime" ,
	 who + ": dblMaxTime by name" );

  // CycleCanceling has three methods, CostScaling three plus the automatic
  // choice, which is its default
  if( algo == "CycleCanceling" || algo == "CostScaling" ) {
   const auto k = slv->int_par_str2idx( "kMethod" );
   const int last = algo == "CostScaling" ? 3 : 2;
   if( algo == "CostScaling" )
    check( slv->get_dflt_int_par( k ) == 3 , who + ": kMethod default" );
   for( int m : { -1 , last + 1 } ) {
    bool thrown = false;
    try { slv->set_par( k , m ); }
    catch( std::invalid_argument & ) { thrown = true; }
    check( thrown , who + ": kMethod " + std::to_string( m ) + " accepted" );
    }
   slv->set_par( k , last );
   check( slv->get_int_par( k ) == last , who + ": kMethod not set" );
   }

  delete slv;
  }
 }

/*--------------------------------------------------------------------------*/
/* Each method of CostScaling, and the automatic choice among them, on the
 * integer instance, whose flows are given as they are, on the fractional
 * flows, which are scaled, and on the negative costs of infinite
 * capacity. */

static void test_cost_scaling_methods( void )
{
 MCFBlock mcf;
 load( mcf , true );
 MCFBlock frac;
 load( frac , true , { -0.4 , 0 , 0 , 0.4 } , { 1 , 4 , 1 , 5 , 1 } ,
       { 0.3 , 0.5 , 0.2 , 0.3 , 0.5 } );
 MCFBlock neg;
 load( neg , false );
 const MCFBlock::Vec_CNumber c = { 1 , 4 , 1 , 5 , 1 };
 neg.chg_costs( c.begin() );
 neg.add_arc( 3 , 2 , -0.5 );

 for( int m = 0 ; m <= 3 ; ++m ) {
  const auto what = "CostScaling kMethod " + std::to_string( m );
  for( auto [ blk , value ] : { std::pair( &mcf , 16.0 ) ,
				std::pair( &frac , 1.6 ) ,
				std::pair( &neg , 12.0 ) } ) {
   auto slv = Solver::new_Solver( name( "CostScaling" , "MCFListDigraph" ) );
   slv->set_par( slv->int_par_str2idx( "kMethod" ) , m );
   blk->register_Solver( slv );
   check_value( slv , value , what );
   blk->unregister_Solver( slv , true );
   }
  }
 }

/*--------------------------------------------------------------------------*/
/// the DMX file is the instance, closed and added arcs included

static void test_dmx( void )
{
 const std::string file = "MCFLemonSolver_test.dmx";

 MCFBlock mcf;
 load( mcf , true );
 mcf.close_arc( 1 );
 mcf.add_arc( 1 , 4 , 7 , 10 );

 for( const auto & gr : GRAPHS ) {
  auto slv = Solver::new_Solver( name( "NetworkSimplex" , gr ) );
  mcf.register_Solver( slv );
  slv->set_par( slv->str_par_str2idx( "strDMXFile" ) , file );
  slv->compute();

  std::ifstream in( file );
  std::string head;
  std::getline( in , head );
  check( head == "p min 4 6" , gr + ": DMX header \"" + head + "\"" );
  int arcs = 0;
  for( std::string line ; std::getline( in , line ) ; )
   arcs += ( line.rfind( "a" , 0 ) == 0 );
  check( arcs == 6 , gr + ": " + std::to_string( arcs ) + " arcs in the DMX" );

  mcf.unregister_Solver( slv , true );
  std::remove( file.c_str() );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 test_unattached();
 test_parameters();
 test_changes();
 test_empty_data();
 test_balance();
 test_fractional_costs();
 test_fractional_costs_paths();
 test_fractional_flows();
 test_negative_costs();
 test_scale_out_of_range();
 test_cost_scaling_methods();
 test_dmx();

 if( failures )
  std::cout << failures << " checks failed" << std::endl;
 else
  std::cout << "All tests passed!!" << std::endl;

 return( failures );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
