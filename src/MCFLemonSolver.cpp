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
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MCFLemonSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register the various LEMONSolver< Alg , GR , V , C > to the Solver factory

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverNetworkSimplex< SmartDigraph , double , double > );

/*!!
SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCycleCanceling< SmartDigraph, int, int> );
	    !!*/

/*!!
SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< SmartDigraph , double , double > );
	    !!*/

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , double , double > );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of MCFLemonSolver ------------------------*/
/*--------------------------------------------------------------------------*/

template< template< typename , typename , typename > class Algo ,
	  typename GR , typename V , typename C >
 requires LEMONGraph< GR >
void MCFLemonSolver< Algo , GR , V , C >::set_Block( Block * block )
{
 if( block == f_Block )  // actually doing nothing
  return;                // cowardly and silently return

 delete f_algo;
 f_algo = nullptr;

 Solver::set_Block( block );  // attach to the new Block

 if( ! block )  // this is just: go sit down in a corner and wait
  return;       // all done

 auto MCFB = dynamic_cast< MCFBlock * >( block );
 if( ! MCFB )
  throw( std::invalid_argument(
                          "MCFSolver:set_Block: block must be a MCFBlock") );

 bool owned = MCFB->is_owned_by(f_id);
 if( ( ! owned ) && ( ! MCFB->read_lock() ) )
  throw( std::logic_error( "cannot acquire read_lock on MCFBlock" ) );

 // create and clear new Graph (ListDigraph or SmartDigraph)
 dgp = new GR;
 dgp->clear();

 // actually reserve get_MaxNNodes() node space in dgp
 dgp->reserveNode( MCFB->get_MaxNNodes() );
 MCFBlock::Index n = MCFB->get_NNodes();

 // add all the nodes, one by one (looks stupid you con't do all in one blow,
 // but who are we to say ...)
 for( MCFBlock::Index i = 0 ; i < n ; ++i )
  dgp->addNode();

 // reserve get_MaxNArcs() arcs space
 dgp->reserveArc( MCFB->get_MaxNArcs() );
 MCFBlock::Index m = MCFB->get_NArcs();

 // get starting node subset and ending node subset
 MCFBlock::c_Subset & sn = MCFB->get_SN();
 MCFBlock::c_Subset & en = MCFB->get_EN();

 // add all arcs in dgp
 // sn[ i ] - 1 and en[ i ] - 1 are used because sn , en nodes start from 1
 for( MCFBlock::Index i = 0 ; i < m ; ++i )
  dgp->addArc( dgp->nodeFromId( sn[ i ] - 1 ) ,
	       dgp->nodeFromId( en[ i ] - 1 ) );

 // new instance of Lemon Algorithm
 f_algo = new Algo< GR , V, C >( * dgp );

 // defining names for types for readability
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;

 // now we are going to fill up ArcMap and NodeMap
 // nhis is the case of upperMap 
 if( ! MCFB->get_U().empty() ) {
  MCFArcMapV um( * dgp );  // create the upper map
  auto & u = MCFB->get_U();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   um.set( dgp->arcFromId( i ) , u[ i ] );

  f_algo->upperMap( um );  // pass it to the Algo
  }

 // this is the case of CostMap
 if( ! MCFB->get_C().empty() ) {
  MCFArcMapV cm( * dgp );  // create the cost map
  auto & c = MCFB->get_C();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   cm.set( dgp->arcFromId( i ) , c[ i ] );

  f_algo->costMap( cm );  // pass it to the Algo
  }

 // this is the case of supplyMap
 if( ! MCFB->get_B().empty() ) {
  MCFNodeMapV bm( * dgp );  // create the supply map
  auto & b = MCFB->get_B();
  for( MCFBlock::Index i = 0 ; i < n ; i++ )
   bm.set( dgp->nodeFromId( i ) , -b[ i ] );  // note tat supply = - deficit

  f_algo->supplyMap( bm );  // pass it to the algo
  }

 if ( ! owned )
  MCFB->read_unlock();
      
  }  // end( MCFLemonSolver< Algo , GR , V , C >set_Block )

/*--------------------------------------------------------------------------*/

template< template< typename , typename , typename > class Algo ,
	  typename GR , typename V , typename C >
 requires LEMONGraph< GR >
int MCFLemonSolver< Algo , GR , V , C >::compute( bool changedvars )
{
 const static std::array< int , 6 > LemonStatus_2_MCFstatus = {
  kErrorStatus, LEMON_sol_type::OPTIMAL, kErrorStatus , 
  LEMON_sol_type::INFEASIBLE,
  LEMON_sol_type::UNBOUNDED, kErrorStatus };

 const static std::array< int , 6 > MCFstatus_2_sol_type = {
  kUnEval , Solver::kOK , kStopTime , kInfeasible , Solver::kUnbounded ,
  Solver::kError };

 lock(); // first of all, acquire self-lock

 if ( ! f_Block )          // there is no [MCFBlock] to solve
  return( kBlockLocked );  // return error

 bool owned = f_Block->is_owned_by( f_id );       // check if already locked
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
  return( kBlockLocked );                         // return error on failure

 // while [read_]locked, process any outstanding Modification
 //TODO: ensure that modification are actually processed for MCFLemonSolver.
 //process_outstanding_Modification();
        
 if( ! f_dmx_file.empty() )  {  // if so required
  // output the current instance (after the changes) to a DMX file
  std::ofstream ProbFile( f_dmx_file , ios_base::out | ios_base::trunc );
  if( ! ProbFile.is_open() )
   throw( std::logic_error( "cannot open DMX file " + f_dmx_file ) );

  writeDimacsMat( ProbFile , *dgp );
  ProbFile.close();
  }

 if( ! owned )             // if the [MCF]Block was actually read_locked
  f_Block->read_unlock();  // read_unlock it

 auto start = chrono::system_clock::now();

 guts_of_compute();  // here the actual magic is done by specialised classes

 auto end = chrono::system_clock::now();

 chrono::duration< double > elapsed = end - start;
 ticks = elapsed.count();
        
 unlock();  // release self-lock

 // now give out the result: note that the vector MCFstatus_2_sol_type[]
 // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
 // hence the returned status has to be shifted by + 1
 return( MCFstatus_2_sol_type[ this->get_status() ] );

 }  // end( MCFLemonSolver< Algo , GR , V , C >::compute )

/*--------------------------------------------------------------------------*/
/*------------------- End File MCFLemonSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
 
