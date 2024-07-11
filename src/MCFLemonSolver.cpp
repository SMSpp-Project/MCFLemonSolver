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

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
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
 #define SMSpp_which_insert_LEMON 7
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
// - three data tyopes ( int , long , double ) in two places ( costs , flows ),
//   i.e., 3^2 = 9 possible data configurations
// so it is 4 * 2 * 8 = 64 versione

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
 delete f_algo; // delete previous f_algo (if present)
 f_algo = nullptr;

 // create and clear new Graph (MCFListDigraph or SmartDigraph)
 dgp = new GR;
 dgp->clear();

 // actually reserve get_MaxNNodes() node space in dgp
 auto maxNNodes = MCFB->get_MaxNNodes();
 dgp->reserveNode( maxNNodes );
 MCFBlock::Index n = MCFB->get_NNodes();
 // add all the nodes, one by one (looks stupid you con't do all in one blow,
 // but who are we to say ...)
 for( MCFBlock::Index i = 0 ; i < n ; ++i )
  dgp->addNode();

 // reserve get_MaxNArcs() arcs space
 auto maxNArcs = MCFB->get_MaxNArcs();
 dgp->reserveArc( maxNArcs );
 MCFBlock::Index m = MCFB->get_NArcs();
 n_arcs = m;
 n_arcs_added = 0;
 // get starting node subset and ending node subset
 MCFBlock::c_Subset& sn = MCFB->get_SN();
 MCFBlock::c_Subset& en = MCFB->get_EN();

 // add all arcs in dgp
 // sn[ i ] - 1 and en[ i ] - 1 are used because sn , en nodes start from 1
 for( MCFBlock::Index i = 0 ; i < m ; ++i )
  if constexpr( std::is_same< GR , SmartDigraph >::value )
   dgp->addArc( dgp->nodeFromId( sn[ i ] - 1 ) ,
                dgp->nodeFromId( en[ i ] - 1 ) );
  else
   dgp->addArc( sn[ i ] - 1 , en[ i ] - 1 );

 // new instance of Lemon Algorithm with no arc mixing
 if constexpr( std::is_same< Algo< GR , V , C > ,
	                     NetworkSimplex< GR , V , C > >::value )
  f_algo = new Algo< GR , V , C >( *dgp , false );
 else
  f_algo = new Algo< GR , V , C >( *dgp );

 // defining names for types for readability
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;

 // now we are going to fill up ArcMap and NodeMap
 if( ! MCFB->get_U().empty() ) {  // this is the case of upperMap
  um = new MCFArcMapV( *dgp );    // create the upper map
  auto& u = MCFB->get_U();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   um->set( dgp->arcFromId( i ) , u[ i ] );

  f_algo->upperMap( *um );        // pass it to the Algo
  }

 if( ! MCFB->get_C().empty() ) {  // this is the case of CostMap
  cm = new MCFArcMapV( *dgp );    // create the cost map
  auto& c = MCFB->get_C();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   cm->set( dgp->arcFromId( i ) , c[ i ] );

  f_algo->costMap( *cm );  // pass it to the Algo
  }

 if( ! MCFB->get_B().empty() ) {  // this is the case of supplyMap
  bm = new MCFNodeMapV( *dgp );   // create the supply map
  auto& b = MCFB->get_B();
  for( MCFBlock::Index i = 0 ; i < n ; i++ )
   bm->set( dgp->nodeFromId( i ) ,
            b[ i ] == 0 ? 0 : -b[ i ] );  // note that supply = - deficit

  f_algo->supplyMap( *bm );       // pass it to the algo
  }
 }  // end( guts_of_set_Block )

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
                         "MCFSolver:set_Block: block must be a MCFBlock" ) );

 bool owned = MCFB->is_owned_by( f_id );
 if( ( ! owned ) && ( ! MCFB->read_lock() ) )
  throw( std::logic_error( "cannot acquire read_lock on MCFBlock" ) );

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

 if( ! f_Block ) // there is no [MCFBlock] to solve
  return( kBlockLocked ); // return error

 bool owned = f_Block->is_owned_by( f_id ); // check if already locked
 if( ( ! owned ) && ( ! f_Block->read_lock() ) ) // if not try to read_lock
  return( kBlockLocked ); // return error on failure

 // while [read_]locked, process any outstanding Modification
 process_outstanding_Modification();

 if( ! f_dmx_file.empty() ) {  // if so required
  // output the current instance (after the changes) to a DMX file
  if constexpr( std::is_same< GR , SmartDigraph >::value ) {
   std::ofstream ProbFile( f_dmx_file , ios_base::out | ios_base::trunc );
   if( ! ProbFile.is_open() )
    throw( std::logic_error( "cannot open DMX file " + f_dmx_file ) );
   writeDimacsMat( ProbFile , *dgp );
   ProbFile.close();
   }
  else
   throw( std::logic_error( "saving to DMX file not supported on ListDigraph"
			    ) );
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
  f_algo->costMap( *cm );
  cost_changed = false;
  }

 if( cap_changed ) {
  f_algo->upperMap( *um );
  cap_changed = false;
  }

 if( supply_changed ) {
  f_algo->supplyMap( *bm );
  supply_changed = false;
  }

 // using Lemon Digraph Writer for debuggging
 // std::ofstream ProbFile("lmn.txt", ios_base::out | ios_base::trunc);
 // if (!ProbFile.is_open()) {
 //     throw(std::logic_error("cannot open file lmn.txt"));
 // }
 // DigraphWriter<GR>(*dgp, ProbFile).nodeMap("supply", *bm).arcMap("cost",
 // *cm).arcMap("upper", *um).run(); ProbFile.close();

 auto start = chrono::system_clock::now();
 guts_of_compute(); // here the actual magic is done by specialised classes
 auto end = chrono::system_clock::now();

 chrono::duration< double > elapsed = end - start;
 ticks = elapsed.count();

 unlock(); // release self-lock

 // now give out the result
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
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is, and call the appropriate
  * method of LEMON. */

 // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto tmod = dynamic_cast< const GroupModification* >( mod ) ) {
  for( const auto& submod : tmod->sub_Modifications() )
   guts_of_poM( submod.get() );

  return;
  }

 // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Note: in the following we can assume that C, B and U are nonempty. This
  * is because they can be empty only if they are so when the object is
  * loaded. But if a Modification has been issued they are no longer empty (a
  * Modification changin nothing from the "empty" state is not issued). */

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
     throw( std::logic_error(
                      "SmartDigraph doesn't support operations on arc" ) );
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
     throw( std::logic_error(
                      "SmartDigraph doesn't support operations on arc" ) );
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

     dgp->addArc( startNode , endNode );
     auto arc = dgp->arcFromId( rng.first );

     n_arcs_added++;
     cm->set( arc , std::isnan( ca ) ? 0 : ca );
     um->set( arc , capacity );

     cost_changed = true;
     cap_changed = true;
     }
    else
     throw( std::logic_error(
                    "SmartDigraph doesn't support operations on arc" ) );
    return;
    }

   // removing arcs- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   case( MCFBlockMod::eRmvArc ) : {
    // only if ListDigraph is used
    if constexpr( std::is_same< GR , MCFListDigraph >::value ) {
     // if the arc was previous deleted, return
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
     throw( std::logic_error(
		    "SmartDigraph doesn't support operations on arc" ) );
    return;
    }

   default :
    throw( std::invalid_argument( "unknown MCFBlockRngdMod type" ) );

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
     if( ! std::isnan( CC[ i ] ) && dgp->valid( dgp->arcFromId( i ) ) )
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
     throw( std::logic_error(
		      "SmartDigraph doesn't support operations on arc" ) );
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
     throw( std::logic_error(
		     "SmartDigraph doesn't support operations on arc" ) );
    return;
    }

   default :
    throw( std::invalid_argument( "unknown MCFBlockSbstMod type" ) );

   }  // end( case )
  }  // end( MCFBlockSbstMod )

 // any remaining Modification is plainly ignored, since it must be an
 // "abstract" Modification, which this Solver does not need to look at

 }  // end( guts_of_poM )

/*--------------------------------------------------------------------------*/
/*------------------- End File MCFLemonSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
