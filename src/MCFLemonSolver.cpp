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

// register the various LEMONSolver GR , V , C > to the Solver factory

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverNetworkSimplex< SmartDigraph , int , int > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverNetworkSimplex< SmartDigraph , double , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverNetworkSimplex< SmartDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverNetworkSimplex< ListDigraph , int , int > );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverNetworkSimplex< ListDigraph , double , double > );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverNetworkSimplex< ListDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCycleCanceling< SmartDigraph, int, int> );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCycleCanceling< SmartDigraph, double, double> );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCycleCanceling< SmartDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverCycleCanceling< ListDigraph, int, int> );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverCycleCanceling< ListDigraph, double, double> );

SMSpp_insert_in_factory_cpp_0_t(
             MCFLemonSolverCycleCanceling< ListDigraph, long, long> );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< SmartDigraph , int , int > );  

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< SmartDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< SmartDigraph , double , double > );   

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< ListDigraph , int , int > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< ListDigraph , double , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCostScaling< ListDigraph , long , long > );
        
SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , int , int > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , double , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , int , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< SmartDigraph , long , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< ListDigraph , int , int > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< ListDigraph , double , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< ListDigraph , long , long > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< ListDigraph , int , double > );

SMSpp_insert_in_factory_cpp_0_t(
            MCFLemonSolverCapacityScaling< ListDigraph , long , double > );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of MCFLemonSolver ------------------------*/
/*--------------------------------------------------------------------------*/

template<template< typename , typename , typename > class Algo ,
   LEMONGraph GR , typename V , typename C >
void MCFLemonSolver<Algo, GR, V, C>::guts_of_set_Block(MCFBlock *MCFB){

 //delete previous graph and f_algo (if present)
 //delete dgp;
 delete f_algo;
 f_algo = nullptr;
 // create and clear new Graph (ListDigraph or SmartDigraph)
if constexpr (std::is_same<GR, ListDigraph>::value){
  dgp = new MCFListDigraph();
}else{
  dgp = new GR;
}
 dgp->clear();



//  auto MCFB = static_cast<MCFBlock *>(f_block);
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
 f_algo = new Algo< GR , V, C >(  *dgp );
  
 // Algo<GR, V, C> * algo = f_algo->reset();
 // defining names for types for readability
 using MCFArcMapV = typename GR::template ArcMap< V >;
 using MCFNodeMapV = typename GR::template NodeMap< V >;

 // now we are going to fill up ArcMap and NodeMap
 // this is the case of upperMap 
 if( ! MCFB->get_U().empty() ) {
  um = new MCFArcMapV( * dgp );  // create the upper map
  auto & u = MCFB->get_U();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   um->set( dgp->arcFromId( i ) , u[ i ] );

  f_algo->upperMap( *um );  // pass it to the Algo
  }

 // this is the case of CostMap
 if( ! MCFB->get_C().empty() ) {
  cm = new MCFArcMapV( * dgp );  // create the cost map
  auto & c = MCFB->get_C();
  for( MCFBlock::Index i = 0 ; i < m ; ++i )
   cm->set( dgp->arcFromId( i ) , c[ i ] );

  f_algo->costMap( *cm );  // pass it to the Algo
  }

 // this is the case of supplyMap
 if( ! MCFB->get_B().empty() ) {
  bm = new MCFNodeMapV( * dgp );  // create the supply map
  auto & b = MCFB->get_B();
  for( MCFBlock::Index i = 0 ; i < n ; i++ )
   bm->set( dgp->nodeFromId( i ) , -b[ i ] );  // note tat supply = - deficit

  f_algo->supplyMap( *bm );  // pass it to the algo
  }


}

template< template< typename , typename , typename > class Algo ,
	  LEMONGraph GR , typename V , typename C >
void MCFLemonSolver< Algo , GR , V , C >::set_Block( Block * block )
{
 if( block == f_Block )  // actually doing nothing
  return;                // cowardly and silently return

 

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

  guts_of_set_Block( MCFB );  // fill up the graph and f_algo

 
 if ( ! owned )
  MCFB->read_unlock();
      
  }  // end( MCFLemonSolver< Algo , GR , V , C >set_Block )

/*--------------------------------------------------------------------------*/

template< template< typename , typename , typename > class Algo ,
	  LEMONGraph GR , typename V , typename C >
int MCFLemonSolver< Algo , GR , V , C >::compute( bool changedvars )
{
 
 const static std::array<int, 3> LEMONstatus_2_sol_type = {
  Solver::kInfeasible, Solver::kOK, Solver::kUnbounded };
 

 lock(); // first of all, acquire self-lock

 if ( ! f_Block )          // there is no [MCFBlock] to solve
  return( kBlockLocked );  // return error

 bool owned = f_Block->is_owned_by( f_id );       // check if already locked
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
  return( kBlockLocked );                         // return error on failure

 // while [read_]locked, process any outstanding Modification
 //TODO: ensure that modification are actually processed for MCFLemonSolver.
 process_outstanding_Modification();
        
 if( ! f_dmx_file.empty() )  {  // if so required
  // output the current instance (after the changes) to a DMX file
  std::ofstream ProbFile( f_dmx_file , ios_base::out | ios_base::trunc );
  if( ! ProbFile.is_open() )
   throw( std::logic_error( "cannot open DMX file " + f_dmx_file ) );
  if constexpr (std::is_same<GR, SmartDigraph>::value){
      writeDimacsMat(ProbFile, *dgp);
  }
  ProbFile.close();
  }else{
   // throw(std::logic_error("ListDigraph doesn't support writeDimacsMat function"));
  }

 if( ! owned )             // if the [MCF]Block was actually read_locked
  f_Block->read_unlock();  // read_unlock it

if( cost_changed ){
  f_algo->costMap(*cm);
  cost_changed = false;
}

if( cap_changed ){
  f_algo->upperMap(*um);
  cap_changed = false;
}

if( supply_changed ){
  f_algo->supplyMap(*bm);
  supply_changed = false;
}


 auto start = chrono::system_clock::now();

 guts_of_compute();  // here the actual magic is done by specialised classes

 auto end = chrono::system_clock::now();

 chrono::duration< double > elapsed = end - start;
 ticks = elapsed.count();
        
 unlock();  // release self-lock

 // now give out the result
 return( LEMONstatus_2_sol_type[ this->get_status() ] );
 //return Solver::kOK;

 }  // end( MCFLemonSolver< Algo , GR , V , C >::compute )

template<template<typename,typename,typename> class Algo, LEMONGraph GR, typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::add_Modification(sp_Mod &mod)
{

  if(std::is_same<GR, ListDigraph>::value){
  

    if (std::dynamic_pointer_cast<const NBModification>(mod))
    {
      // this is the "nuclear option": the MCFBlock has been re-loaded, so
      // the MCFClass solver also has to (immediately)
      auto MCFB = static_cast<MCFBlock *>(f_Block);
      guts_of_set_Block(MCFB);
      // besides, any outstanding modification makes no sense any longer
      mod_clear();
      }
      else
      push_back(mod);

  }else{
    throw (std::logic_error("SmartDigraph doesn't support Modification"));
  }
}

  
  template <template<typename,typename,typename> class Algo, LEMONGraph GR, typename V, typename C>
  void MCFLemonSolver<Algo, GR, V, C>::process_outstanding_Modification(void)
  {

    if(std::is_same<GR, SmartDigraph>::value){
    throw (std::logic_error("SmartDigraph doesn't support Modification"));
    }
    // no-frills loop: do them in order, with no attempt at optimizing
    // note that NBModification have already been dealt with and therefore need
    // not be considered here

    for (;;)
    {
      auto mod = pop();
      if (!mod)
        break;

      guts_of_poM(mod.get());
    }
  } // end( MCFLemonSolver::process_outstanding_Modification )

  /*--------------------------------------------------------------------------*/

  template <template<typename, typename, typename> class Algo, LEMONGraph GR, typename V, typename C>
  void MCFLemonSolver<Algo, GR, V, C>::guts_of_poM(c_p_Mod mod)
  {

    if constexpr (std::is_same<GR, ListDigraph>::value){
    

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      
      // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      /* This requires to patiently sift through the possible Modification types
      * to find what this Modification exactly is, and call the appropriate
      * method of MCFClass. */

      // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (auto tmod = dynamic_cast<const GroupModification *>(mod))
      {
        for (const auto &submod : tmod->sub_Modifications())
          guts_of_poM(submod.get());

        return;
      }

      // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      /* Note: in the following we can assume that C, B and U are nonempty. This
      * is because they can be empty only if they are so when the object is
      * loaded. But if a Modification has been issued they are no longer empty (a
      * Modification changin nothing from the "empty" state is not issued). */

      if (auto tmod = dynamic_cast<const MCFBlockRngdMod *>(mod))
      {
        auto rng = tmod->rng();

        switch (tmod->type())
        {
        case (MCFBlockMod::eChgCost):
          if (rng.second == rng.first + 1)
          {
            
            //Questo controllo può essere sostituito dalla funzione valid di ListDigraph
            if (dgp->valid(dgp->arcFromId(rng.first))){
               //TODO: change MCFC function to Algo function.
              //If arc is not deleted, then change cost
              //Create reference to the ArcMap that represents the cost and change its value
              //MCFC::ChgCost(rng.first, MCFB->get_C(rng.first));
              cm->set( dgp->arcFromId(rng.first), MCFB->get_C(rng.first));
              cost_changed = true;
            }
             
              
          }
          else
          {
            for(MCFBlock::Index i = rng.first; i < rng.second; i++){
                  cm->set(dgp->arcFromId(i), MCFB->get_C(i));
            
              }
            cost_changed = true;
          }

          
          return;

        case (MCFBlockMod::eChgCaps):
          if (rng.second == rng.first + 1)
          {
            if (dgp->valid(dgp->arcFromId(rng.first)))
              um->set(dgp->arcFromId(rng.first), MCFB->get_U(rng.first));
              cap_changed = true;
          }
          else{
            //MCFC::ChgUCaps(MCFB->get_U().data() + rng.first, nullptr,
            //              rng.first, rng.second);
            for(MCFBlock::Index i = rng.first; i < rng.second; i++){
                  um->set(dgp->arcFromId(i), MCFB->get_U(i));
              }
            cap_changed = true;
          }

          return;

        case (MCFBlockMod::eChgDfct):
          if (rng.second == rng.first + 1){
            //MCFC::ChgDfct(rng.first, MCFB->get_B(rng.first));
            bm->set(dgp->nodeFromId(rng.first), -MCFB->get_B(rng.first));
            supply_changed = true;
          }

          else{
            //MCFC::ChgDfcts(MCFB->get_B().data() + rng.first, nullptr,
                          // rng.first, rng.second);

            for(MCFBlock::Index i = rng.first; i < rng.second; i++){
                  bm->set(dgp->nodeFromId(i), -MCFB->get_B(i));
              }
            supply_changed = true;
          }
          return;
        //TODO: change MCFC function to Algo function.
        /*case (MCFBlockMod::eOpenArc):
          for (; rng.first < rng.second; ++rng.first)
            if ((!MCFB->is_deleted(rng.first)) &&
                (!MCFC::IsDeletedArc(rng.first)))
              MCFC::OpenArc(rng.first);
          return;
        //TODO: change MCFC function to Algo function.
        case (MCFBlockMod::eCloseArc):
          for (; rng.first < rng.second; ++rng.first)
            if ((!MCFB->is_deleted(rng.first)) &&
                (!MCFC::IsDeletedArc(rng.first)))
              MCFC::CloseArc(rng.first);
          return;
        */
        case (MCFBlockMod::eAddArc):
        {
          auto ca = MCFB->get_C(rng.first);
          // auto arc = MCFC::AddArc(MCFB->get_SN(rng.first),
          //                         MCFB->get_EN(rng.first),
          //                         MCFB->get_U(rng.first),
          //                         std::isnan(ca) ? 0 : ca);
          if(static_cast<MCFListDigraph*>(dgp)->first_free_arc == static_cast<int>(rng.first)){
            while(!first_free_arcs->empty()) first_free_arcs->pop();
          }else if(!first_free_arcs->empty()){
            static_cast<MCFListDigraph*>(dgp)->first_free_arc = first_free_arcs->front();
            first_free_arcs->pop();
          }
          auto arc = dgp->addArc(dgp->addNode(), dgp->addNode());
          cm->set(arc , ca);
          um->set(arc, MCFB->get_U(rng.first));
          cost_changed = true;
          cap_changed = true;
          if (arc != dgp->arcFromId(rng.first))
            throw(std::logic_error("name mismatch in AddArc()"));
          return;
        }

        case (MCFBlockMod::eRmvArc):
          //TODO: change MCFC function to Algo function.
          // MCFC::DelArc(rng.second - 1);
          dgp->erase(dgp->arcFromId(rng.second - 1));
          first_free_arcs->push(rng.second - 1);
          
          return;

        default:
          throw(std::invalid_argument("unknown MCFBlockRngdMod type"));
        }
      }

      // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (auto tmod = dynamic_cast<const MCFBlockSbstMod *>(mod))
      {/*
        switch (tmod->type())
        {
        case (MCFBlockMod::eOpenArc):
          for (auto arc : tmod->nms())
            if ((!MCFB->is_deleted(arc)) &&
                (!MCFC::IsDeletedArc(arc)))
              //TODO: change MCFC function to Algo function.
              MCFC::OpenArc(arc);
          return;

        case (MCFBlockMod::eCloseArc):
          for (auto arc : tmod->nms())
            if ((!MCFB->is_deleted(arc)) &&
                (!MCFC::IsDeletedArc(arc)))
              //TODO: change MCFC function to Algo function.
              MCFC::CloseArc(arc);
          return;
        }*/

        // have to InINF-terminate the vector of indices (damn!)
        // meanwhile, when appropriate remove the indices of deleted
        // arcs for which the operations make no sense;
        ;
        switch (tmod->type())
        {
        case (MCFBlockMod::eChgCost):
        {
          //Check if the subset of arcs is valid, then change the costs
          //on the ArcMap that represents the costs
          // MCFBlock::Subset nmsI;
          // nmsI.reserve(tmod->nms().size() + 1);
          // MCFBlock::Vec_CNumber NCost;
          // NCost.reserve(tmod->nms().size());
          for (auto i : tmod->nms())
            if (dgp->valid(dgp->arcFromId(i)))
            {
              auto arc = dgp->arcFromId(i);
              // NCost.push_back(ci);
              // nmsI.push_back(i);
              cm->set(arc, MCFB->get_C(i));

            }
          cost_changed = true;

            //f_algo->costMap(*cm);
          //nmsI.push_back(Inf<MCFBlock::Index>());
          //MCFC::ChgCosts(NCost.data(), nmsI.data());
          return;
        }

        case (MCFBlockMod::eChgCaps):
        {
          // MCFBlock::Subset nmsI;
          // nmsI.reserve(tmod->nms().size() + 1);
          // MCFBlock::Vec_FNumber NCap;
          // NCap.reserve(tmod->nms().size());
          auto &CC = MCFB->get_C();
          auto &U = MCFB->get_U();
          for (auto i : tmod->nms())
            if (!std::isnan(CC[i]) && dgp->valid(dgp->arcFromId(i)))
            {
              // NCap.push_back(U[i]);
              // nmsI.push_back(i);
              um->set(dgp->arcFromId(i), U[i]);

            }
          cap_changed = true;
          //nmsI.push_back(Inf<MCFBlock::Index>());
          //MCFC::ChgUCaps(NCap.data(), nmsI.data());
          return;
        }

        case (MCFBlockMod::eChgDfct):
        {
          MCFBlock::Vec_FNumber NDfct(tmod->nms().size());
          MCFBlock::Subset nmsI(tmod->nms().size() + 1);
          *copy(tmod->nms().begin(), tmod->nms().end(), nmsI.begin()) =
              Inf<MCFBlock::Index>();
          auto B = MCFB->get_B();
          for (MCFBlock::Index i = 0; i < NDfct.size(); i++)
            // NDfct[i] = B[nmsI[i]];
            bm->set(dgp->nodeFromId(i), B[nmsI[i]]);
          supply_changed = true;
          //MCFC::ChgDfcts(NDfct.data(), nmsI.data());
          return;
        }

        default:
          throw(std::invalid_argument("unknown MCFBlockSbstMod type"));
        }
      }

    }

    // any remaining Modification is plainly ignored, since it must be an
    // "abstract" Modification, which this Solver does not need to look at

  } // end( guts_of_poM )

/*--------------------------------------------------------------------------*/
/*------------------- End File MCFLemonSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/