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

//#include <math.h>

#include "MCFLemonSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
//using NS_GRVS = NetworkSimplex< SmartDigraph, double, double >;

/*template< typename Algo >
  struct Fields {};
 */  
  template<typename GR, typename V, typename C, typename TR>
  struct Fields< SMSppCapacityScaling< GR, V, C, TR > > {
    int f_factor;
  };

  template<typename GR, typename V, typename C, typename TR>  
  struct Fields< SMSppCostScaling< GR, V, C, TR > > {
    using CSMethod = typename SMSppCostScaling< GR, V, C, TR >::Method;
    int f_factor;
    CSMethod f_method;
  };
  template<typename GR, typename V, typename C>
  struct Fields< NetworkSimplex<GR, V, C> > {
    using NSPivotRule = typename NetworkSimplex<GR, V, C>::PivotRule;
    NSPivotRule f_pivot_rule;
  };
  template<typename GR, typename V, typename C>
  struct Fields< CycleCanceling<GR, V, C> > {
    using CCMethod = typename CycleCanceling<GR, V, C>::Method;
    CCMethod f_method;
  };

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register the various LEMONSolver< Alg , GR , V , C > to the Solver factory



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register MCFSolverState to the State factory

// SMSpp_insert_in_factory_cpp_0( MCFSolverState );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the various static maps

/*--------------------------------------------------------------------------*/
/*--------------------------------- SPECIALIZED CLASSES --------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename GR, typename V, typename C>
class MCFLemonSolver<NetworkSimplex, GR, V, C> : public CDASolver, Fields< NetworkSimplex<GR, V, C> >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

   

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   kUnEval = 0     compute() has not been called yet

   kUnbounded = kUnEval + 1     the model is provably unbounded
 *  @{ */
  const int kErrorStatus = -1;
/*
    
    kUnEval = 0     compute() has not been called yet

    kUnbounded = kUnEval + 1     the model is provably unbounded

    kInfeasible                  the model is provably infeasible

    kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

    kOK = 7         successful compute()
                    Any return value between kUnEval (excluded) and kOK
        (included) means that the object ran smoothly

    kStopTime = kOK + 1          stopped because of time limit

    kStopIter                    stopped because of iteration limit

    kError = 15     compute() stopped because of unrecoverable error
                    Any return value >= kError means that the object was
         forced to stop due to some error, e.g. of numerical nature

    kLowPrecision = kError + 1   a solution found but not provably optimal
    */



/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: does nothing special
 /** Void constructor: does nothing special, except verifying the template
  * arguments. */

 MCFLemonSolver( void ) : CDASolver(), Fields<NetworkSimplex<GR, V, C > >() {
  f_algo = NULL;


 /* static_assert( std::is_same< Algo< GR , V , C > ,
                               SMSppCapacityScaling< GR , V , C > >::value ||
		 std::is_same< Algo< GR , V , C > ,
                               SMSppCostScaling< GR , V , C >::value     ||
		 std::is_same< Algo< GR , V , C > ,
                               CycleCanceling< GR , V , C > >::value       ||
		 std::is_same< Algo< GR , V , C > ,
		               NetworkSimplex< GR , V , C > >::value ,
		 "Algo must be one of the LEMON algorithms");
  */
  }
 
  ~MCFLemonSolver( void ) {
    delete f_algo;
    delete dgp;
    
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/
 /* dblMaxTime = 0    maximum time for the next call to solve()

    dblRelAcc         relative accuracy for declaring a solution optimal

    dblAbsAcc          absolute accuracy for declaring a solution optimal

    dblUpCutOff        upper cutoff for stopping the algorithm

    dblLwCutOff        lower cutoff for stopping the algorithm

    dblRAccSol          maximum relative error in any reported solution

    dblAAccSol          maximum absolute error in any reported solution

    dblFAccSol          maximum constraint violation in any reported solution

    dblRAccDSol         maximum relative error in any dual solution

    dblAAccDSol         maximum absolute error in any dual solution

    dblFAccDSol         maximum absolute error in any dual solution

    dblLastParCDAS      first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/
 /// public enum "extending" int_par_type_CDAS to
 /// MCFLemonSolver< NetworkSimplex, C , V >

 enum int_par_type_LEMON_NS{
  kPivot = intLastParCDAS ///< the pivot rule
                /**< convenience value for easily allow derived classes
		       * to further extend the set of types of return codes */
  };  // end( int_par_type_LEMON_NS )
  
  ///public enum "extending" int_par_type_CDAS to
  /// MCFLemonSolver< CostScaling, C , V >

  enum int_par_type_LEMON_COS{
    kMethod = intLastParCDAS ///< the method for scaling the costs
  };

/* /// public enum "extending" int_par_type_CDAS to
  /// MCFLemonSolver< CapacityScaling, C , V >

  enum int_par_type_LEMON_CAS{
    kMethod = intLastParCDAS ///< the method for scaling the capacities
  };*/

  /// public enum "extending" int_par_type_CDAS to
  /// MCFLemonSolver< CycleCanceling, C , V >

 /* enum int_par_type_LEMON_CC{
    kMethod = intLastParCDAS ///< the method for run() function
  }; */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// public enum "extending" dbl_par_type_CDAS 
 /// to MCFLemonSolver< NetworkSimplex , C , V >

 enum dbl_par_type_LEMON_NS {
  dblBlockSizeFactor = dblLastParCDAS ///< the block size factor

  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };  // end( dbl_par_type_LEMON_NS )

  /// public enum "extending" dbl_par_type_CDAS 
  /// to MCFLemonSolver< CostScaling , C , V >

 /// public enum for the type of the solution
 /// to MCFLemonSolver< CostScaling , C , V >

 enum dbl_par_type_LEMON_COS{
  dblCostScaleFactor = dblLastParCDAS ///< the cost scale factor
 };

/// public enum for the type of the solution
/// to MCFLemonSolver< CapacityScaling , C , V >

enum dbl_par_type_LEMON_CAS{
  dblCapacityScaleFactor = dblLastParCDAS ///< the capacity scale factor
};

/// public enum for the type of the solution
/// to MCFLemonSolver< CycleCanceling , C , V >

enum dbl_par_type_LEMON_CC{
  dblCycleCancelingFactor = dblLastParCDAS ///< the cycle canceling factor
};


 enum sol_type
 {
  UNSOLVED, //= NULL, ///< the problem has not been solved yet
  OPTIMAL,           ///< the problem has been solved
  KSTOPTIME, //= NULL,     ///< the problem has been stopped because of time limit
  INFEASIBLE,   ///< the problem is provably infeasible
  UNBOUNDED,    ///< the problem is provably unbounded
  KERROR //= NULL         ///< the problem has been stopped because of unrecoverable error
  };               // end( sol_type )

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other parameter of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 void set_Block(Block *block) override
    {
     
      
      if (block == f_Block) // actually doing nothing
        return;             // cowardly and silently return
      
      delete f_algo;
      f_algo=NULL;

      Solver::set_Block(block); // attach to the new Block

      if (block)
      { // this is not just resetting everything
        auto MCFB = dynamic_cast<MCFBlock *>(block);
        if (!MCFB)
          throw(std::invalid_argument(
              "MCFSolver:set_Block: block must be a MCFBlock"));

        bool owned = MCFB->is_owned_by(f_id);
        if ((!owned) && (!MCFB->read_lock()))
          throw(std::logic_error("cannot acquire read_lock on MCFBlock"));

        // load the new MCFBlock into the :MCFClass object
        // TODO: change MCFC function to Algo function.
        // TODO: convert array from MCFB functions to Map for Algo functions.
        dgp = new GR;
        dgp->clear();

        dgp->reserveNode( MCFB->get_MaxNNodes() );
        MCFBlock::Index n = MCFB->get_NNodes();

        for( MCFBlock::Index i = 0; i < n; ++i)
          dgp->addNode();

        dgp->reserveArc( MCFB->get_MaxNArcs() );
        MCFBlock::Index m = MCFB->get_NArcs();

        MCFBlock::c_Subset & sn = MCFB->get_SN();
        MCFBlock::c_Subset & en = MCFB->get_EN();

        for( MCFBlock::Index i = 0; i < m; ++i){
          dgp->addArc( dgp->nodeFromId( sn[i] - 1) , dgp->nodeFromId( en[i] - 1 ) );
        }


        
        f_algo = new NetworkSimplex< GR , V, C >(*dgp);

        using MCFArcMapV = typename GR::template ArcMap< V >;
        using MCFNodeMapV = typename GR::template NodeMap< V >;

        if(!MCFB->get_U().empty())
        {
          MCFArcMapV um(*dgp);
          MCFBlock::c_Vec_FNumber & u = MCFB->get_U();
          for( MCFBlock::Index i = 0; i < m; ++i){
          um.set( dgp->arcFromId(i), u[i]);
          }
          f_algo->upperMap(um);
        }

        if(!MCFB->get_C().empty())
        {
          MCFArcMapV cm(*dgp);
          MCFBlock::c_Vec_FNumber & c = MCFB->get_C();
          for( MCFBlock::Index i = 0; i < m; ++i){
            cm.set( dgp->arcFromId(i), c[i]);
          }
          f_algo->costMap(cm);

        }

        if(!MCFB->get_B().empty())
        {

          MCFNodeMapV bm(*dgp);
          MCFBlock::c_Vec_FNumber & b = MCFB->get_B();
          for( MCFBlock::Index i = 0; i < n; i++){
            bm.set( dgp->nodeFromId(i), -b[i]);
          }
          f_algo->supplyMap(bm);

        }

        
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the data
        //       changes, some of the deductions (say, reducing the capacity but
        //       of some arcs) may no longer be correct and they should be undone,
        //       there isn't any proper way to handle this. Thus, PreProcess() has
        //       to be disabled for now; maybe later on someone will take care to
        //       make this work (or maybe not).
        // MCFC::PreProcess();

        // once done, read_unlock the MCFBlock (if it was read-lock()-ed)
        if (!owned)
          MCFB->read_unlock();

        // TODO: maybe log it
        //delete dgp;
        
      }
      
    } // end( set_Block )

    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
   /* template < NetworkSimplex ,  GR,  V,  C >
    void set_par(idx_type par, int value) 
    {
      if(par == kPivot){
        switch(par){
          case PivotRule::FIRST_ELIGIBLE:
              FirstEligiblePivotRule(f_algo);
              break;
          case PivotRule::BEST_ELIGIBLE:
              BestEligiblePivotRule(f_algo);
              break;
          case PivotRule::BLOCK_SEARCH:
              BlockSearchPivotRule(f_algo);
              break;
          case PivotRule::ALTERING_LIST:
              AlteringListPivotRule(f_algo);
              break;
          default:
            break;          
        }

        return;
      }

      CDASolver::set_par(par, value);
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    /*void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
      ;
      //  MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */
    /// (try to) solve the MCF encoded in the MCFBlock 
        int compute( bool changedvars = true ) override {
        const static std::array<int, 6> LemonStatus_2_MCFstatus = {
                kErrorStatus, OPTIMAL, kErrorStatus , INFEASIBLE,
                UNBOUNDED, kErrorStatus};
        
        const static std::array<int, 6> MCFstatus_2_sol_type = {
                kUnEval, Solver::kOK, kStopTime, kInfeasible, Solver::kUnbounded,
                Solver::kError};

        lock(); // first of all, acquire self-lock

        if (!f_Block)            // there is no [MCFBlock] to solve
                return (kBlockLocked); // return error

        bool owned = f_Block->is_owned_by(f_id); // check if already locked
        if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
                return (kBlockLocked);                 // return error on failure

        // while [read_]locked, process any outstanding Modification
        //TODO: ensure that modification are actually processed for MCFLemonSolver.
        //process_outstanding_Modification();

        if (!f_dmx_file.empty())
        { // if so required
                // output the current instance (after the changes) to a DMX file
                std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
                if (!ProbFile.is_open())
                throw(std::logic_error("cannot open DMX file " + f_dmx_file));

                //WriteMCF(ProbFile);
                writeDimacsMat(ProbFile, *dgp);
                ProbFile.close();
        }

        if (!owned)               // if the [MCF]Block was actually read_locked
                f_Block->read_unlock(); // read_unlock it

        auto start = chrono::system_clock::now();
        
        //if(Fields<NetworkSimplex< GR, V, C > >::f_pivot_rule != NULL){
        this->status = f_algo->run(Fields<NetworkSimplex< GR, V, C > >::f_pivot_rule);
       // }else{
                //Build f_method and execute run() method
      //  }
        
        auto end = chrono::system_clock::now();

        chrono::duration< double > elapsed = end - start;
        ticks = elapsed.count();
        if(LemonStatus_2_MCFstatus[this->get_status()] == kErrorStatus){
                return Solver::kError;
        }

        
        
        unlock(); // release self-lock

        // now give out the result: note that the vector MCFstatus_2_sol_type[]
        // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
        // hence the returned status has to be shifted by + 1
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    


    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    double get_elapsed_time(void) const override
    {
      return (this->ticks);
    }



    /*--------------------------------------------------------------------------*/
    //Return the lower bound solution(optimal) for the problem
    OFValue get_lb(void) override { return OFValue(f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //Return the upper bound solution(optimal) for the problem
    OFValue get_ub(void) override { return OFValue(f_algo->totalCost()); }

    /*--------------------------------------------------------------------------*/
    //TODO: change MCFC function to Algo function. DONE
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( NetworkSimplex< GR , V , C >::ProblemType::OPTIMAL ):
      case( NetworkSimplex< GR , V , C >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( NetworkSimplex< GR , V , C >::ProblemType::OPTIMAL ):
      case( NetworkSimplex< GR , V , C >::ProblemType::INFEASIBLE ):
       return( true );
      default:
       return( false );
      }
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool is_var_feasible( void ) override { return( true ); }

     virtual bool is_dual_feasible( void ) override { return( true ); }
    */
    /*--------------------------------------------------------------------------*/
    /// write the "current" flow in the x ColVariable of the MCFBlock
    /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
     * the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the dual
     * solution". In all other cases, the flow solution is saved. */

    void get_var_solution(Configuration *solc = nullptr) override
    {/*
      if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_FNumber X(MCFB->get_NArcs());
      this->MCFGetX(X.data());
      MCFB->set_x(X.begin());
   */}
  
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the "current" dual solution in the Constraint of the MCFBlock
    /** Write the "current" dual solution, i.e., node potentials and flow
     * reduced costs, in the dual variables of the Constraint (respectively,
     * the flow conservation constraints and bound ones) of the MCFBlock. To
     * keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 1, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the primal
     * solution". In all other cases, the flow solution is saved. */
    
    void get_dual_solution(Configuration *solc = nullptr) override
    {
    /*  if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_CNumber Pi(MCFB->get_NNodes());
      this->MCFGetPi(Pi.data());
      MCFB->set_pi(Pi.begin());

      MCFBlock::Vec_FNumber RC(MCFB->get_NArcs());
      this->MCFGetRC(RC.data());
      MCFB->set_rc(RC.begin());
    */}
    
    /*--------------------------------------------------------------------------*/

   // bool new_var_solution(void) override { return (this->HaveNewX()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   // bool new_dual_solution(void) override { return (this->HaveNewPi()); }

    /*--------------------------------------------------------------------------*/
    /*
     virtual void set_unbounded_threshold( const OFValue thr ) override { }
    */

    /*--------------------------------------------------------------------------*/

    bool has_var_direction(void) override { return (true); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_direction(void) override { return (true); }

    /*--------------------------------------------------------------------------*/
    /// write the current direction in the x ColVariable of the MCFBlock
    /** Write the unbounded primal direction, i.e., augmenting cycle with
     * negative cost and unbounded capacity, in the x ColVariable of the
     * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is done*,
     * since the Configuration is meant to say "only save/map the dual
     * information". In all other cases, the direction (cycle) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_var_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      throw(std::logic_error(
          "MCFSolver::get_var_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnbCycl()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular, if
     * solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value == 1,
     * then *nothing is done*, since the Configuration is meant to say "only
     * save/map the primal information". In all other cases, the direction (cut)
     * is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_dual_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      throw(std::logic_error(
          "MCFSolver::get_dual_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnfCut()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool new_var_direction( void ) override { return( false ); }

     virtual bool new_dual_direction( void ) override{ return( false ); }
    */
    /** @} ---------------------------------------------------------------------*/
    /*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
    /*--------------------------------------------------------------------------*/

    /*
     virtual bool is_dual_exact( void ) const override { return( true ); }
    */

    /*--------------------------------------------------------------------------*/
    /// "publicize" MCFClass::WriteMCF
    /** Make the method
     *
     *      void WriteMCF( ostream &oStrm , int frmt = 0 )
     *
     * of the base (private) MCFClass public, so that it can be freely used. */
    //TODO: change MCFC function to Algo function.
    //using MCFC::WriteMCF;

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double parameters. If
     * this is the case, it will have to specialize the following methods to
     * handle them. The general definition just handles the case of the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */


    [[nodiscard]] idx_type get_num_int_par(void) const override
    {/*
      return (CDASolver::get_num_int_par() + 1);
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {/*
      return (CDASolver::get_num_dbl_par());
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {/*
      return (CDASolver::get_num_str_par() + 1);
    */}
    
    /*--------------------------------------------------------------------------*/
    
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {/*
      if (par == intLastParCDAS)
        return (MCFClass::kYes);

      return (CDASolver::get_dflt_int_par(par));
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {/*
      return (CDASolver::get_dflt_dbl_par(par));
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {
      static const std::string _empty;
      if (par == strLastParCDAS)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    [[nodiscard]] int get_int_par(idx_type par) const override
    {/*
      if (Solver_2_MCFClass_int[par] >= 0)
      {
        int val;
        this->GetPar(Solver_2_MCFClass_int[par], val);
        return (val);
      }

      return (get_dflt_int_par(par));
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {/*
      if (Solver_2_MCFClass_dbl[par] >= 0)
      {
        double val;
        this->GetPar(Solver_2_MCFClass_dbl[par], val);
        return (val);
      }

      return (get_dflt_dbl_par(par));
    */}
        
    /*--------------------------------------------------------------------------*/
    
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {/*
      if (name == "kReopt")
        return (kReopt);

      return (CDASolver::int_par_str2idx(name));
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {/*
      return (CDASolver::dbl_par_str2idx(name));
    */}
      
    /*--------------------------------------------------------------------------*/
    
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {/*
      static const std::string my_name = "kReopt";

      if (idx == intLastParCDAS)
        return (my_name);

      return (CDASolver::int_par_idx2str(idx));
    */}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {/*
      return (CDASolver::dbl_par_idx2str(idx));
    */}
    
    /** @} ---------------------------------------------------------------------*/
    /*------------ METHODS FOR HANDLING THE State OF THE MCFLemonSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the State of the MCFSolver
     *  @{ */

    //[[nodiscard]] State *get_State(void) const override;

    /*--------------------------------------------------------------------------*/

    //void put_State(const State &state) override;

    /*--------------------------------------------------------------------------*/

    //void put_State(State &&state) override;

    /** @} ---------------------------------------------------------------------*/
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

    /*void add_Modification(sp_Mod &mod) override
    {
      if (std::dynamic_pointer_cast<const NBModification>(mod))
      {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        //TODO: change MCFC function to Algo function.
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        MCFC::LoadNet(MCFB->get_MaxNNodes(), MCFB->get_MaxNArcs(),
                      MCFB->get_NNodes(), MCFB->get_NArcs(),
                      MCFB->get_U().empty() ? nullptr : MCFB->get_U().data(),
                      MCFB->get_C().empty() ? nullptr : MCFB->get_C().data(),
                      MCFB->get_B().empty() ? nullptr : MCFB->get_B().data(),
                      MCFB->get_SN().data(), MCFB->get_EN().data());
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the
        //       data changes, some of the deductions (say, reducing the capacity
        //       of some arcs) may no longer be correct and they should be undone,
        //       but there isn't any proper way to handle this. Thus, PreProcess()
        //       has to be disabled for now; maybe later on someone will take
        //       care to make this work (or maybe not).
        // MCFC::PreProcess();
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
      }
      else
        push_back(mod);
    }*/

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void process_outstanding_Modification(void);

    void guts_of_poM(c_p_Mod mod);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS  ---------------------------*/
    /*--------------------------------------------------------------------------*/

    const static std::vector<int> Solver_2_MCFClass_int;
    // the (static const) map between Solver int parameters and MCFClass ones

    const static std::vector<int> Solver_2_MCFClass_dbl;
    // the (static const) map between Solver int parameters and MCFClass ones

    std::string f_dmx_file; 
    // string for DMX file output


    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

     /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    NetworkSimplex< GR , V , C > * f_algo;  //Algorithm used by lemon

    GR* dgp;  //Rapresentation of directed graph
    V* value;   //Type of value of node
    C* costs;  //Type of costs of arcs
    int counter=0;
    int status = UNSOLVED;  //Variable used in compute function for getting status
    typename NetworkSimplex< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    
    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver )

  SMSpp_insert_in_factory_cpp_0_t(
 MCFLemonSolver< NetworkSimplex , SmartDigraph , double , double >);




/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/* Managing parameters for MCFSimplex---------------------------------------*/
/*
 * MCFSimplex has the following extra parameters:
 *
 * - kAlgPrimal     parameter to set algorithm (Primal/Dual):
 * - kAlgPricing    parameter to set algorithm of pricing
 * - kNumCandList   parameter to set the number of candidate list for
 *                  Candidate List Pivot method
 * - kHotListSize   parameter to set the size of Hot List
 *
 * (plus, actually, a few about Quadratic MCF that are not relevent here).
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */

#ifdef HAVE_MFSMX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< MCFSimplex >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 MCFSimplex::kAlgPrimal ,
 MCFSimplex::kAlgPricing ,
 MCFSimplex::kNumCandList ,
 MCFSimplex::kHotListSize
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< MCFSimplex >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFSimplex >::get_num_int_par( void ) const
{
 //return( CDASolver::get_num_int_par() + 5 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFSimplex >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< MCFSimplex >::get_dflt_int_par( idx_type par ) const
{/*
 static const std::array< int , 5 > my_dflt_int_par = { MCFClass::kYes ,
		  MCFClass::kYes , MCFSimplex::kCandidateListPivot , 0 , 0 };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< MCFSimplex >::get_dflt_dbl_par( const idx_type par )
 const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFSimplex >::int_par_str2idx(
					     const std::string & name ) const
{/*
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAlgPrimal" )
  return( intLastParCDAS + 1 );
 if( name == "kAlgPricing" )
  return( intLastParCDAS + 2 );
 if( name == "kNumCandList" )
  return( intLastParCDAS + 3 );
 if( name == "kHotListSize" )
  return( intLastParCDAS + 4 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFSimplex >::dbl_par_str2idx(
					const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< MCFSimplex >::int_par_idx2str( idx_type idx )
 const {/*
 static const std::array< std::string , 5 > my_int_pars_str = {
  "kReopt" , "kAlgPrimal" , "kAlgPricing" , "kNumCandList" , "kHotListSize"
  };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< MCFSimplex >::dbl_par_idx2str( idx_type idx )
 const { }
 */

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for RelaxIV -----------------------------------------*/
/*
 * RelaxIV has the following extra parameters:
 *
 * - kAuction     the auction/shortest paths initialization is used
 *
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */

#ifdef HAVE_RELAX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 RelaxIV::kAuction
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_int_par( void ) const
{
 //return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< RelaxIV >::get_dflt_int_par( const idx_type par ) const
{/*
 static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
						        MCFClass::kYes };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< RelaxIV >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::int_par_str2idx(
					     const std::string & name ) const
{/*
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAuction" )
  return( intLastParCDAS + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::dbl_par_str2idx(
					 const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< RelaxIV >::int_par_idx2str( idx_type idx )
 const {
 /*static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
								"kAuction" };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< RelaxIV >::dbl_par_idx2str( idx_type idx )
 const { }
 */

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for CPLEX -------------------------------------------*/

#ifdef HAVE_CPLEX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 MCFCplex::kQPMethod
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_dbl =
{
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_int_par( void ) const {
 //return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< MCFCplex >::get_dflt_int_par( idx_type par ) const {
 /*static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
                                                        MCFClass::kYes };
 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 double MCFSolver< MCFCplex >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::int_par_str2idx(
					   const std::string & name ) const {
 /*if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kQPMethod" )
  return( kQPMethod + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< MCFCplex >::int_par_idx2str( idx_type idx )
 const {
 /*static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
                                                                "kQPMethod" };
 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< MCFCplex >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for SPTree ------------------------------------------*/

#ifdef HAVE_SPTRE

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter,        // intMaxIter
 -1,                        // intMaxSol
 -1,                        // intLogVerb
 -1,                        // intMaxDSol
 MCFClass::kReopt,          // intLastParCDAS
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime,        // dblMaxTime
 -1,                        // dblRelAcc
 MCFClass::kEpsFlw,         // dblAbsAcc
 -1,                        // dblUpCutOff
 -1,                        // dblLwCutOff
 -1,                        // dblRAccSol
 -1,                        // dblAAccSol
 -1,                        // dblFAccSol
 -1,                        // dblRAccDSol
 MCFClass::kEpsCst,         // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_int_par( void ) const {
 /*return( CDASolver::get_num_int_par() + 1 );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< SPTree >::get_dflt_int_par( idx_type par ) const {
 /*static const std::array< int , 1 > my_dflt_int_par = { MCFClass::kYes };

 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< SPTree >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::int_par_str2idx(
					  const std::string & name ) const {
 /*if( name == "kReopt" )
  return( intLastParCDAS );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< SPTree >::int_par_idx2str( idx_type idx )
 const {/*
 static const std::array< std::string , 1 > my_int_pars_str = { "kReopt" };

 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 const std::string & MCFSolver< SPTree >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
 