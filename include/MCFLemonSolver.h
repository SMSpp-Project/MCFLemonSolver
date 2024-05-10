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

//!!
#include <concepts>

#include <lemon/list_graph.h>

#include <lemon/smart_graph.h>

#include <lemon/dimacs.h>

#include <lemon/concepts/digraph.h>


/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/
/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace lemon;
 using namespace lemon::concepts;
 using namespace std;

/*--------------------------------------------------------------------------*/
/*-------------------------- TEMPLATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFLemonSolver_TYPES Template Types in MCFLemonSolver.h
 *
 * Algorithms in the LEMON projects are template over at least three types:
 *
 * - GR, which is the type of graph (DISCUSS THE POSSIBILITIES);
 *
 * - V, which is the type of flows / deficits; typically, double can be used
 *   for maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *   maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
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
 * types of LEMON graphs, i.e., (DISCUSS THE POSSIBILITIES);
 *
 * Similarly, the concept LEMONAlgorithm is defined that only allows all
 * possible types of LEMON algorithms, i.e., [SMSpp]CapacityScaling,
 * [SMSpp]CostScaling, CycleCanceling, and NetworkSimplex.
 *
 *  @{ */

 /// concept for "one of the LEMON graphs"
 template< typename Type >
  concept LEMONGraph =
   std::is_same< Type , SmartDigraph >::value   ||
   std::is_same< Type , StaticDigraph >::value;

 //!! std::is_base_of< Graph , Type >::value;
 /*!!
   std::is_same< Type , CompactDigraph >::value ||
   std::is_same< Type , FullDigraph >::value;
   std::is_same< Type , GridGraph >::value      ||
   std::is_same< Type , HypercubeGraph >::value ||
   !!*/

 /// CapacityScaling algorithm using the default trait
 template< LEMONGraph GR , typename V , typename C, typename TR >
 class SMSppCapacityScaling : public
  CapacityScaling< GR , V , C , CapacityScalingDefaultTraits< GR , V , C > >
  {};

 /// CostScaling algorithm using the default trait
 template < LEMONGraph GR , typename V , typename C, typename TR >
 class SMSppCostScaling : public
  CostScaling< GR , V , C , CostScalingDefaultTraits< GR , V , C > >
  {};

/** @} ---------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFLemonSolver_CLASSES Classes in MCFLemonSolver.h
 *  @{ */

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
 * The MCFLemonSolver is template over four different types:
 *
 * - GR, which is the type of graph (DISCUSS THE POSSIBILITIES);
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
 *   C) implemented in the LEMON class (DISCUSS THE POSSIBILITIES);
 *   Note that scaling-type algorithms may behave in different ways according
 *   to which combination of V and C is used, and there are different
 *   "traits" for this which are another scaling parameter. However, in order
 *   to make MCFLemonSolver class template over always the same number of
 *   template parameters we fix the use of the default trait, which is why
 *   SMSppCapacityScaling and SMSppCostScaling are defined (as template over
 *   < GR , V , C >) that are meant to be used instead of the original
 *   CapacityScaling and CostScaling. */

template< template< typename , typename , typename > typename Algo ,
          typename GR , typename V , typename C >
  requires LEMONGraph< GR >
class MCFLemonSolver : public CDASolver
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

 MCFLemonSolver( void ) : CDASolver() {

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

 
 
 Algo< GR , V , C > * f_algo;  //Algorithm used by lemon

 GR digraph;  //Rapresentation of directed graph
 V* value;   //Type of value of node
 C* costs;  //Type of costs of arcs

 int status = UNSOLVED;  //Variable used in compute function for getting status
 typename Algo< GR , V , C >::ProblemType status_2_pType;
 //Status of compute() method
 
 double ticks;  //Elaped time in ticks for compute() method

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
 /// public enum "extending" int_par_type_CDAS to MCFSolver

 enum int_par_type_MCFS {
  kReopt = intLastParCDAS, ///< whether or not to reoptimize
  intLastParMCF       ///< first allowed parameter value for derived classes
                      /**< convenience value for easily allow derived classes
		       * to further extend the set of types of return codes */
  };  // end( int_par_type_MCFS )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// public enum "extending" dbl_par_type_CDAS to MCFSolver

 enum dbl_par_type_MCFS {
  dblLastParMCF = dblLastParCDAS
  ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };  // end( dbl_par_type_MCFS )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// public enum "extending" str_par_type_CDAS to MCFSolver

 enum str_par_type_MCFS {
  strDMXFile = strLastParCDAS, ///< DMX filename to output the instance
  strLastParMCF        ///< first allowed parameter value for derived classes
                       /**< convenience value for easily allow derived classes
			* to further extend the set of types of return codes */
  };  // end( dbl_par_type_MCFS )

 /// public enum for the type of the solution

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


        digraph.reserveNode(MCFB->get_MaxNNodes());
        digraph.reserveArc(MCFB->get_MaxNArcs());

        for(int i = 0; i < MCFB->get_NNodes();i++){
          //Digraph::Node n;
          digraph.addNode();    
        }

        for(int i = 0; i < MCFB->get_NArcs(); i++){
          digraph.addArc(digraph.nodeFromId(MCFB->get_SN(i)), digraph.nodeFromId(MCFB->get_EN(i)));
        }

        auto dgp = new GR;

        dgp->reserveNode( MCFB->get_MaxNNodes() );
        auto n = MCFB->get_NNodes();

        for( MCFBlock::Index i = 0; i < n; ++i)
          dgp->addNode();

        dgp->reserveArc( MCFB->get_MaxNArcs() );
        auto m = MCFB->get_NArcs();

        MCFBlock::c_Subset & sn = MCFB->get_SN();
        MCFBlock::c_Subset & en = MCFB->get_EN();

        for( MCFBlock::Index i = 0; i < m; ++i){
          dgp->addArc( digraph.nodeFromId( sn[i] ), digraph.nodeFromId( en[i] ));
        }
        
        
        f_algo = new Algo< GR , V, C >(*dgp);

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
          for( MCFBlock::Index i = 0; i < m; ++i){
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
      }
    } // end( set_Block )

    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
    /* //TODO: change MCFC function to Algo function.
    void set_par(idx_type par, int value) override
    {
      if (Solver_2_MCFClass_int[par] >= 0)
        MCFC::SetPar(Solver_2_MCFClass_int[par], int(value));
    }
    */
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /*//TODO: change MCFC function to Algo function.
    void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
        MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }
    */
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    void set_par(idx_type par, const std::string &value) override
    {
      if (par == strDMXFile)
        f_dmx_file = value;
    }
    

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */

    /// (try to) solve the MCF encoded in the MCFBlock
    int compute( bool changedvars = true ) override
    {
      //TODO: map errors from Algo to MCFSolver. DONE
      //TODO: resolve for NULL values in LemonStatus_2_MCFstatus.
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
        writeDimacsMat(ProbFile, digraph);
        ProbFile.close();
      }

      if (!owned)               // if the [MCF]Block was actually read_locked
        f_Block->read_unlock(); // read_unlock it

      // ensure the timer exists (or reset it)
      //TODO: implement timer with ctime. feature not present in Algo.
      
      //this->MCFC::SetMCFTime();

      // then (try to) solve the MCF
      //TODO: change MCFC function to Algo function. DONE
      //Probably running the method run of the Algo class. DONE

      auto start = chrono::system_clock::now();

      //this->MCFC::SolveMCF();
      this->status = f_algo->run();

      auto end = chrono::system_clock::now();

      chrono::duration< double > elapsed = end - start;
      ticks = elapsed.count();
      int status = this->get_status();
      if(LemonStatus_2_MCFstatus[status] == kErrorStatus){
        return Solver::kError;
      }

      

      unlock(); // release self-lock

      // now give out the result: note that the vector MCFstatus_2_sol_type[]
      // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
      // hence the returned status has to be shifted by + 1
      //TODO: change MCFC function to Algo function.
      //Da rivedere, this->get_status() potrebbe non essere corretta.
      return (MCFstatus_2_sol_type[status]);
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
      case( Algo< GR , V , C >::ProblemType::OPTIMAL ):
      case( Algo< GR , V , C >::ProblemType::UNBOUNDED ):  
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
      case( Algo< GR , V , C >::ProblemType::OPTIMAL ):
      case( Algo< GR , V , C >::ProblemType::INFEASIBLE ):
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
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] const std::string &get_str_par(idx_type par)
        const override
    {
      if (par == strDMXFile)
        return (f_dmx_file);

      return (get_dflt_str_par(par));
    }
    
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
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] idx_type str_par_str2idx(const std::string &name)
        const override
    {/*
      if (name == "strDMXFile")
        return (strDMXFile);

      return (CDASolver::str_par_str2idx(name));
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
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    [[nodiscard]] const std::string &str_par_idx2str(idx_type idx)
        const override
    {/*
      static const std::string my_name = "strDMXFile";

      if (idx == strDMXFile)
        return (my_name);

      return (CDASolver::str_par_idx2str(idx));
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

    std::string f_dmx_file; // string for DMX file output

    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

    /*--------------------------------------------------------------------------*/

  }; // end( class MCFSolver )

  /*--------------------------------------------------------------------------*/
  /*------------------------- CLASS MCFSolverState ---------------------------*/
  /*--------------------------------------------------------------------------*/
  /// class to describe the "internal state" of a MCFSolver
  /** Derived class from State to describe the "internal state" of a MCFSolver,
   *  i.e., a MCFClass::MCFState (*). Since MCFClass::MCFState does not allow
   *  serialization, all that part does not work.  */

  class MCFLemonState : public State
  {
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

  public:
    /*------------- CONSTRUCTING AND DESTRUCTING MCFSolverState ----------------*/

    /// constructor, doing everything or nothing.
    /** Constructor of MCFSolverState. If provided with a pointer to a MCFSolver
     * it immediately copies its "internal state", which is the only way in which
     * the MCFSolverState can be initialised out of an existing MCFSolver. If
     * nullptr is passed (as by default), then an "empty" MCFSolverState is
     * constructed that can only be filled by calling deserialize().
     *
     * Note: to avoid having to duplicate the SMSpp_insert_in_factory_cpp call
     *       for every MCFClass, the pointer is directly that of a MCFClass,
     *       since every MCFSolver derives from a :MCFClass and we only need
     *       access to MCFGetState(). */

    /*MCFSolverState(MCFClass *mcfc = nullptr) : State()
    {
      f_state = mcfc ? mcfc->MCFGetState() : nullptr;
    }*/

    /*--------------------------------------------------------------------------*/
    /// de-serialize a MCFSolverState out of netCDF::NcGroup
    /** Should de-serialize a MCFSolverState out of netCDF::NcGroup, but in
     * fact it does not work. */

    /*void deserialize(const netCDF::NcGroup &group) override
    {
      f_state = nullptr;
    }*/

    /*--------------------------------------------------------------------------*/
    /// destructor

    //virtual ~MCFSolverState() { delete f_state; }

    /*---------- METHODS DESCRIBING THE BEHAVIOR OF A MCFSolverState -----------*/

    /// serialize a MCFSolverState into a netCDF::NcGroup
    /** The method should serialize the MCFSolverState into the provided
     * netCDF::NcGroup, so that it can later be read back by deserialize(), but
     * in fact it does not work.*/

    void serialize(netCDF::NcGroup &group) const override {}

    /*-------------------------------- FRIENDS ---------------------------------*/

    template <class MCFC>
    friend class MCFSolver; // make MCFSolver friend

    /*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

  protected:
    /*-------------------------- PROTECTED METHODS -----------------------------*/

    void print(std::ostream &output) const override
    {
      output << "MCFSolverState [" << this << "]";
    }

    /*--------------------------- PROTECTED FIELDS -----------------------------*/

   // MCFClass::MCFStatePtr f_state; ///< the (pointer to) MCFState

    /*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

  private:
    /*---------------------------- PRIVATE FIELDS ------------------------------*/

    SMSpp_insert_in_factory_h;

    /*--------------------------------------------------------------------------*/

  }; // end( class( MCFSolverState ) )

  /** @} end( group( MCFSolver_CLASSES ) ) */
  /*--------------------------------------------------------------------------*/
  /*------------------- inline methods implementation ------------------------*/
  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  State *MCFSolver<Algo>::get_State(void) const
  {
    return (new MCFSolverState(const_cast<MCFSolver<MCFC> *>(this)));
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::put_State(const State &state)
  {
    // if state is not a const MCFSolverState &, exception will be thrown
    auto s = dynamic_cast<const MCFSolverState &>(state);

    this->MCFPutState(s.f_state);
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::put_State(State &&state)
  {
    // if state is not a MCFSolverState &&, exception will be thrown
    auto s = dynamic_cast<MCFSolverState &&>(state);

    this->MCFPutState(s.f_state);
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::process_outstanding_Modification(void)
  {
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
  }*/ // end( MCFSolver::process_outstanding_Modification )

  /*--------------------------------------------------------------------------*/

  /*template <class Algo>
  void MCFSolver<Algo>::guts_of_poM(c_p_Mod mod)
  {
    auto MCFB = static_cast<MCFBlock *>(f_Block);
    */
    // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /* This requires to patiently sift through the possible Modification types
     * to find what this Modification exactly is, and call the appropriate
     * method of MCFClass. */

    // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*if (auto tmod = dynamic_cast<const GroupModification *>(mod))
    {
      for (const auto &submod : tmod->sub_Modifications())
        guts_of_poM(submod.get());

      return;
    }*/

    // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /* Note: in the following we can assume that C, B and U are nonempty. This
     * is because they can be empty only if they are so when the object is
     * loaded. But if a Modification has been issued they are no longer empty (a
     * Modification changin nothing from the "empty" state is not issued). */

    /*if (auto tmod = dynamic_cast<const MCFBlockRngdMod *>(mod))
    {
      auto rng = tmod->rng();

      switch (tmod->type())
      {
      case (MCFBlockMod::eChgCost):
        if (rng.second == rng.first + 1)
        {
          if (!MCFB->is_deleted(rng.first))
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCost(rng.first, MCFB->get_C(rng.first));
        }
        else
        {
          if (std::any_of(MCFB->get_C().data() + rng.first,
                          MCFB->get_C().data() + rng.second,
                          [](auto ci)
                          { return (std::isnan(ci)); }))
          {
            MCFBlock::Vec_CNumber NCost(MCFB->get_C().data() + rng.first,
                                        MCFB->get_C().data() + rng.second);
            for (auto &ci : NCost)
              if (std::isnan(ci))
                ci = 0;
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCosts(NCost.data(), nullptr, rng.first, rng.second);
          }
          else
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCosts(MCFB->get_C().data() + rng.first, nullptr,
                           rng.first, rng.second);
        }
        return;

      case (MCFBlockMod::eChgCaps):
        if (rng.second == rng.first + 1)
        {
          if (!MCFB->is_deleted(rng.first))
            //TODO: change MCFC function to Algo function.
            MCFC::ChgUCap(rng.first, MCFB->get_U(rng.first));
        }
        else
          //TODO: change MCFC function to Algo function.
          MCFC::ChgUCaps(MCFB->get_U().data() + rng.first, nullptr,
                         rng.first, rng.second);
        return;

      case (MCFBlockMod::eChgDfct):
        if (rng.second == rng.first + 1)
          //TODO: change MCFC function to Algo function.
          MCFC::ChgDfct(rng.first, MCFB->get_B(rng.first));
        else
          //TODO: change MCFC function to Algo function.
          MCFC::ChgDfcts(MCFB->get_B().data() + rng.first, nullptr,
                         rng.first, rng.second);
        return;
      //TODO: change MCFC function to Algo function.
      case (MCFBlockMod::eOpenArc):
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

      case (MCFBlockMod::eAddArc):
      {
        auto ca = MCFB->get_C(rng.first);
        //TODO: change MCFC function to Algo function.
        auto arc = MCFC::AddArc(MCFB->get_SN(rng.first),
                                MCFB->get_EN(rng.first),
                                MCFB->get_U(rng.first),
                                std::isnan(ca) ? 0 : ca);
        if (arc != rng.first)
          throw(std::logic_error("name mismatch in AddArc()"));
        return;
      }

      case (MCFBlockMod::eRmvArc):
        //TODO: change MCFC function to Algo function.
        MCFC::DelArc(rng.second - 1);
        return;

      default:
        throw(std::invalid_argument("unknown MCFBlockRngdMod type"));
      }
    }

    // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (auto tmod = dynamic_cast<const MCFBlockSbstMod *>(mod))
    {
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
      }

      // have to InINF-terminate the vector of indices (damn!)
      // meanwhile, when appropriate remove the indices of deleted
      // arcs for which the operations make no sense;
      ;
      switch (tmod->type())
      {
      case (MCFBlockMod::eChgCost):
      {
        MCFBlock::Subset nmsI;
        nmsI.reserve(tmod->nms().size() + 1);
        MCFBlock::Vec_CNumber NCost;
        NCost.reserve(tmod->nms().size());
        auto &C = MCFB->get_C();
        for (auto i : tmod->nms())
          if (auto ci = C[i]; !std::isnan(ci))
          {
            NCost.push_back(ci);
            nmsI.push_back(i);
          }
        nmsI.push_back(Inf<MCFBlock::Index>());
        //TODO: change MCFC function to Algo function.
        MCFC::ChgCosts(NCost.data(), nmsI.data());
        return;
      }

      case (MCFBlockMod::eChgCaps):
      {
        MCFBlock::Subset nmsI;
        nmsI.reserve(tmod->nms().size() + 1);
        MCFBlock::Vec_FNumber NCap;
        NCap.reserve(tmod->nms().size());
        auto &C = MCFB->get_C();
        auto &U = MCFB->get_U();
        for (auto i : tmod->nms())
          if (!std::isnan(C[i]))
          {
            NCap.push_back(U[i]);
            nmsI.push_back(i);
          }
        nmsI.push_back(Inf<MCFBlock::Index>());
        //TODO: change MCFC function to Algo function.
        MCFC::ChgUCaps(NCap.data(), nmsI.data());
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
          NDfct[i] = B[nmsI[i]];
        //TODO: change MCFC function to Algo function.
        MCFC::ChgDfcts(NDfct.data(), nmsI.data());
        return;
      }

      default:
        throw(std::invalid_argument("unknown MCFBlockSbstMod type"));
      }
    }

    // any remaining Modification is plainly ignored, since it must be an
    // "abstract" Modification, which this Solver does not need to look at

  }*/ // end( guts_of_poM )

  /*--------------------------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* MCFLemonSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFLemonSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
