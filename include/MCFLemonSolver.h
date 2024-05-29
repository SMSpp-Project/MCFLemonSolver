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
 * - GR, which is the type of graph
 *  
 *   The possibilities are:
 * 
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
 * 
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
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
 * types of LEMON graphs, i.e.,
 * 
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
 * 
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
 *
 * Similarly, the concept LEMONAlgorithm is defined that only allows all
 * possible types of LEMON algorithms, i.e., [SMSpp]CapacityScaling,
 * [SMSpp]CostScaling, CycleCanceling, and NetworkSimplex.
 *
 *  @{ */

  template< typename Algo >
  struct Fields {};
   

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
 template< LEMONGraph GR , typename V , typename C >
 class SMSppCapacityScaling : public
  CapacityScaling< GR , V , C , CapacityScalingDefaultTraits< GR , V , C > >
  {};

 /// CostScaling algorithm using the default trait
 template < LEMONGraph GR , typename V , typename C>
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
 * - GR, which is the type of graph:
 * 
 *   The possibilities are:
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
 * 
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
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
 *   C) implemented in the LEMON class:
 *   The possibilities are:
 * 
 *   - NetworkSimplex implements the primal Network Simplex algorithm for finding
 *     a minimum cost flow. This algorithm is a highly efficient specialized version
 *     of the linear programming simplex method directly for the minimum cost flow
 *     problem.
 *     
 *   - CycleCanceling implements three different cycle-canceling algorithms for finding 
 *     a minimum cost flow. The most efficent one is the Cancel-and-tighten algorithm,
 *     thus it is the default method. It runs in strongly polynomial time, but in practice, 
 *     it is typically orders of magnitude slower than the scaling algorithms and NetworkSimplex.
 * 
 *   - CostScaling implements a cost scaling algorithm that performs push/augment and
 *     relabel operations for finding a minimum cost flow. It is a highly efficient primal-dual
 *     solution method, which can be viewed as the generalization of the preflow push-relabel
 *     algorithm for the maximum flow problem. It is a polynomial algorithm.
 * 
 *   - CapacityScaling implements the capacity scaling version of the successive shortest path
 *     algorithm for finding a minimum cost flow. It is an efficient dual solution method,
 *     which runs in polynomial time.
 *     In special case it can be more efficient than CostScaling and NetworkSimplex algorithms.
 * 
 *   In general, NetworkSimplex and CostScaling are the fastest implementations available in LEMON
 *   for solving this problem.
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
          typename GR , typename V , typename C >
  requires LEMONGraph< GR >
 class MCFLemonSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

  using ThisAlgo = Algo< GR , V , C >;
 
  enum str_par_type_LEMON {
  strDMXFile = strLastParCDAS ,  ///< DMX filename to output the instance
  strLastParLEMON    ///< first allowed parameter value for derived classes
                   /**< convenience value for easily allow derived classes
                    * to further extend the set of types of return codes */
  };

  enum LEMON_sol_type
  {
  UNSOLVED, //= NULL, ///< the problem has not been solved yet
  OPTIMAL,           ///< the problem has been solved
  KSTOPTIME, //= NULL,     ///< the problem has been stopped because of time limit
  INFEASIBLE,   ///< the problem is provably infeasible
  UNBOUNDED,    ///< the problem is provably unbounded
  KERROR //= NULL         ///< the problem has been stopped because of unrecoverable error
  }; 
  
  void set_Block( Block * block) override;


  protected:

  void guts_of_constructor( void ) {
   f_algo = nullptr;
   }

  void guts_of_destructor( void ) {
   delete f_algo;
   delete dgp;
   }

  ThisAlgo * f_algo;
  ///< the actual LEMON algorithm for solving the MCFBlock

  GR * dgp;
  /**< represents the directed graph implemented by two classes by Lemon
   * (ListDigraph and SmartDigraph) */

 };// End of MCFLemonSolver
          

/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- NETWORKSIMPLEX --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolver<NetworkSimplex, GR, V, C> that contains
 * specialized  compute() method, enums for indexing algorithimc parameters
 * and function  set/get_*_par for manage them.
 * 
 *  Template parameters are:
 * 
 *  - NetworkSimplex, implements the primal Network Simplex algorithm for
 *    finding a minimum cost flow. This algorithm is a highly efficient
 *    specialized version of the linear programming simplex method directly
 *    for the minimum cost flow problem.
 * 
 *  - GR that represents the directed graph, the possibilities are described 
 *    in the file MCFLemonSolver.h
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 */

template < typename GR , typename V , typename C >
class MCFLemonSolverNetworkSimplex :
 public MCFLemonSolver< NetworkSimplex , GR , V , C >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using BaseClass = MCFLemonSolver< NetworkSimplex , GR , V , C >;
 
 using BaseClass::intLastParCDAS;
 //using BaseClass::idx_type;

 using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
 using SMSpp_di_unipi_it::Solver::OFValue;
 using SMSpp_di_unipi_it::ThinComputeInterface::kUnEval;
 using SMSpp_di_unipi_it::CDASolver::kStopTime;
 using SMSpp_di_unipi_it::CDASolver::kInfeasible;
 using SMSpp_di_unipi_it::Solver::lock;
 using SMSpp_di_unipi_it::Solver::unlock;
 using SMSpp_di_unipi_it::Solver::f_Block;
 using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
 using SMSpp_di_unipi_it::Solver::f_id;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::strLastParLEMON;
 using BaseClass::str_par_type_LEMON::strDMXFile;
 using BaseClass::LEMON_sol_type::UNSOLVED;

 using typename BaseClass::ThisAlgo;
 using NSPivotRule = typename ThisAlgo::PivotRule;

/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: Initializes algorithm parameters
 /** Void constructor. Define f_pivot_rule to the default algorithm
  * parameters used by NetworkSimplex.
  */

 MCFLemonSolverNetworkSimplex( void ) {
  BaseClass::guts_of_constructor();
  f_pivot_rule = NSPivotRule::BLOCK_SEARCH;
  }
 
  /// destructor: delete algorithm parameter
  /** Does nothing special, delete Fields f_pivot_rule, an algorithmic
   * parameter of  NetworkSimplex */

  ~MCFLemonSolverNetworkSimplex() {
   BaseClass::guts_of_destructor();
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
 enum LEMON_NS_dbl_par_type{
  dblLastParLEMON_NS ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };

 enum LEMON_NS_int_par_type{
  kPivot = intLastParCDAS,
  intLastParLEMON_NS ///< first allowed parameter value for derived classes
  /**< convenience value for easily allow derived classes
   * to further extend the set of types of return codes */
  };
 
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
 *  @{ */

 static constexpr int kErrorStatus = -1;

/*--------------------------------------------------------------------------*/
 // set the ostream for the Solver log
 // not really, MCFClass objects are remarkably silent
 //
 // virtual void set_log( std::ostream *log_stream = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// @brief set the parameter par with value
 /// @param par 
 /// @param value 

 void set_par( idx_type par , int value ) override {      
  if( par == kPivot ) {      
   if( ( value < 0 ) || ( value > 4 ) )
    throw( std::invalid_argument( "Error: invalid kPivot " +
				  std::to_string( value ) ) );

   if( value == f_pivot_rule )
    return;  // nothing is changed
               
   f_pivot_rule = NSPivotRule( value );
   return;
   }

  CDASolver::set_par( par , value );
  return;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the MCF encoded by the current MCFBlock
 *  @{ */

 /// (try to) solve the MCF encoded in the MCFBlock 

 int compute(bool changedvars = true) override;
/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */
    
 int get_status( void ) const { return( this->status ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 ///Get elapsed time for run() method

 double get_elapsed_time( void ) const override {
  return( this->ticks );
  }

/*--------------------------------------------------------------------------*/
 //Return the lower bound solution(optimal) for the problem

 OFValue get_lb( void ) override { return( OFValue( f_algo->totalCost() ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 //Return the upper bound solution(optimal) for the problem

 OFValue get_ub( void ) override { return( OFValue( f_algo->totalCost() ) ); }

    /*--------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/

 ///@return number of int algorithmic parameters
 [[nodiscard]] idx_type get_num_int_par( void ) const override
 {
  return( intLastParLEMON_NS );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
 /// @return number of dbl algorithmic parameters
 [[nodiscard]] idx_type get_num_dbl_par(void) const override
 {
  return (dblLastParLEMON_NS);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// MOVE THIS INTO BASE CLASS
 /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (strLastParLEMON);
    }
    
/*--------------------------------------------------------------------------*/
 /// @brief used for obtain default value of an int parameter
 /// @param par 
 /// @return default value of parameter par, if exists
 [[nodiscard]] int get_dflt_int_par(idx_type par) const override
 {
  if( par > intLastParLEMON_NS)
   throw std::invalid_argument("Invalid int parameter: out_of_range "
			       + std::to_string(par));

  if( par == kPivot )
   return( NSPivotRule::BLOCK_SEARCH );

  return( CDASolver::get_dflt_int_par( par ) );
  }
    
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @brief used for obtain default value of a double parameter
 /// @param par 
 /// @return default value of parameter par, if exists
 /*
 [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if( par > dblLastParLEMON_NS){
        throw std::invalid_argument("Invalid dbl parameter: out_of_range " + std::to_string(par));
      }
      switch(par){
        default: return(CDASolver::get_dflt_dbl_par(par)); 
      }   
    }
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// MOVE THIS TO BASE CLASS
 /// @brief used for obtain default value of a string parameter
 /// @param par
 /// @return default valule of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {
      if(par > strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }
      static const std::string _empty;
      if (par == strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
        
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// MOVE THIS TO BASE CLASS
 /// @brief used for get the value of string parameters
 /// @param par 
 /// @return value of parameter indexed by par
 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {       
  if( par == strDMXFile )
   return( this->f_dmx_file );

  return( get_dflt_str_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 /// @brief used for get the value of int parameters
 /// @param par 
 /// @return value of parameter indexed by par

 [[nodiscard]] int get_int_par( idx_type par ) const override
 {
  if( par == kPivot )
   return( f_pivot_rule );

  return( get_dflt_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @brief used for get the value of dbl parameters
 /// @param par 
 /// @return value of parameter indexed by par
 /*
 [[nodiscard]] double get_dbl_par( idx_type par ) const override
 {
  //Da finire parametri algoritmici dbl

  return( get_dflt_dbl_par( par ) );
  }
 */
/*--------------------------------------------------------------------------*/
 /// @brief used for convert string to int parameter's index
 /// @param name 
 /// @return an index that denotes parameter name

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
 const override {
  if( name == "kPivot" )
   return( kPivot );

  return( CDASolver::int_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @brief used for convert string to dbl parameter's index
 /// @param name 
 /// @return an index that denotes parameter name

 /* not useful so far
 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
 const override {
  return (CDASolver::dbl_par_str2idx(name));
  }
 */

/*--------------------------------------------------------------------------*/
 /// @brief used for convert int index idx to phrasal rapresentation
 /// @param idx 
 /// @return a string that denotes index idx parameter

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
    
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// @brief used for convert dbl index idx to phrasal rapresentation
 /// @param idx 
 /// @return a string that denotes index idx parameter

 /*
 [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_NS){
        throw std::invalid_argument("Invalid index: out_of_range " + std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
 */

/** @} ---------------------------------------------------------------------*/
/*-------- METHODS FOR HANDLING THE State OF THE MCFLemonSolver -----------*/
/*--------------------------------------------------------------------------*//** @name Handling the State of the MCFSolver
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

 // friend class MCFLemonState; // make MCFSolverState friend

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

 std::string f_dmx_file; 
 ///< string for DMX file output

 NSPivotRule f_pivot_rule;

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

 int status = UNSOLVED;
 // Variable used in compute function for getting status

 typename NetworkSimplex< GR , V , C >::ProblemType status_2_pType;
 // Status of compute() method

 double ticks;  //Elaped time in ticks for compute() method

/*--------------------------------------------------------------------------*/

 }; // end( class MCFLemonSolver<NetworkSimplex, GR, V, C>  Specialization)


/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CYCLECANCELING --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CycleCanceling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 * 
 *  Template parameters are:
 * 
 *  -  CycleCanceling implements three different cycle-canceling algorithms for finding 
 *     a minimum cost flow. The most efficent one is the Cancel-and-tighten algorithm,
 *     thus it is the default method. It runs in strongly polynomial time, but in practice, 
 *     it is typically orders of magnitude slower than the scaling algorithms and NetworkSimplex.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/
  template< typename GR,  typename V, typename C>
class MCFLemonSolverCycleCanceling:  
public MCFLemonSolver<CycleCanceling, GR, V, C>
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 using BaseClass = MCFLemonSolver< CycleCanceling , GR , V , C >;
 
 using BaseClass::intLastParCDAS;
 //using BaseClass::idx_type;

 using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
 using SMSpp_di_unipi_it::Solver::OFValue;
 using SMSpp_di_unipi_it::BoxConstraint::kUnEval;
 using SMSpp_di_unipi_it::CDASolver::kStopTime;
 using SMSpp_di_unipi_it::CDASolver::kInfeasible;
 using SMSpp_di_unipi_it::CDASolver::dblLastParCDAS;
 using SMSpp_di_unipi_it::Solver::lock;
 using SMSpp_di_unipi_it::Solver::unlock;
 using SMSpp_di_unipi_it::Constraint::f_Block;
 using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
 using SMSpp_di_unipi_it::Solver::f_id;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::strLastParLEMON;
 using BaseClass::str_par_type_LEMON::strDMXFile;
 using BaseClass::LEMON_sol_type::UNSOLVED;

 using typename BaseClass::ThisAlgo;
 using CCMethod = typename ThisAlgo::Method;
 


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

 /// constructor: assign the default parameter
 /** Void constructor: Build f_method with default parameter */

 MCFLemonSolverCycleCanceling( void ) : MCFLemonSolver<CycleCanceling, GR, V, C>() {
  BaseClass::guts_of_constructor();
  f_method = CycleCanceling<GR, V, C>::Method::CANCEL_AND_TIGHTEN;
  }
 
 ///destructor:
 /**Void destructor: Delete f_method from memory */
  ~MCFLemonSolverCycleCanceling( void ) {
    BaseClass::guts_of_destructor();
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

 enum LEMON_CC_int_par_type{
        kMethod = intLastParCDAS,
        intLastParLEMON_CC ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };

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

/// public enum for the type of the solution
/// to MCFLemonSolver< CycleCanceling , C , V >

enum dbl_par_type_LEMON_CC{
  dblCycleCancelingFactor = dblLastParCDAS ,
  dblLastParLEMON_CC ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
};

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

 

    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
    /// @brief set the parameter par with value
    /// @param par 
    /// @param value 
    void set_par(idx_type par, int value){

        if( par == kMethod ) {      
        if( ( value < 0 ) || ( value > 4 ) )
          throw( std::invalid_argument( "Error: invalid kPivot " +
                std::to_string( value ) ) );

        if( value == f_method )
          return;  // nothing is changed
                    
        f_method = CCMethod( value );
        return;
        }
        

        CDASolver::set_par( par , value );
        return;
    } 
 
    
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
        int compute( bool changedvars = true) override;
        

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    ///Return the status of the run() in compute method
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    //Return the elapsed time for the run() in compute method   
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
    
    /// @brief Analyze the return status of run() 
    /// @param  
    /// @return if the solution found is OPTIMAL or UNBOUNDED returns true else false
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( CycleCanceling< GR , int , int >::ProblemType::OPTIMAL ):
      case( CycleCanceling< GR , int , int >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief Analyze the return status of run()
    /// @param  
    /// @return if the solution found is OPTIMAL or INFEASIBLE returns true else false
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( CycleCanceling< GR , int , int >::ProblemType::OPTIMAL ):
      case( CycleCanceling< GR , int , int >::ProblemType::INFEASIBLE ):
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


     ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    {
        return (intLastParLEMON_CC);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return (dblLastParLEMON_CC);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CC){
        throw(std::invalid_argument(std::to_string(par)));
      }

      switch(par){
        case kMethod: return CycleCanceling<GR, int, int>::Method::CANCEL_AND_TIGHTEN;
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of a double parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if(par > dblLastParLEMON_CC){
        throw(std::invalid_argument(std::to_string(par)));
      }

      switch(par){
        default: return (CDASolver::get_dflt_dbl_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of a string parameter
    /// @param par
    /// @return default valule of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

       if(par > strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {
      
      if(par == kMethod){
        return f_method;
      }

      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      //Da finire parametri algoritmici dbl
      return (get_dflt_dbl_par(par));
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
      if( par == strDMXFile )
        return( this->f_dmx_file );

      return( get_dflt_str_par( par ) );
    }
   

    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      if (name == "kMethod")
        return (kMethod);

      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {

      if(idx > intLastParLEMON_CC){
        throw std::invalid_argument(std::to_string(idx));
      }
      static const std::string par = "kMethod";
      switch(idx){
        case kMethod: return (par);
        default: break;
      }


      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CC){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
    
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
    CCMethod f_method;

    int status = UNSOLVED;  //Variable used in compute function for getting status
    typename CycleCanceling< GR , int , int >::ProblemType status_2_pType;
    //Status of compute() method
    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<CycleCanceling> specialization )

  /*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CAPACITYSCALING -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CapacityScaling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 * 
 *  Template parameters are:
 * 
 *   - CapacityScaling implements the capacity scaling version of the successive shortest path
 *     algorithm for finding a minimum cost flow. It is an efficient dual solution method,
 *     which runs in polynomial time.
 *     In special case it can be more efficient than CostScaling and NetworkSimplex algorithms.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/
template< typename GR, typename V, typename C>
class MCFLemonSolverCapacityScaling: 
public MCFLemonSolver<SMSppCapacityScaling, GR, V, C>
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


 using BaseClass = MCFLemonSolver< SMSppCapacityScaling , GR , V , C >;
 
 using BaseClass::intLastParCDAS;
 //using BaseClass::idx_type;

 using SMSpp_di_unipi_it::BoxConstraint::idx_type;
 using SMSpp_di_unipi_it::Solver::OFValue;
 using SMSpp_di_unipi_it::BoxConstraint::kUnEval;
 using SMSpp_di_unipi_it::CDASolver::kStopTime;
 using SMSpp_di_unipi_it::CDASolver::kInfeasible;
 using SMSpp_di_unipi_it::CDASolver::dblLastParCDAS;
 using SMSpp_di_unipi_it::Solver::lock;
 using SMSpp_di_unipi_it::Solver::unlock;
 using SMSpp_di_unipi_it::Constraint::f_Block;
 using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
 using SMSpp_di_unipi_it::Solver::f_id;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::strLastParLEMON;
 using BaseClass::str_par_type_LEMON::strDMXFile;
 using BaseClass::LEMON_sol_type::UNSOLVED;

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

 /// constructor: does nothing
 MCFLemonSolverCapacityScaling( void ) : BaseClass() {
 BaseClass::guts_of_constructor();
 }
 
  /// destructor: doesnothing  
  ~MCFLemonSolverCapacityScaling( void ) {   
    BaseClass::guts_of_destructor();
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

  enum LEMON_CS_int_par_type{
    intLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
  };


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

enum LEMON_CS_dbl_par_type{
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
};



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



    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
  
    void set_par(idx_type par, int value){

        CDASolver::set_par(par,value);
        return;
        
    }
   
    
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
        int compute( bool changedvars = true) override;
        
    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    /// @brief getter for status field
    /// @param  
    /// @return value of this->status
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    ///Get elapsed time in run() method
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
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::UNBOUNDED ):  
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
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::INFEASIBLE ):
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

     ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    {
        return(intLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return(dblLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }

      switch(par){
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an dbl parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
       if(par > dblLastParLEMON_CS){
          throw std::invalid_argument(std::to_string(par));
       }    
       
       switch(par){
        default: return (CDASolver::get_dflt_dbl_par(par));
       }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an string parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

      if(par > strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {
      //No int parameters
      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      //No dbl parameters
      return (get_dflt_dbl_par(par));
    }
        
        
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
        if( par == strDMXFile )
        return( this->f_dmx_file );

        return( get_dflt_str_par( par ) );
        }
   

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {
      if(idx > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }

    
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
    
    int status = UNSOLVED;  //Variable used in compute function for getting status
    typename SMSppCapacityScaling< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    
    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<SMSppCapacityScaling, GR, V, C>  Specialization)



/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- COSTSCALING -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CostScaling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 * 
 *  Template parameters are:
 * 
 *   - CostScaling implements a cost scaling algorithm that performs push/augment and
 *     relabel operations for finding a minimum cost flow. It is a highly efficient primal-dual
 *     solution method, which can be viewed as the generalization of the preflow push-relabel
 *     algorithm for the maximum flow problem. It is a polynomial algorithm.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/

template< typename GR, typename V, typename C>
class MCFLemonSolverCostScaling: 
public MCFLemonSolver<SMSppCostScaling, GR, V, C>
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

 using BaseClass = MCFLemonSolver< SMSppCostScaling , GR , V , C >;
 
 using BaseClass::intLastParCDAS;
 //using BaseClass::idx_type;

 using SMSpp_di_unipi_it::BoxConstraint::idx_type;
 using SMSpp_di_unipi_it::Solver::OFValue;
 using SMSpp_di_unipi_it::BoxConstraint::kUnEval;
 using SMSpp_di_unipi_it::CDASolver::kStopTime;
 using SMSpp_di_unipi_it::CDASolver::kInfeasible;
 using SMSpp_di_unipi_it::Solver::lock;
 using SMSpp_di_unipi_it::Solver::unlock;
 using SMSpp_di_unipi_it::Constraint::f_Block;
 using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
 using SMSpp_di_unipi_it::Solver::f_id;
 using BaseClass::dgp;
 using BaseClass::f_algo;
 using BaseClass::strLastParLEMON;
 using BaseClass::str_par_type_LEMON::strDMXFile;
 using BaseClass::LEMON_sol_type::UNSOLVED;

 using typename BaseClass::ThisAlgo;
 using CSMethod = typename ThisAlgo::Method;
   

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

 /// constructor: Build f_method
 /** Void constructor: Initialize f_method to default algorithmic parameter */

 MCFLemonSolverCostScaling( void ) {
  BaseClass::guts_of_constructor();
  f_method = SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;


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
  
  /// @brief delete f_method algorithmic parameter pointer
  /// @param  
  ~MCFLemonSolverCostScaling( void ) {
    BaseClass::guts_of_constructor();
    delete Fields<SMSppCostScaling<GR, V, C>>::f_method;    
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

   enum LEMON_CS_int_par_type{
        kMethod = intLastParCDAS,
        intLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
   };

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
 enum LEMON_CS_dbl_par_type{
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };

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



    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
    
    /// @brief set the parameter par with value
    /// @param par 
    /// @param value 
    void set_par(idx_type par, int value){

        if( par == kMethod ) {      
        if( ( value < 0 ) || ( value > 4 ) )
          throw( std::invalid_argument( "Error: invalid kPivot " +
                std::to_string( value ) ) );

        if( value == f_method )
          return;  // nothing is changed
                    
        f_method = CSMethod( value );
        return;
        }
        

        CDASolver::set_par( par , value );
        return;
    } 
    
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
    int compute( bool changedvars = true ) override;
        
    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    /// @brief getter of the status
    /// @param  
    /// @return this->status
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief getter of ticks (elapsed time of run())
    ///@return this->ticks
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
      case( SMSppCostScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCostScaling< GR , V , C >::ProblemType::UNBOUNDED ):  
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
      case( SMSppCostScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCostScaling< GR , V , C >::ProblemType::INFEASIBLE ):
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

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    { 
       return (intLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    ///@return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return (dblLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    ///@return number of str algorithmic parameters
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {

      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }

      switch(par){
        case kMethod: return SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an dbl parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }
      
      switch(par){
        default: return(CDASolver::get_dflt_dbl_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an string parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

      if(par > strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {

      if(par == kMethod){
        return f_method;
      }
      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      return (get_dflt_dbl_par(par));
    }
        
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
        if( par == strDMXFile )
        return( this->f_dmx_file );

        return( get_dflt_str_par( par ) );
    }  
    
    
    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      if (name == "kMethod")
        return (kMethod);

      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return ( name == "kMethod" ? kMethod : CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {
      if(idx > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }
      
      static const std::string par = "kMethod";
      switch(idx){
        case kMethod: return (par);
        default: break;
      }

      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
    
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
    CSMethod f_method;
    int status = UNSOLVED;  //Variable used in compute function for getting status
    typename SMSppCostScaling< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    double ticks;  //Elapsed time in ticks for compute() method

    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<SMSppCostScaling, GR, V, C>  Specialization)


/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFLemonState ---------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe the "internal state" of a MCFLemonSolver
/** Derived class from State to describe the "internal state" of a
 * MCFLemonSolver.  */

  class MCFLemonState : public State
  {
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

  public:
    /*------------- CONSTRUCTING AND DESTRUCTING MCFLemonState ----------------*/

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

    //virtual ~MCFLemonState() { delete f_state; }

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
