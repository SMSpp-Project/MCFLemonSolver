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

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register the various LEMONSolver< Alg , GR , V , C > to the Solver factory

SMSpp_insert_in_factory_cpp_0_t(
 MCFLemonSolver< NetworkSimplex , SmartDigraph , double , double >);

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register MCFSolverState to the State factory

// SMSpp_insert_in_factory_cpp_0( MCFSolverState );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the various static maps

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

template<GR, V, C >
int MCFLemonSolver<NetworkSimplex<GR, V, C>, GR, V, C>::compute( bool changedvars = true ) override
{
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

      // then (try to) solve the MCF

      auto start = chrono::system_clock::now();

      if(f_pivot_rule != NULL){
          this->status = f_algo->run(f_pivot_rule);
      }else{
        //Build f_pivot_rule and execute run() method
      }
       
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
      //TODO: change MCFC function to Algo function.
      //Da rivedere, this->get_status() potrebbe non essere corretta.
      return (MCFstatus_2_sol_type[this->get_status()]);
}

    /// @brief partial specialization for 1st parameter CapacityScaling 
    /// @param changedvars 
    /// @return type of solution
        template<GR, V, C>
        int MCFLemonSolver<SMSppCapacityScaling>::compute( bool changedvars = true ) override
        {
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
        //TODO: check if factor parameter needs to be changed.
        this->status = f_algo->run();
        
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
        //TODO: change MCFC function to Algo function.
        //Da rivedere, this->get_status() potrebbe non essere corretta.
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    /// @brief partial specialization of 1st parameter for CostScaling
    /// @param changedvars 
    /// @return type of solution
        template<GR, V, C>    
        int MCFLemonSolver<SMSppCostScaling>::compute( bool changedvars = true ) override
        {.
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
        
        if(f_method != NULL){
        this->status = f_algo->run(f_method);
        }else{
                //Build f_method and execute run() method
        }
        
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


    /// @brief partial specialization of 1st parameter CycleCanceling
    /// @tparam GR 
    /// @tparam C 
    /// @tparam V   
    /// @param changedvars 
    /// @return type of solution
        template<GR, V, C>
        int MCFLemonSolver<CycleCanceling>::compute( bool changedvars = true ) override
        {
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
        
        if(f_method != NULL){
        this->status = f_algo->run(f_method);
        }else{
                //Build f_method and execute run() method
        }
        
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
 