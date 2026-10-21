/*--------------------------------------------------------------------------*/
/*------------------------ File MCFLemonSolver.cppm ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.mcflemonsolver, which exports every
 * name that the headers of MCFLemonSolver declare at namespace scope,
 * together with the modules it depends on; a translation unit can then write
 *
 *     import smspp.mcflemonsolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
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
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "MCFLemonSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.mcflemonsolver;

export import smspp.mcfblock;

export namespace SMSpp_di_unipi_it {

 // MCFLemonSolver.h
 using SMSpp_di_unipi_it::MCFLemonSolver;
 using SMSpp_di_unipi_it::MCFLemonSolverCapacityScaling;
 using SMSpp_di_unipi_it::MCFLemonSolverCostScaling;
 using SMSpp_di_unipi_it::MCFLemonSolverCycleCanceling;
 using SMSpp_di_unipi_it::MCFLemonSolverNetworkSimplex;
 using SMSpp_di_unipi_it::MCFListDigraph;
 using SMSpp_di_unipi_it::SMSppCapacityScaling;
 using SMSpp_di_unipi_it::SMSppCostScaling;
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File MCFLemonSolver.cppm ----------------------*/
/*--------------------------------------------------------------------------*/
