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

  template< typename Algo >
  struct Fields {};
   
  /*template<typename GR, typename V, typename C, typename TR>
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
  };*/

template< template< typename , typename , typename > typename Algo ,
          typename GR , typename V , typename C >
  requires LEMONGraph< GR >
class MCFLemonSolver: public CDASolver, Fields< Algo<GR, V, C> > {

  int compute(bool changedvars = true) override;
};
          

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
