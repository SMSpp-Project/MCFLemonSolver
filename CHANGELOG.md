# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `kFactor`, the scaling factor of `MCFLemonSolverCapacityScaling`, i.e., the
  base of the geometric sequence of the deltas of the successive
  approximations: it was the only algorithmic parameter of the four algorithms
  that no `ComputeConfig` could reach, the pivot rule of the network simplex
  and the method of cycle canceling and of cost scaling being there already:
  a configuration can now ask for the factor the instance at hand wants, as it
  can already ask for a pivot rule. The default is the 4 of LEMON and a value
  below 2, which `CapacityScaling` refuses, is refused here

- the network simplex re-optimizes: after a change of the costs alone the
  basis of the previous run is still primal feasible, so the potentials are
  recomputed on its tree and the simplex goes on from there instead of
  starting from the artificial basis, which is the change a Lagrangian or a
  Frank-Wolfe decomposition makes at each iteration. LEMON has no warm start
  of its own, the `run()` of every algorithm calling its `init()`, hence
  `runWarm()` is added by the shim of the build [see shim/README.md]. On the
  instances of the `MCFBlock` suite a re-solve after a change of the costs
  takes between one seventh and one fiftieth of a solution from scratch, the
  larger the instance the larger the gain; a change of the capacities, of the
  supplies or of the graph drops the basis, and there the cost is the one it
  was

- the network simplex re-optimizes after a change of the capacities as well:
  the basis of the last run stays primal feasible if every arc whose capacity
  changes keeps the flow it holds within the new bound and, when it is out of
  the basis tree, sits at its lower bound, which is the case of an arc that is
  closed while it carries no flow and of one that is opened again, i.e. of
  what a decomposition that opens and closes arcs does at each iteration; an
  arc that sits at its upper bound, or one whose flow no longer fits, drops
  the basis as before. On the instances of the `MCFBlock` suite a re-solve
  after closing and re-opening idle arcs takes between one sixteenth and one
  three-hundredth of a solution from scratch

- `has_var_direction()`, `get_var_direction()` and a `Solution` that says it
  holds a direction for the network simplex: the cycle of negative cost and
  infinite capacity that proves the instance unbounded is the one of the
  pivot the algorithm cannot perform when it answers `UNBOUNDED`, which the
  shim reads out of it [see `unbCycle()`]; the other three algorithms give no
  certificate and say so

- the test of the module, `MCFLemonSolver_test`: the four algorithms on the
  two graphs against a small instance whose optimal value is known after each
  of the changes a MCFBlock can undergo, plus the Solver destroyed without a
  MCFBlock, the data given only by a Modification, the deficits that sum to
  zero only up to rounding, the parameters and the DMX file

### Changed

- an arc the `MCFBlock` closes is given zero capacity instead of being taken
  out of the graph, and one it re-opens is given back the capacity the
  `MCFBlock` has: the graph is left alone, hence the algorithm need not be
  reset, which it has to be whenever the graph changes and which throws away
  everything it knows, and `SmartDigraph`, which cannot take arcs out at all,
  now follows the closing and the re-opening of arcs as well. On the flow
  relaxation of the capacitated facility location, where every round of the
  test closes and re-opens facilities, the battery takes 1464 seconds against
  the 1849 it took, i.e., a fifth less; where few arcs of many are closed
  there is nothing to gain, the closed ones being scanned by the algorithm
  anyway

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole
- the DMX file of `strDMXFile` is written by the MCFBlock in the complete
  DIMACS format, on both graphs, rather than as the adjacency matrix of the
  SmartDigraph only, which had no costs, capacities or deficits

- deficits that sum to zero only up to rounding are balanced before each run,
  and deficits that do not sum to zero make the problem infeasible: LEMON
  alone reads the former as infeasible and the latter, when the demand
  exceeds the supply, as demands that may be left unmet

- the flows and the potentials are available only at an optimal solution

### Fixed

- a change of an arc that is removed from the `MCFBlock` before the Solver is
  asked to solve again made it throw "invalid arc name": the two Modification
  reach it in one batch and the first one names an arc that is no longer
  there, so the loops that go by arc name now stop at the arcs the `MCFBlock`
  has, leaving the graph to the Modification that removes the arc

- the makefiles of the module put the library of the core SMS++ before the
  objects that use it, so that a tester linking against `libSMS++.a` was left
  with unresolved symbols; it goes last, as in every other module

- CostScaling, CycleCanceling and CapacityScaling could report as infeasible,
  or not terminate on, fractional capacities or supplies: these are given to
  them scaled by a power of 2 and rounded, which are integer, and the flows
  are brought back to the true scale, each capacity being given as at most a
  bound that no optimal flow of least total flow exceeds, so that a large
  capacity no optimal flow uses does not set the rounding of the other data.
  CostScaling runs by default the `AUGMENT` method on scaled flows and the
  `PARTIAL_AUGMENT` method, which is much slower on them, on the others,
  through the new value 3 of `kMethod`

- CostScaling, CycleCanceling and CapacityScaling answered "unbounded" on a
  negative cost of infinite capacity even when the optimum is finite: the
  infinite capacity is given to them as a finite bound that no optimal flow
  reaches, and a flow that reaches it makes the problem unbounded

- CostScaling could read out of its own vectors, CycleCanceling end away from
  the optimum, and CapacityScaling report as infeasible a problem that is not,
  on fractional costs: the three of them are given the costs scaled by a power
  of 2 and rounded, which are integer, CostScaling with an integer large cost
  type, and the value and the potentials are those of the true costs

- the network simplex of LEMON could pivot for ever on fractional costs, a
  reduced cost that is zero being computed a few ulp below it, and could report
  as infeasible a problem that is not, a few ulp of flow being left on its
  artificial arcs: the shim of the build compiles it with tolerances on both,
  zero for integer types. The multicommodity Lagrangian dual, whose costs are
  fractional, and the flow relaxation of the capacitated facility location,
  whose supplies are, are now solved by it

- adding an arc to the MCFListDigraph did not grow the maps of the costs and
  of the capacities, which were then written past their end: the heap was
  corrupted as soon as an arc was added, and the run crashed later

- a Solver attached to a MCFBlock that had closed or deleted arcs built its
  graph with those arcs open, and with the NaN cost of the deleted ones

- the capacities, costs and deficits missing from the MCFBlock left their map
  unallocated, so a Modification giving them, or any change of the arcs,
  dereferenced a null pointer; each map now holds what the empty vector means

- re-attaching the Solver, or reloading the MCFBlock, leaked the previous
  graph, algorithm and maps; a Solver never attached to a MCFBlock destroyed
  an uninitialised graph; MCFLemonSolverCostScaling never destroyed its
  algorithm and graph

- compute() left the Solver locked when there was no MCFBlock or its lock
  could not be acquired

- get_lb() and get_ub() returned the cost of whatever flow LEMON left when
  the problem was infeasible or unbounded, instead of +Inf and -Inf

- the string parameter `strDMXFile` could not be set, since set_par() for
  strings was not implemented, nor found by name, since its lookup
  overrode `dbl_par_str2idx()` and so broke the lookup of every double
  parameter by name

- `kMethod` of CycleCanceling and CostScaling accepted five values out of
  the three the algorithms have

- the messages of the exceptions name the class and the method

## [0.2.0] - 2026-09-12

### Added

- `MCFLemonSolver::get_Solution()`, which builds the MCFSolution out of the
  data structures of the LEMON algorithm, without writing anything into the
  MCFBlock and therefore without requiring any Variable to exist

### Changed

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- the flows and the potentials were indexed by the order in which ArcIt and
  NodeIt enumerate them, which is not the order the arcs and the nodes have
  in the MCFBlock: the solution came out permuted. They are now indexed by
  the LEMON id, which is the position of creation

- `get_dual_solution()` called `set_pi()` with the arguments the other way
  round, the node index as the potential and the potential as the node index,
  which threw "invalid node name"

- the package configuration file finds the libraries the module links, so that
  a project using the installed module needs nothing more than find_package(),
  LEMON comprised

[Unreleased]: https://gitlab.com/smspp/mcflemonsolver/-/compare/0.2.0...develop
[0.2.0]: https://gitlab.com/smspp/mcflemonsolver/-/compare/0.1.0...0.2.0

