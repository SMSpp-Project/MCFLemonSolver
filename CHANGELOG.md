# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the test of the module, `MCFLemonSolver_test`: the four algorithms on the
  two graphs against a small instance whose optimal value is known after each
  of the changes a MCFBlock can undergo, plus the Solver destroyed without a
  MCFBlock, the data given only by a Modification, the deficits that sum to
  zero only up to rounding, the parameters and the DMX file

### Changed

- the DMX file of `strDMXFile` is written by the MCFBlock in the complete
  DIMACS format, on both graphs, rather than as the adjacency matrix of the
  SmartDigraph only, which had no costs, capacities or deficits

- deficits that sum to zero only up to rounding are balanced before each run,
  and deficits that do not sum to zero make the problem infeasible: LEMON
  alone reads the former as infeasible and the latter, when the demand
  exceeds the supply, as demands that may be left unmet

- the flows and the potentials are available only at an optimal solution

### Fixed

- CostScaling could read out of its own vectors, and CycleCanceling end away
  from the optimum, on fractional costs: both are given the costs scaled by a
  power of 2 and rounded, which are integer, CostScaling with an integer large
  cost type, and the value and the potentials are those of the true costs

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

