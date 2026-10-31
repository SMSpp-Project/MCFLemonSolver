# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- lemon-ng, the maintained fork of LEMON, is used as it is: the CMake build
  and the makefiles recognise it by `lemon/bits/mcf_feasible_flow.h`, which
  only the fork has, and rewrite none of its headers, which already have
  what the shim writes into the ones of LEMON; `MCFListDigraph` names the
  storage of `ListDigraphBase` with a using-declaration, which the ranges
  `nodes()` and `arcs()` of the graphs of lemon-ng would hide otherwise

- the makefiles build the C++20 module `smspp.mcflemonsolver` when
  `SMSppMOD` is not empty (e.g., `make SMSppMOD=1`), with clang++ 16 or
  later [see `lib/makefile-mod` of SMS++]

- the C++20 module `smspp.mcflemonsolver` (`src/MCFLemonSolver.cppm`), which
  exports every name that the headers declare at namespace scope and
  re-exports `smspp.mcfblock`, so that a source file can write `import
  smspp.mcflemonsolver;` in place of the `#include` of the headers; it is
  built with the CMake option `SMSpp_BUILD_MODULE`, as the target
  `SMS++::MCFLemonSolver_module`

### Changed

- the certificate of unboundedness of the network simplex is
  `unboundedCycle()`, the name it has in lemon-ng, also in the shim, which
  gives the cycle in the direction of the flow as lemon-ng does

- the pipeline builds the module also with its makefiles, after CMake, and
  runs on a change of the makefiles too [see `ci/makefile` of the umbrella]

- `SMSpp_BUILD_MODULE` no longer stops the configuration when the C++20
  module cannot be built, i.e., with a CMake older than 3.28, a generator
  other than Ninja or Visual Studio, a compiler other than Clang 16 (with
  `clang-scan-deps`) or MSVC 19.34, or a module it imports that is not
  there: the option is turned off with a warning that says why, and the
  build goes on with the headers [see `cmake/CheckModule.cmake` of SMS++]

- the licence is the GNU Lesser General Public License version 3.0

- `SMSppCapacityScaling` runs CapacityScaling with a 4-ary heap in its
  shortest path searches, rather than the binary one of the default trait,
  through the trait `SMSppCapacityScalingTraits`

- the network simplex of LEMON updates the potentials of a pivot on the
  smaller of the subtree that the pivot moves and the rest of the tree,
  rather than always on the subtree, and brings the potential of the root
  back to zero at the end [see `shim/ns_potential.inc`]

- CostScaling and CycleCanceling find the feasible flow they start from with
  a push-relabel algorithm on their own residual network, rather than with the
  `Circulation` of LEMON on the graph and its maps, and their re-optimization
  after a change of the capacities or of the supplies starts it from the flow
  of the last run put within the new bounds, rather than finding a feasible
  flow anew [see `shim/cs_feasible.inc`]

### Fixed

- `cmake/FindLEMON.cmake` keeps the version of a found library in the cache,
  with the include directory it is read from, and reads it again when that
  directory changes: the check of the cache wanted the version, which was not
  cached, so that every configure searched the library again

- `cmake/FindLEMON.cmake` drops the cached paths of LEMON when one of them
  does not exist any more (e.g., in a build directory configured once outside
  and once inside a container, or after an upgrade), so that all of them are
  searched again; the stale include directory was kept, since `find_path()`
  does not search a cached one again, and the configuration went on with a
  LEMON that is not there

- CostScaling and CycleCanceling on integer flows are given the capacities
  anew before each run: they replace an infinite capacity by the sum of the
  deficits in the vector of the capacities they are given, so that the run
  after the deficits had grown past that bound answered "infeasible"

- `get_var_solution()` writes the flows into the `MCFBlock` all at once, with
  `set_x()` of the whole vector, instead of one arc at a time: the dynamic
  arcs of the `MCFBlock` are a list, which setting the flow of one arc scans
  up to it, so that writing a solution took time quadratic in the number of
  arcs; the closed arcs, which LEMON does not enumerate, get flow 0

- the makefiles of the tests and of the tools build also with g++: they
  take the option that stops at the first error and the one of the
  debugging information, which clang++ and g++ spell differently
  (`-ferror-limit=1` and `-glldb` against `-fmax-errors=1` and `-ggdb`),
  from `SMSppSWERR` and `SMSppSWDBG` of `lib/makefile-OS` of SMS++

- FindLEMON looks in `LEMON_ROOT` before the directories of the system, and
  for both names of the library (`lemon` of the distributions, `emon` of the
  CMake build of LEMON) in a directory before the next one, so that the
  LEMON of `LEMON_ROOT` is not passed over for a `liblemon` of the system

## [0.3.0] - 2026-10-09

### Added

- the tester of this directory carries the label of the module, so that the
  pipeline, which selects with `ctest -L <module>`, runs it: it was built and
  never run, although it is what walks the four algorithms on both graphs

- `kFactor`, the scaling factor of `MCFLemonSolverCapacityScaling`, i.e., the
  base of the geometric sequence of the deltas of the successive
  approximations: it was the only algorithmic parameter of the four algorithms
  that no `ComputeConfig` could reach, the pivot rule of the network simplex
  and the method of cycle canceling and of cost scaling being there already:
  a configuration can now ask for the factor the instance at hand wants, as it
  can already ask for a pivot rule. The default is the 4 of LEMON and a value
  below 2, which `CapacityScaling` refuses, is refused here

- the network simplex re-optimizes: after a change of the costs the basis of
  the previous run is still primal feasible, so the potentials are recomputed
  on its tree and the simplex goes on from there instead of starting from the
  artificial basis, which is the change a Lagrangian or a Frank-Wolfe
  decomposition makes at each iteration; after a change of the capacities or
  of the supplies the flow of the tree is recomputed from the leaves up, and
  the subtree below an arc whose flow falls outside its bounds is hung from
  the root by its artificial arc, as `init()` hangs every node, so that the
  basis is again primal feasible and the simplex drives the artificial arcs
  out by their cost; the basis of a run that ends by proving the instance
  infeasible is kept as well, being optimal for the problem with the
  artificial arcs. LEMON has no warm start of its own, the `run()` of every
  algorithm calling its `init()`, hence `runWarm()` is added by the shim of
  the build [see shim/README.md]. On the instances of the `MCFBlock` suite a
  re-solve after a change of the costs takes between one seventh and one
  fiftieth of a solution from scratch, the larger the instance the larger the
  gain, and on a GOTO instance of 1024 nodes and 65536 arcs the instructions
  of 20 re-solves after changes of the capacities, of the deficits or of the
  closed arcs are between 15 and 42 times fewer than from scratch. The basis
  is dropped after a change of the graph, with lower bounds, and when the
  supplies do not sum to zero, the artificial arcs being built differently
  then

- cost scaling, capacity scaling and cycle canceling re-optimize, each with
  its own `kReopt`: capacity scaling, a primal-dual algorithm, starts from the
  flow and the potentials of the last run after any change of the costs, of
  the capacities or of the supplies, the flow put within the new bounds and
  at the one the reduced cost says and the excess this leaves sent by its
  shortest augmenting paths; cost scaling and cycle canceling start from the
  flow of the last run when it is still feasible, i.e., after a change of
  the costs, the former from the largest violation of the optimality
  conditions under the potentials of the last run rather than from the
  largest cost, the latter without the circulation that finds a feasible
  flow anew [see shim/README.md]

- `kReopt` for the network simplex, whether it re-optimizes from the basis of
  the previous run (1, the default) or starts from scratch at each solve (0),
  with the name the `:MCFSolver` give to the same parameter, so that the two
  can be compared against their own solve from scratch in one configuration

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

- the network simplex of LEMON declared infeasible some feasible instances
  with negative costs of a floating-point type: the cost of its artificial
  arcs is larger than that of any path of nonnegative cost, not of any path,
  and with negative costs a cycle through the root could cost less than
  zero and keep flow on the artificial arcs; the cost is now computed on the
  absolute values, both at the start and in the warm start (11 wrong answers
  out of 3000 random instances before, none after)

- the basis that the warm start of the network simplex repairs is strongly
  feasible, which the anti-cycling rule of its pivots relies upon: a tree
  arc left at the bound that would make the tree not strongly feasible (at
  its capacity if directed towards the root, at zero if away from it) leaves
  the tree as one outside its bounds does, and so does a tree arc closed at
  zero flow

- the warm start of cost scaling starts from the phase after the largest
  violation of the optimality conditions, the pair of the last run being
  already that close to optimal, rather than refining to where it already is

- the warm start of capacity scaling shifts the potentials of the last run
  so that the smallest is zero, which leaves the reduced costs as they are
  and keeps their values from drifting over a long sequence of runs

- cost scaling, cycle canceling and capacity scaling are given integer costs
  as they are, with scale 1, as they are given integer capacities and
  supplies: they were always multiplied by the largest power of 2 that the
  bounds of the arithmetic allow, which on integer data only added phases to
  cost scaling and to the cancel-and-tighten method of cycle canceling

- the build stops, naming what is missing, if a rewrite of a header of LEMON
  has not found its anchor [see shim/markers.txt], rather than going on
  with a header rewritten in part

- `get_dual_solution()` wrote the potentials only, leaving in the bound
  constraints the dual values of whatever had written there before, so that
  the dual solution read out of the abstract representation was not one;
  it now writes the reduced costs as well, as the `:MCFSolver` do

- cost scaling could read and write past the end of its buckets: in the
  price refinement the rank of a node is a sum of ranks along a path, each
  of them checked against the number of buckets but not their sum, which
  LEMON as released gets wrong on small random instances; the shim bounds
  it, the rank being a lower bound on how far the potential has to go down

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- a project built against the installed package did not compile: the
  rewritten LEMON headers were given to the build tree alone, so that
  `MCFLemonSolver.h`, which includes the LEMON headers and calls what only
  the rewritten ones have (e.g., `runWarm()`), found the system ones; they
  are now installed next to the headers of the module, in `lemon-shim`, and
  the exported target puts them before the system ones

- the tester did not compile against an installed LEMON whose
  `lemon/capacity_scaling.h` misses the include of `lemon/maps.h`, since the
  shim of those headers was given to the library alone: it now reaches
  whatever is built in the same tree against the library, as
  `MCFLemonSolver.h` includes the LEMON headers

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

[Unreleased]: https://gitlab.com/smspp/mcflemonsolver/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/mcflemonsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/mcflemonsolver/-/compare/0.1.0...0.2.0
