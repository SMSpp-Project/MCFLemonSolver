# shim

The code that the build inserts into its private copies of the headers of
LEMON, written once and in C++ rather than inside the build files.

LEMON is taken from the platform, and the build rewrites private copies of
three of its headers, which it puts on the include path before the system
ones [see the *Requirements* of [../README.md](../README.md)]. Most of the
rewrites are of one token, e.g. a sign compared against a tolerance instead
of against zero, and those live in the build files, one line each; what is a
block of code lives here, and both the CMake build and the makefiles insert
it at the same anchor, so that the two agree by construction:

| file | anchor | what it is |
|---|---|---|
| `ns_eps.inc` | `namespace lemon {` | the two tolerances of `NetworkSimplex`, of the pivot rules and of the flow left on the artificial arcs, zero for an exact type so that on integer data the algorithm is the one of LEMON, and how far the update of the potentials of a pivot may move the potential of the root from zero, and the shift that brings it back [see `ns_potential.inc`] |
| `ns_runwarm.inc` | the `return start( pivot_rule );` of `run()` | the tail of `run()`, then `runWarm()` and `unbCycle()`: the warm start, i.e., the simplex goes on from the basis of the last run, whose potentials are recomputed on its tree after a change of the costs and whose tree flow is recomputed after a change of the capacities or of the supplies, a subtree whose flow no longer fits being hung from the root by its artificial arc [see `repairBasis()`], and the certificate of unboundedness, i.e., the cycle of the pivot the algorithm cannot perform when it answers `UNBOUNDED` |
| `ns_potential.inc` | the `updatePotential()` of `network_simplex.h`, which it replaces whole | the potentials of a pivot updated on the smaller of the subtree the pivot moves and the rest of the tree, i.e., the nodes that follow the subtree in the thread up to its first one, whose potentials go the other way; the reduced costs are the same, and the potential of the root, which this moves away from zero by no more than a bound set with the cost of the artificial arcs [see `ns_pi_lim()` in `ns_eps.inc`], is brought back to zero at the end of `start()` |
| `capsc_runwarm.inc` | the declaration of `_sum_supply` in `capacity_scaling.h` | `runWarm()` of `CapacityScaling`: the flow of the last run is put within the new bounds and at the one the reduced cost under the potentials of the last run says, and the shortest augmenting paths send the excess this leaves, as after any change a primal-dual algorithm allows |
| `cs_runwarm.inc` | the declaration of `_sum_supply` in `cost_scaling.h` | `runWarm()` of `CostScaling`: the flow of the last run is put within the new bounds and made to meet the new supplies [see `cs_feasible.inc`], the potentials are scaled back and the scaling starts from the largest violation of the optimality conditions under them rather than from the largest cost |
| `cc_runwarm.inc` | the declaration of `_sum_supply` in `cycle_canceling.h` | `runWarm()` of `CycleCanceling`: the flow of the last run is put within the new bounds and made to meet the new supplies [see `cs_feasible.inc`], and the negative cycles are canceled from it |
| `cs_feasible.inc` | the declaration of `_sum_supply` in `cost_scaling.h` and in `cycle_canceling.h`, after the warm start; a call to it replaces the lines of `init()` from the maps of the `Circulation` to the residual capacities set from its flow | `feasibleFlow()`: the flow within the bounds that meets the supplies, which `init()` of LEMON finds with `Circulation` on the graph and its maps, found by the same push-relabel algorithm on the residual network the algorithm works on, from the flow the arcs into the nodes with a deficit take greedily or, for `runWarm()`, from the flow of the last run put within the new bounds |

The anchors are lines of LEMON as it is released: if a new version of LEMON
moves one of them the build finds nothing to replace and the rewrite does
not happen. `markers.txt` lists, for each rewritten header, what it has to
contain once rewritten, and both builds check it and stop, naming the text
that is missing, rather than go on with a header rewritten in part: an
inserted block would make the compilation fail anyway, but a rewrite of one
token (a tolerance, a hook that keeps the warm start informed of a change)
would not, and would leave the algorithm silently wrong.

The rewrites are idempotent: applied to an already patched LEMON they are
no-ops, since each of them looks for what only the unpatched header has.
