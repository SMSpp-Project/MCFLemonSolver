# shim

The warm start of the minimum cost flow algorithms of abLEMON, which
abLEMON does not have and the build writes into private copies of five of
its headers, put on the include path before the installed ones.

`rewrites.txt` lists the rewrites, one a line as `header|anchor|what`: the
anchor is a whole line of the header, which has to be there exactly once,
and what is `<file` (the file inserted before the line), `+text` (appended
to the line) or `=text` (in place of the line). The CMake build applies
them in `CMakeLists.txt` and the makefiles with `rewrite.awk`, so that the
two agree by construction; both stop, naming the line, when an anchor is
not there exactly once, rather than go on with a header rewritten in part.

| file | what it is |
|---|---|
| `ns_runwarm.inc` | `runWarm()` of `NetworkSimplex`: the simplex goes on from the basis of the last run, whose potentials are recomputed on its tree after a change of the costs and whose tree flow is recomputed after a change of the capacities or of the supplies, a subtree whose flow no longer fits being hung from the root by its artificial arc [see `repairBasis()`] |
| `cs_runwarm.inc` | `runWarm()` of `CostScaling`: the flow of the last run is put within the new bounds and made to meet the new supplies by `feasibleFlow()` of `lemon/bits/mcf_feasible_flow.h`, which starts from it, the potentials are scaled back and the scaling starts from the largest violation of the optimality conditions under them |
| `cc_runwarm.inc` | `runWarm()` of `CycleCanceling`: the flow of the last run is put within the new bounds and made to meet the new supplies, as for `CostScaling`, and the negative cycles are canceled from it |
| `capsc_runwarm.inc` | `runWarm()` of `CapacityScaling`: the flow of the last run is put within the new bounds and at the one the reduced cost under the potentials of the last run says, and the shortest augmenting paths send the excess this leaves |

The one-line rewrites keep the warm start informed: a run records whether
it has left a basis or a flow to start from, and a change of the lower
bounds, of the supply type, of the digraph or of all the parameters drops
it; in `NetworkSimplex`, a change of an upper bound or of a supply marks
the basis to be repaired.
