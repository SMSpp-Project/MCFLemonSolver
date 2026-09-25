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
| `ns_eps.inc` | `namespace lemon {` | the two tolerances of `NetworkSimplex`, of the pivot rules and of the flow left on the artificial arcs, zero for an exact type so that on integer data the algorithm is the one of LEMON |
| `ns_runwarm.inc` | the `return start( pivot_rule );` of `run()` | the tail of `run()`, then `runWarm()` and `unbCycle()`: the warm start, i.e., after a change of the costs alone the basis of the last run is still primal feasible, so the potentials are recomputed on its tree and the simplex goes on from there, and the certificate of unboundedness, i.e., the cycle of the pivot the algorithm cannot perform when it answers `UNBOUNDED` |

The anchors are lines of LEMON as it is released: if a new version of LEMON
moves them the build finds nothing to replace, the rewrite does not happen
and the compilation fails on the name that is missing, which is the loud
failure one wants rather than a silent one.

The rewrites are idempotent: applied to an already patched LEMON they are
no-ops, since each of them looks for what only the unpatched header has.
