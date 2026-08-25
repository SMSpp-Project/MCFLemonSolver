# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `MCFLemonSolver::get_Solution()`, which builds the MCFSolution out of the
  data structures of the LEMON algorithm, without writing anything into the
  MCFBlock and therefore without requiring any Variable to exist

### Fixed

- the flows and the potentials were indexed by the order in which ArcIt and
  NodeIt enumerate them, which is not the order the arcs and the nodes have
  in the MCFBlock: the solution came out permuted. They are now indexed by
  the LEMON id, which is the position of creation

- `get_dual_solution()` called `set_pi()` with the arguments the other way
  round, the node index as the potential and the potential as the node index,
  which threw "invalid node name"

