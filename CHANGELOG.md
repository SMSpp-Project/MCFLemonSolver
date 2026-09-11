# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

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

