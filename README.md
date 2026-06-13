# MCFLemonSolver

This project implements `MCFLemonSolver`, a SMS++ `:Solver` for
[MCFBlock](https://gitlab.com/smspp/mcfblock) based on interfacing
solvers from the [LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon).
In fact, `MCFLemonSolver` is not a single solver but up to 72 classes
obtained instantiating 4 different "base" solvers on 2 different kinds
of graphs and with all possible combinations of 3 different types of
flows and costs (`double`, `long`, `int`).

## MCFLemonSolver* variants

Upon compiling, a number of `Solver` variants is added to the `Solver`
factory. These depend on the four main algorithms implemented by the
LEMON project

- `NetworkSimplex`
- `CostScaling`
- `CycleCanceling`
- `CapacityScaling`

Each solver is instantiated on two different types of graphs:
`SmartDigraph`, which is more efficient but static (it does not allow
to add/remove, open/close arcs) and `MCFListDigraph` (a small ad-hoc
improvement of the original `ListDigraph`), which can be less efficient
but allows to add/remove, open/close arcs. These already make 8 variants,
each of which can be implemented with different of costs and capacities
on the arcs, which can be `int`, `double` or `long`. Although `int` is
supported, for large graphs (or large costs/capacities/deficits on small
graphs) it risks overflows.

Thus, up to 72 variants can be inserted in the `Solver` factory, with the
general form

    A<G,C,V>

with A chosen in

    MCFLemonSolverNetworkSimplex ,
    MCFLemonSolverCycleCanceling ,
    MCFLemonSolverCostScaling ,
    MCFLemonSolverCapacityScaling ,

G chosen in `SmartDigraph`, `MCFListDigraph`, and C, V chosen in `int`,
`double` or `long` (the first being the cost type, the second the flows
type). Since these can be many, the macro `SMSpp_which_insert_LEMON`
in `MCFLemonSolver.cpp` allows to restrict them to only a subset of the
supported cost / flow types. The value is numeric and coded bitwise:

 - bit 0 ( + 1 ): double costs and/or flows

 - bit 1 ( + 2 ): long costs and/or flows

 - bit 2 ( + 3 ): int costs and/or flows

With only one bit set to 1, the 8 variants having only that type as flows
and costs are inserted. With two bits set to 1, 8 variants are inserted for
each of the 4 possible combinations of the two types as flows and costs
(i.e., 32 variants). With all three bits set to 1, 8 variants are inserted
for each of the 9 possible combinations of the three types as flows and
costs (i.e., 72 variants. The macro can be changed in the `makefile`.

The solvers provided by the LEMON project are used through a thin interface in
`MCFLemonSolver.h`. LEMON has not had a release since 1.3.1 (2014), so its
headers need a few portability fixes:

- **C++20**: `lemon/bits/array_map.h` and `lemon/path.h` still call the
  `std::allocator::construct`/`destroy` members that were *removed* in C++20;
- **MSVC**: `lemon/adaptors.h` selects an `LEMON_SCOPE_FIX` definition that does
  not compile under MSVC;
- **missing include**: `lemon/capacity_scaling.h` uses `RangeMap` from
  `lemon/maps.h`, which some packagings (e.g. MacPorts `coinor-liblemon`) fail
  to include.

Rather than redistribute a patched copy of LEMON, the build consumes the LEMON
package provided by your platform (see *Requirements* below) and applies a tiny,
self-contained **compatibility shim** at build time: it rewrites private copies
of just those headers and puts them on the include path *before* the system ones
so they shadow them. The rewrites are idempotent (no-ops on an already-patched
LEMON) and are applied identically by both the CMake build and the makefiles, so
nothing has to be done by hand on any platform — no manual download of LEMON,
and no manual edit of its headers.

## Getting started

These instructions will let you build `MCFLemonSolver` on your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and
  its requirements.

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its
  requirements (but no actual `MCFSolver` are needed, since
  `LEMONSolver` provides its own).

- The [LEMON graph library](https://lemon.cs.elte.hu/trac/lemon), installed as a
  system package — no manual download is needed, and its C++20 incompatibility is
  handled automatically by the build (see above):
  - **Linux** (apt): `liblemon-dev`
  - **macOS** (MacPorts): `coinor-liblemon` — note that Homebrew's `lemon` is the
    unrelated LALR parser generator; the community tap
    `donn/lemon-graph/lemon-graph` also exists but is unmaintained
  - **Windows** (vcpkg): `liblemon` (already listed in the project `vcpkg.json`)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

You can also choose the following configuration options:

| Variable                   | Description                                                          | Default value |
|----------------------------|----------------------------------------------------------------------|---------------|
| `SMSpp_which_insert_LEMON` | 1, i.e., double costs<br/>2, i.e., long costs<br/>3, i.e., int costs | 3             |

You can set them with:

```sh
cmake <source-path> -D<var>=<value>
```

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(MCFLemonSolver)
target_link_libraries(<my_target> SMS++::MCFLemonSolver)
```


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Like CMake, the makefiles use the installed LEMON package directly (no manual
download) and generate the same C++20 shim. By default they look for LEMON under
the `/usr` prefix; for other locations override `LEMON_ROOT` — or, more finely,
`LEMON_INCDIR`/`LEMON_LIBDIR` — on the `make` command line or in the environment,
for example:

```sh
make -f makefile LEMON_ROOT=/opt/local
```

Each executable using `LEMONSolver` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the `test/MCF_MILP`
tester that compares different solvers for `MCFBlock` (be them specialised
or general-purpose). The makefiles in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/lemonsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](../CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Daniele Caliandro**  
  Dipartimento di Informatica  
  Università di Pisa

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

The SMS++ code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details. LEMON itself is not redistributed:
it is consumed as an external system package (available under the Boost Software
License, Version 1.0), and the build only rewrites a private copy of two of its
headers in-place for C++20 compatibility, as described above.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

