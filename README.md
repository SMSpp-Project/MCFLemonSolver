# MCFLemonSolver

[![Pipeline](https://gitlab.com/smspp/mcflemonsolver/badges/develop/pipeline.svg)](https://gitlab.com/smspp/mcflemonsolver/-/pipelines?ref=develop)
[![Release](https://gitlab.com/smspp/mcflemonsolver/-/badges/release.svg)](https://gitlab.com/smspp/mcflemonsolver/-/releases)
[![Documentation](https://img.shields.io/badge/docs-Doxygen-blue)](https://smspp.gitlab.io/smspp-project/)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

This project implements `MCFLemonSolver`, a SMS++ `:Solver` for
[MCFBlock](https://gitlab.com/smspp/mcfblock) based on interfacing
solvers from the [LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon), as
maintained in [abLEMON](https://github.com/dmeoli/abLEMON).
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
to add/remove arcs) and `MCFListDigraph` (a small ad-hoc
improvement of the original `ListDigraph`), which can be less efficient
but allows to add/remove arcs. Closing and re-opening an arc, which the
`MCFBlock` does by fixing and un-fixing its flow `Variable`, both graphs
take: a closed arc is one of zero capacity, hence the graph is left alone
and the algorithm need not be reset. These already make 8 variants,
each of which can be instantiated with different types of flows and costs
on the arcs, which can be `int`, `double` or `long`. Although `int` is
supported, for large graphs (or large costs/capacities/deficits on small
graphs) it risks overflows.

Thus, up to 72 variants can be inserted in the `Solver` factory, with the
general form

    A<G,V,C>

with A chosen in

    MCFLemonSolverNetworkSimplex ,
    MCFLemonSolverCycleCanceling ,
    MCFLemonSolverCostScaling ,
    MCFLemonSolverCapacityScaling ,

G chosen in `SmartDigraph`, `MCFListDigraph`, and V, C chosen in `int`,
`double` or `long` (the first being the flow type, the second the cost
type). Since these can be many, the macro `SMSpp_which_insert_LEMON`
in `MCFLemonSolver.cpp` allows to restrict them to only a subset of the
supported cost / flow types. The value is numeric and coded bitwise:

 - bit 0 ( + 1 ): double costs and/or flows

 - bit 1 ( + 2 ): long costs and/or flows

 - bit 2 ( + 4 ): int costs and/or flows

With only one bit set to 1, the 8 variants having only that type as flows
and costs are inserted. With two bits set to 1, 8 variants are inserted for
each of the 4 possible combinations of the two types as flows and costs
(i.e., 32 variants). With all three bits set to 1, 8 variants are inserted
for each of the 9 possible combinations of the three types as flows and
costs (i.e., 72 variants). The macro is set by the CMake variable of the
same name, or by `which_insert_LEMON` in the `makefile`.

The solvers of LEMON are used through a thin interface in `MCFLemonSolver.h`,
from [abLEMON](https://github.com/dmeoli/abLEMON), the maintained fork of
LEMON, whose headers compile with the compilers and standards of today and
whose algorithms take floating-point data. What abLEMON does not have is the
re-optimization: the `run()` of every algorithm starts from scratch, and
nothing of the previous call is reused. The build writes `runWarm()` into
private copies of the headers of the four algorithms, put on the include path
before the installed ones: `NetworkSimplex` re-optimizes from the basis of the
last run, whose potentials it recomputes on its tree after a change of the
costs and whose tree flow it recomputes after a change of the capacities or of
the supplies, hanging from the root by its artificial arc the subtree below an
arc whose flow falls outside its bounds, a change of the graph dropping that
basis; `CapacityScaling` starts from the flow and the potentials of the last
run after any change, `CostScaling` and `CycleCanceling` from the flow of the
last run, which after a change of the capacities or of the supplies is put
within the new bounds and made to meet the new supplies by sending on only
what the change displaces. The code lives in [shim](shim), with the rules that
both the CMake build and the makefiles apply.

CostScaling, CycleCanceling and CapacityScaling are exact only on integer costs
too, and for them no header is rewritten: `MCFLemonSolver` gives them the
costs scaled by a power of 2 and rounded, which are integer, with an integer
large cost type for CostScaling, and gives back the value and the potentials
of the true costs; CapacityScaling is also given a 4-ary heap in its
shortest path searches, rather than the binary one of its default trait.
In the same way CostScaling, CycleCanceling and CapacityScaling, which are
exact only on integer capacities and supplies, are given fractional ones
scaled and rounded, and an infinite capacity as a finite bound that no optimal
flow reaches, since the three of them saturate the arcs of negative cost. By
default CostScaling runs its `AUGMENT` method on scaled flows, on which its
`PARTIAL_AUGMENT` method is much slower, and the latter on the others.

## Getting started

These instructions will let you build `MCFLemonSolver` on your system.

The module also comes ready-made: on Ubuntu as `libsmspp-mcflemon-dev`, in
vcpkg as `smspp[core,mcflemon]`, and in the other channels that the
[installation
guide](https://gitlab.com/smspp/smspp-project/-/wikis/Installing-SMS++) of
SMS++ lists. What follows is about building it yourself.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and
  its requirements.

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its
  requirements (but no actual `MCFSolver` are needed, since
  `MCFLemonSolver` provides its own).

- [abLEMON](https://github.com/dmeoli/abLEMON) 1.7.0 or later, the maintained
  fork of the [LEMON graph library](https://lemon.cs.elte.hu/trac/lemon), as a
  package (`ablemon` in conda-forge, vcpkg, Conan and the Homebrew tap
  `dmeoli/ablemon`) or built from its sources; LEMON itself is not enough,
  and the build says so.

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
| `SMSpp_which_insert_LEMON` | bit-wise: + 1 double, + 2 long, + 4 int flows and costs              | 3             |

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

Like CMake, the makefiles use the installed abLEMON directly and write the
same `runWarm()` into private copies of its headers. By default they look for
it under the `/usr` prefix; for other locations override
`LEMON_ROOT` (or, more finely, `LEMON_INCDIR`/`LEMON_LIBDIR`) on the `make`
command line or in the environment, for example:

```sh
make -f makefile LEMON_ROOT=/opt/local
```

Each executable using `MCFLemonSolver` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprising the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the `tests/MCFBlock`
tester of the umbrella project, which compares different solvers for
`MCFBlock` (be them specialised or general-purpose). The makefiles in turn recursively include all the
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
[open a new issue](https://gitlab.com/smspp/mcflemonsolver/-/issues/new).


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
see the [LICENSE](LICENSE) file for details. abLEMON is not redistributed: it
is consumed as an external package (under the Boost Software License, Version
1.0), and the build only rewrites private copies of five of its headers, as
described above.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

