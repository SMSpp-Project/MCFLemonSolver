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
type). Since these can be many, the macro `SMSpp\_which\_insert\_LEMON'
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

The interface with the solvers provided by the LEMON project works thanks to the
changes made to the `include/lemon/bits/array_map.h` file. Since it has not been
updated for some time, it was not compatible with C++20 versions, causing
problems, especially with the compilation of Concepts. This file is provided
within the include/lemon directory, modified to ensure correctness, and a patch
has been proposed to the LEMON developers, who I assume will modify this file to
make it available in their repositories.

The files are provided "working", so the user does not need to download LEMON
themselves. If they did, the `array_map.h` file downloaded from LEMON might not
have been updated yet, potentially compromising the compilation of
`MCFLemonSolver`.

If LEMON provides a new release and the user wants to update, it is necessary
to check the correctness of the `include/lemon/bits/array_map.h` file,
particularly the correct use of `std::allocator`.

It is also necessary to fix the following line of code in the
`include/lemon/adaptors.h` file:

```c++
#ifdef _MSC_VER
#define LEMON_SCOPE_FIX(OUTER, NESTED) OUTER::NESTED
#else
#define LEMON_SCOPE_FIX(OUTER, NESTED) typename OUTER::template NESTED
#endif
```

and change it to:

```c++
#define LEMON_SCOPE_FIX(OUTER, NESTED) typename OUTER::template NESTED
```

otherwise it will not compile on Windows with MSVC.

## Getting started

These instructions will let you build `MCFLemonSolver` on your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and
  its requirements.

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its
  requirements (but no actual `MCFSolver` are needed, since
  `LEMONSolver` provides its own).

- The [LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon), currently contained
  in this repository: do not download the LEMON Project from their site if there
  aren't new releases, for the reasons discussed above, but use what is provided
  as "functional" in this repository.


## TODO: UPDATE

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

Optionally, install the library in the system with:

```sh
cmake --install .
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(LEMONSolver)
target_link_libraries(<my_target> SMS++::LEMONSolver)
```


### Running the tests with CMake

A unit test will be built with the library.
To disable it, set the option `BUILD_TESTING` to `OFF`.

The test takes an instance of a MCF in DIMACS or NC4 format. The MCF problem
is then repeatedly solved with several changes in costs/capacities/deficits,
arcs openings/closures and arcs additions/deletions. The same operations are
performed on the two solvers, and the results are compared.


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

To build code with LEMONSolver dependencies you will first need to separately
build LEMON with:

```sh
cd lemon-development
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../lib
cmake --build .
cmake --install .
cd ..
rm -Rf build
```

This will build a LEMON instance in
`lemon-development/lib` (getting rid of the build folder).

Then configure and build the library with:

```sh
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../lib
cmake --build .
```

Each executable using `LEMONSolver` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the
[tester comparing MCFBlock + MCFLemonSolver with direct usage of the
original solver](test/test.cpp) alluded to in the previous section.
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries are not at their default location, it should only be
necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/mcfblock/-/issues/new).


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


## License

The SMS++ code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details. Due to the above-mentioned
issues we are also redistributing a (lightly patched) version of the latest
LEMON version available at the time of released, which is separately
available under the (more permissive) Boost Software License, Version 1.0;
see the [LICENSE](lemon-development/LICENSE) file in the
[lemon-development](lemon-development/) folder.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

