# LEMONSolver

This project implements `LEMONSolver`, a SMS++ :Solver
for [MCFBlock](https://gitlab.com/smspp/mcfblock) based
on interfacing solvers from the
[LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon).
In fact, `LEMONSolver` is not a single solver but no less
than 32 classes obtained instantiating 4 different "base"
solvers ...

## TODO: COMPLETE DESCRIPTION


## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and
  its requirements

- [MCFBlock](https://gitlab.com/smspp/mcfblock) and its
  requirements (but no actual `MCFSolver` are needed, since
  `LEMONSolver` provides its own)

- The [LEMON PROJECT](https://lemon.cs.elte.hu/trac/lemon)


## TODO: UPDATE

### Build and install with CMake

First you have to separately build LEMON with

```sh
cd lemon-development
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../lib
make
make install
cd ..
rm -Rf build
```

This will build a LEMON instance in
lemon-development/lib (getting rid of the build folder).

Then configure and build the library with:

```sh
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../lib
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(MCFBlock)
target_link_libraries(<my_target> SMS++::MCFBlock)
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

Each executable using `MCFBlock` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the
[tester comparing MCFBlock + MCFSolver with direct usage of the
original :MCFClass solver](test/test.cpp) alluded to in the previous section.
The makefiles in turn recursively include all the required other makefiles,
hence one should only need to edit the "main makefile" for compilation type
(C++ compiler and its options) and it all should be good to go. In case some
of the external libraries are not at their default location, it should only be
necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

Note thar the [MCFClass
project](https://github.com/frangio68/Min-Cost-Flow-Class) has a similar
arrangement with its own extlib/ folder, but due to some magic it is noy
necessary to that must be independently edit it in an analogous way.


## Tools

We provide a simple tool that converts MCF instances written in the DIMACS
standard into netCDF files. Optionally it hacks into the netCDF file to
change the number of static and dynamic nodes and arcs, as well as the
maximum number of nodes and arcs.

You can run the tool from the `<build-dir>/tools` directory or install it
with the library (see above). Run the tool without arguments for info on
its usage:

```sh
dmx2nc4
```


## Data

We provide a small sample of small-to-mid-size MCF problems in the
[data](data) folder. The instances comes compressed in the `dmx.tgz` file
in [data/dmx](data/dmx). Once this is decompressed and `data/nc4` is
created, the netCDF versions of the instances can be created in there by
running the `batch` file in the `tools` folder.


## Tests

The [test](test) folder contains a tester that reads an instance of a MCF
from a file (in either DIMACS or netCDF format) in an `MCFBlock`, and from
there in an object of a class MCFC derived from `MCFClass`, as decided by
the macro `WHICH_MCF`. Then, a `MCFSolver< MCFC >` is attached to the
`MCFBlock`. The MCF problem is then repeatedly solved with several changes in
costs / capacities / deficits, arcs openings / closures and arcs additions /
deletions. The same operations are performed on the two solvers, and the
results are compared. This mostly tests `MCFBlock` and `MCFSolver`, since
the actual `MCFClass` solved is the same, and so it can easily be wrong in
the same way for both the objects. The `batch` file tests basically only one
instance but in many different configurations (there can actually be two
`MCFBlock`, one of which is modified and the other solved, in all possible
combinations) and repeatedly, while the `batch-l` tests only the simplest
case but on several different problems of the [data](data) folder.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/mcfblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
