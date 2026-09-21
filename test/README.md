# test

A tester for `MCFLemonSolver`.

It builds a small Min-Cost Flow instance in memory, whose optimal value is
known at every step, and attaches to its `MCFBlock` the four LEMON algorithms
on the `MCFListDigraph`. The `MCFBlock` is then changed in every way a
`MCFBlock` can change (costs, capacities, closing and opening arcs, adding and
removing arcs), and after each change the value of every attached `Solver` is
compared with the known one, and so is that of a `Solver` of each of the eight
kinds (four algorithms on two graphs) attached anew to the changed `MCFBlock`,
which has to build its graph with the arcs closed and deleted as they are.

The cases that the solution of an instance does not reach are checked on their
own: a `Solver` destroyed without having ever been attached to a `MCFBlock`, an
instance with no capacities and no costs that a `Modification` provides later,
deficits that sum to zero only up to rounding and deficits that do not, the
parameters looked up by their name, the invalid values of the parameters of
the algorithms, and the DMX file of the instance.

The exit code is the number of failed checks, zero printing
`All tests passed!!`. The `makefile` builds the executable including the
`MCFLemonSolver` and `MCFBlock` modules and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
