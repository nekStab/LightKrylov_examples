# A collection of examples for LightKrylov

This repository contains a collection of examples illustrating various features of [`LightKrylov`](https://github.com/nekStab/LightKrylov).
A list of the different examples with a quick description is given below.

- [`spmv`](https://github.com/nekStab/LightKrylov_examples/tree/main/spmv) : Showcases the low overhead incurred by `LightKrylov`'s abstract object-oriented programming using a large sparse matrix-vector product as illustration.

- [`poisson`](https://github.com/nekStab/LightKrylov_examples/tree/main/poisson) : Showcases how to implement a simple MPI-based Poisson solver for the regular 5-points finite-difference discretization of the Laplace operator on the unit-square.

- [`neklab`](https://github.com/nekStab/LightKrylov_examples/tree/main/neklab) : Showcases how `LightKrylov` can be integrated into an already existing relatively large code base. The example computes the (unstable) steady solution of the nonlinear incompressible Navier-Stokes equations as well as the leading eigenpairs of the associated linearized operator.
