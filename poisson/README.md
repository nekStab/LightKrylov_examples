# Solving the Poisson equation on the unit-square with preconditioned conjugate gradient

This project demonstrates a minimal example of solving the two‑dimensional Poisson equation with a block‑Jacobi preconditioner, using the **LightKrylov** library as a linear‑algebra abstraction layer.

### Key Components

| Section | Purpose |
|---------|---------|
| **Derived types** | *`vector`* – extends `abstract_vector_rdp` and stores a distributed array `u`. <br>*`Laplacian`* – extends `abstract_sym_linop_rdp`; implements the matrix‑vector product via a stencil. <br>*`blk_jacobi_precond`* – extends `abstract_precond_rdp`; contains a tridiagonal matrix that is solved in each block. |
| **MPI/halo handling** | Global and local indices (`istart`, `iend`, `jstart`, `jend`) together with `exchange_halo` ensure that halo regions of the distributed arrays are kept consistent across MPI ranks. |
| **Vector operations** | Implements the required `abstract_vector_rdp` interface: `zero`, `dot`, `scal`, `axpby`, `rand`, and `get_size`. All heavy work is off‑loaded to small, pure kernels (`dot_kernel`, `scal_kernel`, `axpby_kernel`) that operate over the local sub‑domain. |
| **Laplacian kernel** | `spmv_kernel` performs a 5‑point stencil for the discrete Laplacian, applying boundary conditions afterwards. |
| **Preconditioner** | `apply_precond` solves a local tridiagonal system with the `solve` routine from `specialmatrices`. |
| **Initialization helpers** | `initialize_vector` and `initialize_preconditioner` provide convenient factory functions that set up the data structures and exchange halos immediately. |
| **Utility routines** | `apply_boundary_conditions` enforces Dirichlet conditions on the physical boundaries. |

### How it fits into LightKrylov

1. **Abstract interfaces** – LightKrylov defines abstract base types (`abstract_vector_rdp`, `abstract_sym_linop_rdp`, `abstract_precond_rdp`).  
2. **Derived types** – The module creates concrete implementations that inherit from these bases, providing the concrete Fortran code for all required operations.  
3. **Generic kernels** – The lightweight kernels keep the heavy numeric work in pure functions, which can be vectorised or accelerated independently.  
4. **MPI‑aware design** – By including halo exchange in all mutating routines, the code remains correct in a distributed environment without requiring LightKrylov to manage communication explicitly.

When combined, this module can be passed to LightKrylov solvers (e.g., CG, MINRES, or GMRES) as `Laplacian` and `blk_jacobi_precond` objects, and the solver will automatically use the vector operations defined here. The result is a self‑contained, MPI‑ready Poisson solver that illustrates how to plug a Fortran implementation into the LightKrylov framework.
