# Sparse matrix-vector product

**Warning -** To run this example, you first need to `unzip data/matrix.zip`.

**High‑level overview**

- `app/main.f90`
    - Declares a program that loads a sparse matrix from a Harwell–Boeing file, creates random dense vectors, and wraps them into LightKrylov data structures (`dense_vector_rdp` and a custom `sparse_linop`).
    - Uses the stdlib routine `spmv` to perform a reference sparse matrix‑vector multiplication.
    - Benchmarks the stdlib implementation against the LightKrylov implementation (`linop%matvec`) to illustrate the low overhead introduced by LightKrylov’s object‑oriented abstractions.
    - Results are written to a `results/` directory via the `benchmark` type from `forbenchmark`.

- `src/utils.f90`
    - Provides the `utils` module that contains:
        - An interface `load` that reads a Harwell–Boeing file and returns a `csc_dp_type` matrix. The implementation (`load_harwell_boeing`) parses the file header, allocates the CSC arrays, and reads column pointers, row indices, and values.
        - A type `sparse_linop` that extends `abstract_linop_rdp`. It holds a `csc_dp_type` matrix and implements two procedures:
        - `matvec`: calls `spmv` on the stored matrix for a dense input vector, writing the result to a dense output vector.
        - `rmatvec`: performs the transpose multiplication via `spmv` with the `sparse_op_transpose` flag.
        - Both procedures use `select type` blocks to ensure that the passed vectors are of type `dense_vector_rdp`, zeroing the output vector before the multiplication.

**Intended message**

The two files together demonstrate that wrapping the standard `spmv` call inside LightKrylov’s linear‑operator interface adds negligible overhead, while providing a clean, object‑oriented API for matrix‑vector products.
