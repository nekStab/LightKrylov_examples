program main
   use stdlib_linalg_constants, only: dp
   use stdlib_sparse, only: csc_dp_type, spmv, sparse_op_transpose
   use LightKrylov, only: dense_vector_rdp, dense_vector
   use forbenchmark
   use utils, only: load, sparse_linop
   implicit none(type, external)

   !> Sparse matrix.
   character(len=32) :: fname
   type(csc_dp_type) :: A
   !> Test vectors.
   real(dp), allocatable :: x(:), y(:)
   !> LightKrylov specifics.
   type(dense_vector_rdp) :: vec_in, vec_out
   type(sparse_linop) :: linop
   !> Benchmarking.
   type(benchmark) :: bench
   integer :: i

   !--------------------------------------------------------
   !-----     Load matrix in Harwell-Boeing format     -----
   !--------------------------------------------------------
   fname = "data/s3dkq4m2.dat"
   A = load(fname)

   print *, "-------------------------------------------"
   print *, "-----     Matrix characteristics      -----"
   print *, "-------------------------------------------"
   print *, ""
   print *, "Dimension :", A%nrows, "x", A%ncols
   print *, "nnz(A)    :", A%nnz
   print *, ""

   !> Create random vectors.
   allocate (x(A%ncols), y(A%ncols), source=0.0_dp)
   call random_number(x)

   !> Wrap into LightKrylov-specific types.
   linop = sparse_linop(A)
   vec_in = dense_vector(x)
   vec_out = dense_vector(y)

   !--------------------------------
   !-----     Benchmarking     -----
   !--------------------------------
   call bench%init(2, "Sparse matrix-vector product", "results/")
   bench%nloops = 10000

   !> Reference implementation for the sparse matrix-vector product from stdlib.
   call bench%start_benchmark(1, "Fortran stdlib")
   do i = 1, bench%nloops
      call spmv(A, x, y)
   end do
   call bench%stop_benchmark()

   !> Equivalent LightKrylov implementation.
   !> Goal: evaluate overhead from the abstraction level.
   call bench%start_benchmark(2, "LightKrylov")
   do i = 1, bench%nloops
      call linop%matvec(vec_in, vec_out)
   end do
   call bench%stop_benchmark()
end program main
