program main
   use mpi_f08
   use stdlib_linalg_constants, only: dp
   use stdlib_io_npy, only: save_npy
   use LightKrylov, only: cg, cg_dp_opts
   use LightKrylov_Logger, only: logger_setup, check_info
   use params
   use poisson
   implicit none

   character(len=128), parameter :: this_module = "Poisson solver"

   !> 2D Laplacian operator.
   type(Laplacian) :: L
   !> Preconditioner.
   type(blk_jacobi_precond) :: P
   !> Forcing term.
   type(vector) :: f
   !> Solution vector.
   type(vector) :: u

   !> Miscellaneous.
   integer :: info
   type(cg_dp_opts) :: opts

   !----- INITIALIZATION -----
   !> Initialize MPI.
   info = initialize()

   !> Logging with LightKrylov.
   call logger_setup()

   !> Initialize right-hand side vector.
   f = create_rhs()

   !> Set the initial guess to zero.
   call u%zero()

   !----- SOLVE POISSON EQUATION -----
   !> Sets options for the conjugate gradient solver.
   opts = cg_dp_opts(maxiter=nx*ny, if_print_metadata=.false.)
   !> Solve the linear system with block Jacobi preconditioning.
   call cg(L, f, u, info, preconditioner=P, rtol=epsilon(1.0_dp), options=opts)
   !> Sanity check.
   call check_info(info, "cg", module=this_module, procedure="main")

   ! ----- CHECK CONVERGENCE -----
   block
      type(vector) :: wrk
      real(dp) :: residual
      call L%matvec(u, wrk)
      call wrk%sub(f)
      residual = wrk%norm()
      if (nid == 0) then
         print *
         print *, "----------------------------------------"
         print *, "Euclidean norm of the residual vector : ", residual
         print *, "----------------------------------------"
      end if
   end block

   info = finalize()

   ! call save_npy("solution.npy", u%u)

contains
   type(vector) function create_rhs() result(f)
      !> Allocate constant forcing field.
      allocate (f%u(istart - 1:iend + 1, jstart - 1:jend + 1), source=1.0_dp)
      !> Boundary points are not actual degrees of freedom.
      if (istart == 1) f%u(istart - 1, :) = 0.0_dp
      if (iend == nx) f%u(iend + 1, :) = 0.0_dp
      if (jstart == 1) f%u(:, jstart - 1) = 0.0_dp
      if (jend == ny) f%u(:, jend + 1) = 0.0_dp
      !> Normalize.
      call f%scal(1.0_dp/f%norm())
      call exchange_halo(f%u)
   end function create_rhs
end program main
