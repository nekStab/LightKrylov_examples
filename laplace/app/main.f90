program main
   use stdlib_linalg_constants, only: dp
   use stdlib_io_npy, only: save_npy
   use LightKrylov, only: cg, cg_dp_opts
   use LightKrylov_Logger, only: logger_setup, check_info
   use laplace
   implicit none

   character(len=128), parameter :: this_module = "Poisson solver"

   !> 2D Laplacian operator.
   type(Laplacian) :: L
   !> Forcing term.
   type(vector) :: f
   !> Solution vector.
   type(vector) :: u

   !> Miscellaneous.
   integer :: info
   type(cg_dp_opts) :: opts

   !----- INITIALIZATION -----
   !> Logging with LightKrylov.
   call logger_setup()

   !> Initialize right-hand side vector.
   f = create_rhs()

   !> Set the initial guess to zero.
   call u%zero()

   !----- SOLVE POISSON EQUATION -----
   !> Sets options for the conjugate gradient solver.
   opts = cg_dp_opts(maxiter=nx*ny, if_print_metadata=.false.)
   !> Solve the linear system.
   call cg(L, f, u, info, options=opts)
   !> Sanity check.
   call check_info(info, "cg", module=this_module, procedure="main")

   ! ----- CHECK CONVERGENCE -----
   block
      type(vector) :: wrk
      call L%matvec(u, wrk)
      call wrk%sub(f)
      print *
      print *, "----------------------------------------"
      print *, "Euclidean norm of the residual vector : ", wrk%norm()
      print *, "----------------------------------------"
   end block

   call save_npy("solution.npy", u%u)

contains
   type(vector) function create_rhs() result(f)
      integer :: i, j
      !> Forcing at the interior of the domain.
      call random_number(f%u)
      !> Normalize.
      call f%scal(1.0_dp/f%norm())
      !> Boundary points are not actual degrees of freedom.
      f%u(0, :) = 0.0_dp; f%u(nx + 1, :) = 0.0_dp
      f%u(:, 0) = 0.0_dp; f%u(:, ny + 1) = 0.0_dp
   end function create_rhs
end program main
