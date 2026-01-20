module poisson
   use mpi_f08
   !> Fortran Standard Library.
   use stdlib_linalg_constants, only: dp
   use stdlib_optval, only: optval
   !> SpecialMatrices package.
   use specialmatrices, only: tridiagonal, solve            ! For block Jacobi preconditioner.
   !> LightKrylov
   use LightKrylov, only: abstract_vector_rdp, &
                          abstract_sym_linop_rdp, &
                          abstract_precond_rdp
   !> Current solver.
   use params, only: dx, dy, inv_dx2, inv_dy2, &            ! Grid spacing.
                     nx, ny, nx_, ny_, &                    ! Global and local grid size.
                     istart, iend, jstart, jend, &          ! Local indices
                     code, world, dp_type, exchange_halo    ! MPI-related
   implicit none
   private

   !> Derived-type for the vector.
   type, extends(abstract_vector_rdp), public :: vector
      real(dp), allocatable :: u(:, :)
   contains
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type vector

   interface vector
      type(vector) module function initialize_vector(u) result(vec)
         real(dp), intent(in) :: u(:, :)
      end function initialize_vector
   end interface vector

   !> Derived-type for the Laplacian.
   type, extends(abstract_sym_linop_rdp), public :: Laplacian
   contains
      procedure, pass(self), public :: matvec
   end type Laplacian

   !> Derived-type for the preconditionner.
   type, extends(abstract_precond_rdp), public :: blk_jacobi_precond
      type(tridiagonal) :: M
   contains
      procedure, pass(self), public :: apply => apply_precond
   end type
   interface preconditioner
      type(blk_jacobi_precond) pure module function initialize_preconditioner() result(P)
      end function initialize_preconditioner
   end interface
   public :: preconditioner

contains

   module procedure initialize_vector
   allocate (vec%u(istart - 1:iend + 1, jstart - 1:jend + 1), source=u)
   call exchange_halo(vec%u)
   end procedure initialize_vector

   module procedure initialize_preconditioner
   P%M = tridiagonal(-1/dx**2, 4/dx**2, -1/dx**2, nx_)
   end procedure

   pure subroutine apply_boundary_conditions(u)
      real(dp), intent(inout) :: u(istart - 1:iend + 1, jstart - 1:jend + 1)
      !> Top-bottom boundary conditions.
      if (jstart == 1) u(:, jstart - 1) = 0.0_dp
      if (jend == ny) u(:, jend + 1) = 0.0_dp
      !> Left-right boundary conditions.
      if (istart == 1) u(istart - 1, :) = 0.0_dp
      if (iend == nx) u(iend + 1, :) = 0.0_dp
   end subroutine apply_boundary_conditions

   !---------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR DERIVED TYPE     -----
   !---------------------------------------------------------

   subroutine zero(self)
      class(vector), intent(inout) :: self
      if (allocated(self%u)) then
         self%u = 0.0_dp
      else
         allocate (self%u(istart - 1:iend + 1, jstart - 1:jend + 1), source=0.0_dp)
      end if
   end subroutine zero

   real(dp) function dot(self, vec) result(alpha)
      class(vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      real(dp) :: alpha_local
      integer :: i, j
      select type (vec)
      type is (vector)
         alpha_local = dot_kernel(self%u, vec%u)
         call mpi_allreduce(alpha_local, alpha, 1, dp_type, mpi_sum, world, code)
      end select
   end function dot

   real(dp) pure function dot_kernel(u, v) result(alpha)
      real(dp), intent(in) :: u(istart - 1:iend + 1, jstart - 1:jend + 1)
      real(dp), intent(in) :: v(istart - 1:iend + 1, jstart - 1:jend + 1)
      integer :: i, j
      alpha = 0.0_dp
      do concurrent(i=istart:iend, j=jstart:jend)
         alpha = alpha + u(i, j)*v(i, j)
      end do
      alpha = alpha*dx*dy
   end function dot_kernel

   subroutine scal(self, alpha)
      class(vector), intent(inout) :: self
      real(dp), intent(in) :: alpha
      call scal_kernel(alpha, self%u)
   end subroutine scal

   pure subroutine scal_kernel(alpha, u)
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: u(istart - 1:iend + 1, jstart - 1:jend + 1)
      integer :: i, j
      do concurrent(i=istart - 1:iend + 1, j=jstart - 1:jend + 1)
         u(i, j) = alpha*u(i, j)
      end do
   end subroutine scal_kernel

   subroutine axpby(alpha, vec, beta, self)
      real(dp), intent(in) :: alpha, beta
      class(abstract_vector_rdp), intent(in) :: vec
      class(vector), intent(inout) :: self
      select type (vec)
      type is (vector)
         call axpby_kernel(alpha, vec%u, beta, self%u)
      end select
   end subroutine axpby

   pure subroutine axpby_kernel(alpha, x, beta, y)
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: x(istart - 1:iend + 1, jstart - 1:jend + 1)
      real(dp), intent(inout) :: y(istart - 1:iend + 1, jstart - 1:jend + 1)
      integer :: i, j
      do concurrent(i=istart - 1:iend + 1, j=jstart - 1:jend + 1)
         y(i, j) = alpha*x(i, j) + beta*y(i, j)
      end do
   end subroutine axpby_kernel

   subroutine rand(self, ifnorm)
      class(vector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      normalize = optval(ifnorm, .false.)
      if (.not. allocated(self%u)) allocate (self%u(istart - 1:iend + 1, jstart - 1:jend + 1))
      call random_number(self%u(istart:iend, jstart:jend))
      if (normalize) call self%scal(1.0_dp/self%norm())
      call exchange_halo(self%u)
   end subroutine rand

   integer function get_size(self) result(n)
      class(vector), intent(in) :: self
      n = nx*ny
   end function get_size

   !----------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR THE LAPLACIAN     -----
   !----------------------------------------------------------

   subroutine matvec(self, vec_in, vec_out)
      class(Laplacian), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_in)
      type is (vector)
         select type (vec_out)
         type is (vector)
            !> Allocate return vector.
            call vec_out%zero()
            !> Matrix-vector product.
            call spmv_kernel(vec_in%u, vec_out%u)
            !> Exchange halos.
            call exchange_halo(vec_out%u)
         end select
      end select
   end subroutine matvec

   pure subroutine spmv_kernel(u, v)
      real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1), intent(in) :: u
      real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1), intent(out) :: v
      real(dp) :: tmp
      integer :: i, j
      !> Interior domain.
      do concurrent(i=istart:iend, j=jstart:jend)
         tmp = 2*u(i, j)
         v(i, j) = (-u(i + 1, j) + tmp - u(i - 1, j))*inv_dx2 &
                   + (-u(i, j + 1) + tmp - u(i, j - 1))*inv_dy2
      end do
      call apply_boundary_conditions(v)
   end subroutine spmv_kernel

   !----------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR THE PRECONDITIONNER     -----
   !----------------------------------------------------------------

   subroutine apply_precond(self, vec, iter, current_residual, target_residual)
      class(blk_jacobi_precond), intent(inout) :: self
      class(abstract_vector_rdp), intent(inout) :: vec
      integer, optional, intent(in) :: iter
      real(dp), optional, intent(in) :: current_residual, target_residual
      !> Solve block system.
      select type (vec)
      type is (vector)
         vec%u(istart:iend, :) = solve(self%M, vec%u(istart:iend, :))
         call exchange_halo(vec%u)
      end select
   end subroutine apply_precond

end module poisson
