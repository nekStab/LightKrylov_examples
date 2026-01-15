module laplace
   use mpi_f08
   use stdlib_linalg_constants, only: dp
   use stdlib_optval, only: optval
   use specialmatrices, only: tridiagonal, solve
   use LightKrylov, only: abstract_vector_rdp, &
                          abstract_sym_linop_rdp, &
                          abstract_precond_rdp

   use params
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
   contains
      procedure, pass(self), public :: apply => apply_precond
   end type

contains

   module procedure initialize_vector
   allocate (vec%u(istart - 1:iend + 1, jstart - 1:jend + 1), source=u)
   call exchange_halo(vec%u)
   end procedure initialize_vector

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
      do concurrent(i=istart:iend, j=jstart:jend)
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
      do concurrent(i=istart:iend, j=jstart:jend)
         y(i, j) = alpha*x(i, j) + beta*y(i, j)
      end do
   end subroutine axpby_kernel

   subroutine rand(self, ifnorm)
      class(vector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      normalize = optval(ifnorm, .false.)
      if (.not. allocated(self%u)) allocate (self%u(istart - 1:iend + 1, jstart - 1:jend + 1))
      call random_number(self%u)
      call exchange_halo(self%u)
      if (normalize) call self%scal(1.0_dp/self%norm())
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
      real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1) :: u
      integer :: i, j
      select type (vec_in)
      type is (vector)
         select type (vec_out)
         type is (vector)
            !> Allocate return vector.
            call vec_out%zero()
            !> Local variable (circumventing the intent(in) statement).
            u = vec_in%u
            !> Exchange halos.
            call exchange_halo(u)
            !> Matrix-vector product.
            call spmv_kernel(u, vec_out%u)
         end select
      end select
   end subroutine matvec

   pure subroutine spmv_kernel(u, v)
      real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1), intent(in) :: u
      real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1), intent(out) :: v
      integer :: i, j
      !> Interior domain.
      do concurrent(i=istart:iend, j=jstart:jend)
         v(i, j) = (-u(i + 1, j) + 2*u(i, j) - u(i - 1, j))/dx**2 &
                   + (-u(i, j + 1) + 2*u(i, j) - u(i, j - 1))/dy**2
      end do
      !> Top-bottom boundary conditions.
      if (jstart == 1) v(:, jstart - 1) = 0.0_dp
      if (jend == ny) v(:, jend + 1) = 0.0_dp
      !> Left-right boundary conditions.
      if (istart == 1) v(istart - 1, :) = 0.0_dp
      if (iend == nx) v(iend + 1, :) = 0.0_dp
   end subroutine spmv_kernel

   !----------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR THE PRECONDITIONNER     -----
   !----------------------------------------------------------------

   subroutine apply_precond(self, vec, iter, current_residual, target_residual)
      class(blk_jacobi_precond), intent(inout) :: self
      class(abstract_vector_rdp), intent(inout) :: vec
      integer, optional, intent(in) :: iter
      real(dp), optional, intent(in) :: current_residual, target_residual
      !> Internal variables.
      type(tridiagonal) :: M

      !> Initialize local matrix.
      M = tridiagonal(-1.0/dx**2, 4.0/dx**2, -1.0/dx**2, nx_)
      !> Solve block system.
      select type (vec)
      type is (vector)
         vec%u(istart:iend, :) = solve(M, vec%u(istart:iend, :))
      end select
   end subroutine apply_precond

end module laplace
