module laplace
   use stdlib_linalg_constants, only: dp
   use stdlib_optval, only: optval
   use LightKrylov, only: abstract_vector_rdp, &
                          abstract_sym_linop_rdp
   implicit none
   private

   !> Number of points per direction.
   integer, parameter, public :: nx = 128
   integer, parameter, public :: ny = 128
   real(dp), parameter, public :: Lx = 1.0_dp
   real(dp), parameter, public :: Ly = 1.0_dp
   real(dp), parameter, public :: dx = Lx/(nx + 1)
   real(dp), parameter, public :: dy = Ly/(ny + 1)

   !> Derived-type for the vector.
   type, extends(abstract_vector_rdp), public :: vector
      real(dp), dimension(0:nx + 1, 0:ny + 1) :: u
   contains
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
      procedure, pass(self), public :: get_size
   end type vector

   interface vector
      pure module function initialize_vector(u) result(vec)
         real(dp), intent(in) :: u(:, :)
         type(vector) :: vec
      end function initialize_vector
   end interface vector

   !> Derived-type for the Laplacian.
   type, extends(abstract_sym_linop_rdp), public :: Laplacian
   contains
      procedure, pass(self), public :: matvec
   end type Laplacian

contains

   module procedure initialize_vector
   vec%u = u
   end procedure initialize_vector

   !---------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR DERIVED TYPE     -----
   !---------------------------------------------------------

   subroutine zero(self)
      class(vector), intent(inout) :: self
      self%u = 0.0_dp
   end subroutine zero

   real(dp) function dot(self, vec) result(alpha)
      class(vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      integer :: i, j
      select type (vec)
      type is (vector)
         alpha = dot_kernel(nx, ny, self%u, vec%u)
      end select
   end function dot

   real(dp) pure function dot_kernel(m, n, u, v) result(alpha)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: u(0:m + 1, 0:n + 1), v(0:m + 1, 0:n + 1)
      integer :: i, j
      alpha = 0.0_dp
      do concurrent(i=1:m, j=1:n)
         alpha = alpha + u(i, j)*v(i, j)
      end do
      alpha = alpha*dx*dy
   end function dot_kernel

   subroutine scal(self, alpha)
      class(vector), intent(inout) :: self
      real(dp), intent(in) :: alpha
      call scal_kernel(nx, ny, alpha, self%u)
   end subroutine scal

   pure subroutine scal_kernel(m, n, alpha, u)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: u(0:m + 1, 0:n + 1)
      integer :: i, j
      do concurrent(i=1:m, j=1:n)
         u(i, j) = alpha*u(i, j)
      end do
   end subroutine scal_kernel

   subroutine axpby(alpha, vec, beta, self)
      real(dp), intent(in) :: alpha, beta
      class(abstract_vector_rdp), intent(in) :: vec
      class(vector), intent(inout) :: self
      integer :: i, j
      select type (vec)
      type is (vector)
         call axpby_kernel(nx, ny, alpha, vec%u, beta, self%u)
      end select
   end subroutine axpby

   pure subroutine axpby_kernel(m, n, alpha, x, beta, y)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: x(0:m + 1, 0:n + 1)
      real(dp), intent(inout) :: y(0:m + 1, 0:n + 1)
      integer :: i, j
      do concurrent(i=0:m + 1, j=0:n + 1)
         y(i, j) = alpha*x(i, j) + beta*y(i, j)
      end do
   end subroutine axpby_kernel

   subroutine rand(self, ifnorm)
      class(vector), intent(inout) :: self
      logical, optional, intent(in) :: ifnorm
      logical :: normalize
      normalize = optval(ifnorm, .false.)
      call random_number(self%u)
      if (normalize) self%u = self%norm()
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
      integer :: i, j
      select type (vec_in)
      type is (vector)
         select type (vec_out)
         type is (vector)
            call spmv_kernel(nx, ny, vec_in%u, vec_out%u)
         end select
      end select
   end subroutine matvec

   pure subroutine spmv_kernel(m, n, u, v)
      integer, intent(in) :: m, n
      real(dp), dimension(0:m + 1, 0:n + 1), intent(in) :: u
      real(dp), dimension(0:m + 1, 0:n + 1), intent(out) :: v
      integer :: i, j
      !> Interior domain.
      do concurrent(i=1:m, j=1:n)
         v(i, j) = (-u(i + 1, j) + 2*u(i, j) - u(i - 1, j))/dx**2 &
                   + (-u(i, j + 1) + 2*u(i, j) - u(i, j - 1))/dy**2
      end do
      !> Top-bottom boundary conditions.
      v(:, 0) = 0.0_dp; v(:, n + 1) = 0.0_dp
      !> Left-right boundary conditions.
      v(0, :) = 0.0_dp; v(m + 1, :) = 0.0_dp
   end subroutine spmv_kernel

end module laplace
