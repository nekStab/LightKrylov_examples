module params
   use stdlib_linalg_constants, only: dp
   use stdlib_optval, only: optval
   implicit none
   private

   !---------------------------------
   !-----     GLOBAL DOMAIN     -----
   !---------------------------------

   !> Number of points in the horizontal direction.
   integer, parameter, public :: nx = 128
   !> Number of points in the vertical direction.
   integer, parameter, public :: ny = 128
   !> Domain length in the horizontal direction.
   real(dp), parameter, public :: Lx = 1.0_dp
   !> Domain length in the vertical direction.
   real(dp), parameter, public :: Ly = 1.0_dp
   !> Grid spacing in the horizontal direction.
   real(dp), parameter, public :: dx = Lx/(nx + 1)
   !> Grid spacing in the vertical direction.
   real(dp), parameter, public :: dy = Ly/(ny + 1)

   !--------------------------------------
   !-----     LOCAL DOMAIN (MPI)     -----
   !--------------------------------------

   !> Number of points in the horizontal subdomain.
   integer, protected, public :: nx_
   !> Number of points in the vertical subdomain.
   integer, protected, public :: ny_
   !> Starting and final index in the horizontal subdomain.
   integer, protected, public :: istart, iend
   !> Starting and final index in the vertical subdomain.
   integer, protected, public :: jstart, jend

   interface
      integer module function initialize() result(status)
      end function initialize
   end interface
   public :: initialize
contains
   module procedure initialize
   !> Sets the number of points in each direction of the subdomain.
   nx_ = nx; ny_ = ny
   !> Set the starting and final index in the horizontal direction.
   istart = 1; iend = nx_
   !> Set the starting and final index in the vertical direction.
   jstart = 1; jend = ny_
   !> Return flag.
   status = 0
   end procedure initialize
end module params
