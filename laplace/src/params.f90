module params
   use mpi_f08
   use stdlib_linalg_constants, only: dp
   use stdlib_optval, only: optval
   use LightKrylov_Constants, only: set_rank, set_comm_size, set_io_rank
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

   !---------------------------------
   !-----     MPI VARIABLES     -----
   !---------------------------------

   !> Communicateur.
   type(mpi_comm), public :: world
   !> Number of mpi processes and rank of each process.
   integer, public :: nprocs, nid
   !> Topology information.
   integer, parameter :: ndim = 2
   integer, dimension(ndim) :: dims
   integer, dimension(ndim) :: coords
   logical, dimension(ndim), parameter :: periods = .false.
   logical, parameter :: reorganisation = .true.
   integer, parameter, public :: N = 1, E = 2, S = 3, W = 4
   !> Neighbourhood.
   integer, dimension(4), public :: neighbors
   !> Derived-types for halo exchanges.
   type(mpi_datatype), public :: dp_type
   type(mpi_datatype) :: row_type, col_type
   !> MPI constants.
   integer, public :: code

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

      integer module function finalize() result(status)
      end function finalize

      module subroutine exchange_halo(u)
         real(dp), dimension(istart - 1:iend + 1, jstart - 1:jend + 1), intent(inout) :: u
      end subroutine exchange_halo
   end interface
   public :: initialize
   public :: finalize
   public :: exchange_halo
contains
   module procedure initialize
   !---------------------------------------------
   !-----     GLOBAL MPI INITIALIZATION     -----
   !---------------------------------------------

   !> Initialize MPI.
   call mpi_init()
   !> Size of the MPI communicator.
   call mpi_comm_size(mpi_comm_world, nprocs)
   !> Number of processes along each direction.
   dims = 0; call mpi_dims_create(nprocs, ndim, dims)
   !> Create cartesian topology.
   call mpi_cart_create(mpi_comm_world, ndim, dims, periods, reorganisation, world)
   !> Coordinates of each process in the topology.
   call mpi_comm_rank(world, nid)
   call mpi_cart_coords(world, nid, ndim, coords)
   !> Find the neighbors of each process.
   call mpi_cart_shift(world, 0, 1, neighbors(W), neighbors(E))
   call mpi_cart_shift(world, 1, 1, neighbors(S), neighbors(N))

   !--------------------------------------------------
   !-----     LET LIGHTKRYLOV KNOW ABOUT MPI     -----
   !--------------------------------------------------

   call set_comm_size(nprocs)
   call set_rank(nid)
   call set_io_rank(0)

   !----------------------------------------
   !-----     LOCAL SUBDOMAIN DATA     -----
   !----------------------------------------

   !> Sets the number of points in each direction of the subdomain.
   nx_ = nx/dims(1); ny_ = ny/dims(2)
   !> Set the starting and final index in the horizontal direction.
   istart = coords(1)*nx/dims(1) + 1
   iend = (coords(1) + 1)*nx/dims(1)
   !> Set the starting and final index in the vertical direction.
   jstart = coords(2)*ny/dims(2) + 1
   jend = (coords(2) + 1)*ny/dims(2)

   !> Derived-type for communications.
   call mpi_type_create_f90_real(precision(1.0_dp), range(1.0_dp), dp_type)
   call mpi_type_vector(ny_, 1, nx_ + 2, dp_type, row_type)
   call mpi_type_commit(row_type)
   call mpi_type_contiguous(nx_, dp_type, col_type)
   call mpi_type_commit(col_type)

   !> Return flag.
   status = 0
   end procedure initialize

   module procedure finalize
   call mpi_comm_free(world)
   call mpi_type_free(row_type)
   call mpi_type_free(col_type)
   call mpi_finalize()
   status = 0
   end procedure finalize

   module procedure exchange_halo
   integer, parameter :: tag = 100
   type(mpi_status) :: status
   !> Send to north and receive from south.
   call mpi_sendrecv(u(istart, jstart), 1, row_type, neighbors(N), tag, &
                     u(iend + 1, jstart), 1, row_type, neighbors(S), tag, &
                     world, status)
   !> Send to south and receive from north.
   call mpi_sendrecv(u(iend, jstart), 1, row_type, neighbors(S), tag, &
                     u(istart - 1, jstart), 1, row_type, neighbors(N), tag, &
                     world, status)
   !> Send to west and receive from east.
   call mpi_sendrecv(u(istart, jstart), 1, col_type, neighbors(W), tag, &
                     u(istart, jend + 1), 1, col_type, neighbors(E), tag, &
                     world, status)
   !> Send to east and receive from west.
   call mpi_sendrecv(u(istart, jend), 1, col_type, neighbors(E), tag, &
                     u(istart, jstart - 1), 1, col_type, neighbors(W), tag, &
                     world, status)
   end procedure exchange_halo
end module params
