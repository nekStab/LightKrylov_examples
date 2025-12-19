module utils
   use stdlib_linalg_constants, only: dp
   use stdlib_sparse, only: csc_dp_type, spmv, sparse_op_transpose
   use LightKrylov_AbstractLinops, only: abstract_linop_rdp
   use LightKrylov_AbstractVectors, only: dense_vector_rdp, abstract_vector_rdp
   implicit none(type, external)
   private

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface load
      module function load_harwell_boeing(fname) result(A)
         character(len=*), intent(in) :: fname
        !! Name of the file to load.
         type(csc_dp_type) :: A
      end function load_harwell_boeing
   end interface
   public :: load

   !------------------------------------------
   !-----     LightKrylov data types     -----
   !------------------------------------------

   type, extends(abstract_linop_rdp), public :: sparse_linop
      type(csc_dp_type) :: A
   contains
      procedure, pass(self) :: matvec => sparse_matvec
      procedure, pass(self) :: rmatvec => sparse_rmatvec
   end type

contains

   !--------------------------------------
   !-----     Utility functions     ------
   !--------------------------------------

   module procedure load_harwell_boeing
   character(len=72) :: title
   character(len=8)  :: key
   character(len=3)  :: mxtype
   character(len=16) :: ptrfmt, indfmt
   character(len=20) :: valfmt, rhsfmt
   integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol, nnz, neltvl
   integer :: i, io

   !------------------------------------------------
   !-----     READ THE HARWELL-BEOING FILE     -----
   !------------------------------------------------
   open (newunit=io, file=fname, status="old", action="read")
   ! > Read in header block.
   read (io, 1000) title, key, &
      totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
      mxtype, nrow, ncol, nnz, neltvl, &
      ptrfmt, indfmt, valfmt, rhsfmt
1000 format(A72, A8/5I14/A3, 11X, 4I14/2A16, 2A20)
   !> Allocate data.
   call A%malloc(num_rows=nrow, num_cols=ncol, nnz=nnz)
   !> Read matrix structure.
   read (io, ptrfmt) (A%colptr(i), i=1, ncol + 1)
   read (io, indfmt) (A%row(i), i=1, nnz)
   !> Read values.
   if (valcrd > 0) read (io, valfmt) (A%data(i), i=1, nnz)
   close (io)
   end procedure load_harwell_boeing

   !------------------------------------------
   !-----     LightKrylov data types     -----
   !------------------------------------------

   subroutine sparse_matvec(self, vec_in, vec_out)
      class(sparse_linop), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_in)
      type is (dense_vector_rdp)
         select type (vec_out)
         type is (dense_vector_rdp)
            call vec_out%zero() ! Easy fix to allocate the data structure. (Not cool though)
            call spmv(self%A, vec_in%data, vec_out%data)
         end select
      end select
   end subroutine

   subroutine sparse_rmatvec(self, vec_in, vec_out)
      class(sparse_linop), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_in)
      type is (dense_vector_rdp)
         select type (vec_out)
         type is (dense_vector_rdp)
            call vec_out%zero() ! Easy fix to allocate the data structure. (Not cool though)
            call spmv(self%A, vec_in%data, vec_out%data, op=sparse_op_transpose)
         end select
      end select
   end subroutine

end module utils

