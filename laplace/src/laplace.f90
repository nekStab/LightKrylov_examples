module laplace
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, laplace!"
  end subroutine say_hello
end module laplace
