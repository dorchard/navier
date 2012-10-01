program foo
  implicit none

  integer, parameter :: c = 0
  call set(c)

  contains

  subroutine set(x)
     integer :: x
     x = 1
   end subroutine set

end program foo
