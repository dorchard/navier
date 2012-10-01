module helpers
  implicit none

  save

  integer, parameter :: C_B = 0, C_F = 16 ! flags for obstacle/boundary and fluid cells
  integer, parameter :: B_N = 1, B_S = 2, B_W = 4, B_E = 8 ! obstacle with fluid cell to the _X
  integer, parameter :: B_NW = B_N + B_W, B_SW = B_S + B_W, B_NE = B_N + B_E, B_SE = B_S + B_E
  integer, parameter :: B_NSEW = B_N + B_S + B_E + B_W 

  contains

  logical function toLogical(x)
    integer :: x
    if (x .eq. 0) then 
       toLogical = .false.
    else
       toLogical = .true.
    end if
  end function toLogical

end module helpers
