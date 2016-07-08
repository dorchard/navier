module helpers_mod
  implicit none

  save

  ! flags for obstacle/boundary and fluid cells
  integer, parameter :: C_B = 0, C_F = 16
  ! obstacle with fluid cell to the _X
  integer, parameter :: B_N = 1, B_S = 2, B_W = 4, B_E = 8
  integer, parameter :: B_NW = ior(B_N, B_W), B_SW = ior(B_S, B_W),            &
                        B_NE = ior(B_N, B_E), B_SE = ior(B_S, B_E)
  integer, parameter :: B_NSEW = ior(B_N, ior(B_S, ior(B_E, B_W)))

  real, parameter ::  tau = 0.5

  real, parameter :: t_end = 40
  real, parameter :: eps = 0.001, omega = 1.7, gamma = 0.9
  real, parameter :: Re = 150.0, ui = 1.0, vi = 0.0

  integer, parameter :: imax = 660, jmax = 120, outputFlag = 1,                &
                        output_freq = 10, itermax = 10
  real, parameter :: output_freqt = 0.05
  real, parameter :: xlength = 22.0, ylength = 4.1
  real, parameter :: delx = xlength/imax, dely = ylength/jmax

  contains

  logical function toLogical(x)
    integer, intent(in) :: x
    if (x .eq. 0) then
       toLogical = .false.
    else
       toLogical = .true.
    end if
  end function

  integer function fromLogical(x)
    logical, intent(in) :: x
    if (x) then
       fromLogical = 0
    else
       fromLogical = 1
    end if
  end function

  integer function toMask(x)
    integer, intent(in) :: x
    if (x == 0) then
       toMask = 0
    else
       toMask = 1
    end if
  end function

end module helpers_mod
