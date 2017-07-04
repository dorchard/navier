module init_mod
  use helpers_mod
  implicit none

  contains

  subroutine init_flag(delx, dely, flag, ibound)
    real, intent(in) :: delx, dely
    integer flag(0:imax+1, 0:jmax+1)

    real :: mx, my, x, y, rad1, x1, y1
    integer :: ibound;

    integer :: i, j

    ! mask of a circular obstacle
    mx = 20.0/41.0*jmax*dely
    my = mx
    rad1 = 5.0/41.0*jmax*dely

    flag = C_F

    ! insert circular mask
    do i = 1, imax
      do j = 1, jmax
        x = (i-0.5)*delx - mx
        y = (j-0.5)*dely - my
        if (x*x + y*y <= rad1*rad1) then
          flag(i,j) = C_B
        end if

        ! Some code for adding a few more circular masks near-by
        !x1 = (i-0.5)*delx - mx*4
        !y1 = (j-0.5)*dely - my*1.25
        !if (x1*x1 + y1*y1 <= rad1*rad1) then
        !   flag(i,j) = C_B
        !end if
        !
        !x1 = (i-0.5)*delx - mx*3
        !y1 = (j-0.5)*dely - my*0.5
        !if (x1*x1 + y1*y1 <= rad1*rad1) then
        !   flag(i,j) = C_B
        !end if
      end do
    end do

    ! mark east & west bounary cells
    flag(0:imax+1,0) = C_B
    flag(0:imax+1,jmax+1) = C_B
    ! mark north & south boundary cells as obstacles
    flag(0,1:jmax) = C_B
    flag(imax+1,1:jmax) = C_B


    do i = 1, imax
      do j = 1, jmax
        if (.not. (toLogical (iand(flag(i,j),C_F)))) then
          ibound = ibound + 1
          if (toLogical (iand(flag(i-1,j), C_F))) flag(i,j) = ior(flag(i,j), B_W)
          if (toLogical (iand(flag(i+1,j), C_F))) flag(i,j) = ior(flag(i,j), B_E)
          if (toLogical (iand(flag(i,j-1), C_F))) flag(i,j) = ior(flag(i,j), B_S)
          if (toLogical (iand(flag(i,j+1), C_F))) flag(i,j) = ior(flag(i,j), B_N)
        end if
      end do
    end do
  end subroutine init_flag
end module init_mod
