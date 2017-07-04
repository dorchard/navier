module boundary_mod
  use helpers_mod
  implicit none

  contains

  ! Given the boundary conditions defined by the flag matrix, update
  ! the u and v velocities. Also enforce the boundary conditions at the
  ! edges of the matrix.
  subroutine apply_boundary_conditions (u, v, flag, t)
    real u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1)
    integer flag(0:imax+1, 0:jmax+1)
    real t

    integer:: i, j

    ! Fluid freely flows in from the west
    u(0,0:jmax+1) = u(1,0:jmax+1)
    v(0,0:jmax+1) = v(1,0:jmax+1)
    ! Fluid freely flows out to the east
    u(imax,0:jmax+1) = u(imax-1,0:jmax+1)
    v(imax+1,0:jmax+1) = v(imax, 0:jmax+1)

    ! The vertical velocity approaches 0 at the north and south
    ! boundaries, but fluid flows freely in the horizontal direction
    v(0:imax+1, jmax) = 0.0
    u(0:imax+1, jmax+1) = u(0:imax+1, jmax)
    v(0:imax+1, 0) = 0.0
    u(0:imax+1, 0) = u(0:imax+1, 1)

    ! Apply no-slip boundary conditions to cells that are adjacent to
    ! internal obstacle cells. This doces the u and v velocity to
    ! tend towards zero in these cells.

    do i = 1, imax
      do j = 1, jmax
        if (toLogical(iand(flag(i,j), B_NSEW))) then
          select case (flag(i, j))
            case (B_N)
              u(i, j) = -u(i, j+1)
            case (B_E)
              u(i, j) = 0.0
            case (B_NE)
              u(i, j) = 0.0
            case (B_SE)
              u(i, j) = 0.0
            case (B_NW)
              u(i, j) = -u(i, j+1)
            case (B_S)
              u(i, j) = -u(i, j-1)
            case (B_SW)
              u(i, j) = -u(i, j-1)
          end select
        end if
      end do
    end do

    do i = 0, (imax-1), 1
      do j = 1, jmax, 1
        if (toLogical(iand(flag(i+1,j), B_NSEW))) then
          select case (flag(i+1,j))
            case (B_N)
              u(i,j) = -u(i,j+1)
            case (B_W)
              u(i,j) = 0.0
            case (B_NE)
              u(i,j) = -u(i,j+1)
            case (B_SW)
              u(i,j) = 0.0
            case (B_NW)
              u(i,j) = 0.0
            case (B_S)
              u(i,j) = -u(i,j-1)
            case (B_SE)
              u(i,j) = -u(i,j-1)
          end select
        end if
      end do
    end do

    do i = 1, imax, 1
      do j = 1, jmax, 1
        if (toLogical(iand(flag(i,j), B_NSEW))) then
          select case (flag(i,j))
            case (B_N)
              v(i,j) = 0.0
            case (B_E)
              v(i,j) = -v(i+1,j)
            case (B_NE)
              v(i,j) = 0.0
            case (B_SE)
              v(i,j) = -v(i+1,j)
            case (B_NW)
              v(i,j) = 0.0
            case (B_W)
              v(i,j) = -v(i-1,j)
            case (B_SW)
              v(i,j) = -v(i-1,j)
          end select
        end if
      end do
    end do

    do i = 1, imax, 1
      do j = 0, jmax-1, 1
        if (toLogical(iand(flag(i,j+1), B_NSEW))) then
          select case (flag(i,j+1))
            case (B_E)
              v(i,j) = -v(i+1,j)
            case (B_S)
              v(i,j) = 0.0
            case (B_NE)
              v(i,j) = -v(i+1,j)
            case (B_SE)
              v(i,j) = 0.0
            case (B_SW)
              v(i,j) = 0.0
            case (B_W)
              v(i,j) = -v(i-1,j)
            case (B_NW)
              v(i,j) = -v(i-1,j)
          end select
        end if
      end do
    end do

    ! Finally, fix the horizontal velocity at the  western edge to have
    ! a continual flow of fluid into the simulation.

    v(0,0) = 2*vi-v(1,0)
    u(0, 1:jmax) = ui
    v(0, 1:jmax) = 2*vi-v(1,1:jmax)

    ! code for a conflicting flow coming in from the opposite side
    !if (t < 45) then
    !  u(imax, 1:jmax) = -ui
    !  u(imax-1, 1:jmax) = -ui
    !  u(imax+1, 1:jmax) = -ui
    !else
    !  u(imax,0:jmax+1) = u(imax-1,0:jmax+1)
    !end if


    ! silly tricks to drop pieces of high flow in at certain points during the sim - for fun
    !if (mod(t, 4.0) < 0.4 .and. t > 0.5) then
    !   do i = -1, 1, 1
    !      do j = -1, 1, 1
    !         u(mod(floor(t*10.0+i),imax),mod(floor(t*5+j),jmax)) = 1
    !         v(mod(floor(t*10.0+i),imax),mod(floor(t*5+j),jmax)) = 1
    !      end do
    !   end do
    !end if
  end subroutine apply_boundary_conditions
end module boundary_mod
