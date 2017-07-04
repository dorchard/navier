module output_mod
  use helpers_mod
  implicit none

  contains

  subroutine write_ppm (u, v, p, flag, iters, outname)
    real p(0:imax+1, 0:jmax+1), u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1)
    integer flag(0:imax+1, 0:jmax+1)
    integer, intent(in) :: iters
    character(10), intent(in) :: outname

    integer i, j

    real :: zmax = -1e10, zmin = 1e10
    real :: pmax = -1e10, pmin = 1e10
    real :: umax = -1e10, umin = 1e10
    real :: vmax = -1e10, vmin = 1e10

    character(64) :: outpath
    character(54) :: outfile

    real zeta(0:imax+1, 0:jmax+1)
    real zval
    integer r, g, b

    write(outfile, '(1i6.6)') (iters/output_freq)
    outpath = (trim(outname))//"/"//(trim(outfile))//".ppm"

    call calc_zeta(u, v, flag, zeta)

    open(unit=8,file=outpath,status='REPLACE')
    ! "P6 %d %d 255\n", imax, jmax
    write(8,'("P6 ", i3, " ", i3, " 255")') imax, jmax

    do j = 1, jmax
      do i = 1, imax
        if (.not.(toLogical(iand(flag(i, j), C_F)))) then
          r = 0
          b = 0
          g = 255
        else
          if (i < imax .and. j < jmax) then
            zval = zeta(i, j)
          else
            zval = 0.0
          end if

          r = (abs(zval/12.6)**0.4) * 255
          g = r
          b = r
          if (r > 255) then
            if (r < 510) then
              g = 255 - (r - 255)
              b = 255 - (r - 255)
              r = 255
            else
              b = 255 - (r - 510)
              g = 255
              r = 255
            end if
          end if
        end if
        write(8,"(A1,A1,A1,$)") char(r), char(g), char(b)
      end do
    end do

    close(8)
  end subroutine write_ppm

  subroutine calc_zeta(u, v, flag, zeta)
    real p(0:imax+1, 0:jmax+1), u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1)
    real zeta(0:imax+1, 0:jmax+1)
    integer flag(0:imax+1, 0:jmax+1)
    integer :: i, j

    ! Computation of the vorticity zeta at the upper right corner
    ! of cell (i,j) (only if the corner is surrounded by fluid cells)
    do i = 1, (imax-1)
      do j = 1, (jmax-1)
        if (toLogical(iand(flag(i,j), C_F)) .and.     &
            toLogical(iand(flag(i+1,j), C_F)) .and.   &
            toLogical(iand(flag(i,j+1), C_F)) .and.   &
            toLogical(iand(flag(i+1,j+1), C_F))) then
          zeta(i,j) = (u(i,j+1)-u(i,j))/dely - (v(i+1,j)-v(i,j))/delx
        else
          zeta(i,j) = 0.0
        end if
      end do
    end do
  end subroutine calc_zeta
end module output_mod
