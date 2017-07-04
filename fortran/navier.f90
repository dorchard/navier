program navier
  use init_mod
  use helpers_mod
  use boundary_mod
  use simulation_mod
  use output_mod

  ! by D. Orchard (2012) based on the classic code from:
  !
  ! Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
  ! Numerical Simulation in Fluid Dynamics,
  ! SIAM, 1998.
  !
  ! http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

  implicit none

  character(10) :: outname = "output"

  real :: t = 0, res = 0, del_t = 0.003, t_count = 0.0
  integer :: i, j, itersor = 0, ifluid = 0, ibound = 0, init_case, iters = 0
  integer :: output_counter = 0

  real u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1), p(0:imax+1, 0:jmax+1)
  real rhs(0:imax+1, 0:jmax+1), f(0:imax+1, 0:jmax+1), g(0:imax+1, 0:jmax+1)
  integer flag(0:imax+1, 0:jmax+1)

  u = ui
  v = vi
  p = 0.0

  ! init flags
  call init_flag(delx, dely, flag, ibound)

  call apply_boundary_conditions(u, v, flag, 0.0)

  !do t = 0.0, t_end, del_t
  do while (t < t_end)
    del_t = set_timestep_interval(u, v, del_t)
    ifluid = (imax * jmax) - ibound

    call compute_tentative_velocity(u, v, f, g, flag, del_t)

    call compute_rhs(f, g, rhs, flag, del_t)

    if (ifluid > 0) then
      itersor = poisson(p, rhs, flag, res, ifluid)
    else
      itersor = 0
    end if

    print '(i0.1, " t:", f0.7, ", del_t:", f0.7, ", SOR iters:", i3.3, ", res:", f0.7, ", &
             bcells:", i0)', iters, t+del_t, del_t, itersor, res, ibound

    !("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
    !                iters, t+del_t, del_t, itersor, res, ibound);

    call update_velocity(u, v, f, g ,p, flag, del_t)

    call apply_boundary_conditions(u, v, flag, t)

    ! An alternate way to print output using the simulation time
    !if (toLogical(outputFlag) .and. (t >= t_count)) then
    !    call write_ppm(u, v, p, flag, output_counter, outname)
    !    output_counter = output_counter + 1
    !    t_count = t_count + output_freqt
    !end if

    if (toLogical(outputFlag) .and. (mod(iters, output_freq) == 0)) then
      call write_ppm(u, v, p, flag, iters, outname)
      output_counter = output_counter + 1
      t_count = t_count + output_freqt
    end if

    t = t + del_t
    iters = iters + 1
  end do
end program navier
