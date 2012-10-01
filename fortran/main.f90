program navier
   use init
   use helpers
   use boundary

   implicit none
         
   real, parameter :: xlength = 22.0, ylength = 4.1
   real, parameter :: t_end = 40, del_t = 0.003
   real, parameter :: eps = 0.001, omega = 1.7, gamma = 0.9
   real, parameter :: Re = 150.0, ui = 1.0, vi = 0.0


   integer, parameter :: imax = 660, jmax = 120, output = 1, output_freq = 100, itermax = 100
   character(10) :: outname

   real :: t, delx = xlength/imax, dely = ylength/jmax, res
   integer :: i, j, itersor = 0, ifluid = 0, ibound = 0, init_case, iters = 0
   
   real u(0:imax+1, 0:jmax+1), v(0:imax+1, 0:jmax+1), p(0:imax+1, 0:jmax+1)
   real rhs(0:imax+1, 0:jmax+1), f(0:imax+1, 0:jmax+1), g(0:imax+1, 0:jmax+1)
   integer flag(0:imax+1, 0:jmax+1)

   !do i = 0, imax+1, 1
   !   do j = 0, jmax+1, 1
   !      u
   !   end do
   !end do
   u = 0.0
   v = 0.0
   p = 0.0
   
   ! init flags
   call initFlag(imax,jmax,delx,dely,flag,ibound)
  
   call apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi)
   
   contains
     

end program navier
