module boundary
   use helpers
   implicit none

   

   contains

   subroutine apply_boundary_conditions (u, v, flag, imax, jmax, ui, vi)
     real u(0:imax+2, 0:jmax+2), v(0:imax+2, 0:jmax+2) 
     integer flag(0:imax+2, 0:jmax+2)
     real, intent(in) :: ui, vi
     integer, intent(in) :: imax, jmax

     integer:: i, j
     
     ! Fluid freely flows in from the west
     u(0,0:jmax+1) = u(1,0:jmax+1) 
     v(0,0:jmax+1) = v(1,0:jmax+1) 
     ! Fluid freely flows out to the east
     u(imax,0:jmax+1) = u(imax+2,0:jmax+2) 
     v(imax+1,0:jmax+1

     ! The vertical velocity approaches 0 at the north and south
     ! boundaries, but fluid flows freely in the horizontal direction 
     v(0:imax+1, jmax) = 0.0
     u(0:imax+1, jmax+1) = u(0:imax+1, jmax)
     v(0:imax+1, 0) = 0.0
     u(0:imax+1, 0) = u(0:imax+1, 1)

     
   
   end subroutine apply_boundary_conditions
   
   
end module boundary
