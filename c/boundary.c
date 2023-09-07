#include <stdio.h>
#include <string.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(double **u, double **v, char **flag,
    int imax, int jmax, double ui, double vi)
{
    int i, j;

    for (j=0; j<=jmax+1; j++) {
        /* Fluid freely flows in from the west */
        u[0][j] = u[1][j];
        v[0][j] = v[1][j];

        /* Fluid freely flows out to the east */
        u[imax][j] = u[imax-1][j];
        v[imax+1][j] = v[imax][j];
    }

    for (i=0; i<=imax+1; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v[i][jmax] = 0.0;
        u[i][jmax+1] = u[i][jmax];

        v[i][0] = 0.0;
        u[i][0] = u[i][1];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */

   for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & B_NSEW) {
                switch (flag[i][j]) {
                    case B_N: 
                        u[i][j]   = -u[i][j+1];
                        break;
                    case B_E: 
                        u[i][j]   = 0.0;
                        break;
                    case B_NE:
                        u[i][j]   = 0.0;
                        break;
                    case B_SE:
                        u[i][j]   = 0.0;                      
                        break;
                    case B_NW:
                        u[i][j]   = -u[i][j+1];
                        break;
                    case B_S:
                        u[i][j]   = -u[i][j-1];
                        break;
                    case B_SW:
                        u[i][j]   = -u[i][j-1];
                        break;
                }
            }
	}
    } 
    for (i=0; i<=(imax-1); i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i+1][j] & B_NSEW) {
                switch (flag[i+1][j]) {
                    case B_N: 
                        u[i][j] = -u[i][j+1];
                        break;
                    case B_W: 
                        u[i][j] = 0.0;
                        break;
                    case B_NE:
                        u[i][j] = -u[i][j+1];
                        break;
                    case B_SW:
                        u[i][j] = 0.0;
                        break;
                    case B_NW:
                        u[i][j] = 0.0;
                        break;
                  case B_S:
                        u[i][j] = -u[i][j-1];
                        break;
                    case B_SE:
                        u[i][j] = -u[i][j-1];
                        break;
                }
            }
	}
    } 


    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & B_NSEW) {
                switch (flag[i][j]) {
                    case B_N: 
                        v[i][j]   = 0.0;
                        break;
                    case B_E: 
                        v[i][j]   = -v[i+1][j];
                        break;
                    case B_NE:
                        v[i][j]   = 0.0;
                        break;
                    case B_SE:
                        v[i][j]   = -v[i+1][j];
                        break;
                    case B_NW:
                        v[i][j]   = 0.0;
                        break;
                    case B_W: 
                        v[i][j]   = -v[i-1][j];
                        break;
                    case B_SW:
                        v[i][j]   = -v[i-1][j];
                        break;
                }
            }
	 }
      } 

    for (i=1; i<=imax; i++) {
      for (j=0; j<=(jmax-1); j++) {
            if (flag[i][j+1] & B_NSEW) {
                switch (flag[i][j+1]) {
                    case B_E: 
                        v[i][j] = -v[i+1][j];
                        break;
                    case B_S:
                        v[i][j] = 0.0;
                        break;
                    case B_NE:
                        v[i][j] = -v[i+1][j];
                        break;
                    case B_SE:
                        v[i][j] = 0.0;
                        break;
                    case B_SW:
                        v[i][j] = 0.0;
			break;
                    case B_W: 
                        v[i][j] = -v[i-1][j];
                        break;
                    case B_NW:
                        v[i][j] = -v[i-1][j];
                        break;
                }
            }
	}
     } 

    /* for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & B_NSEW) {
                switch (flag[i][j]) {
                    case B_N: 
                        v[i][j]   = 0.0;
                        u[i][j]   = -u[i][j+1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_E: 
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        v[i][j-1] = -v[i+1][j-1];
                        break;
                    case B_S:
                        v[i][j-1] = 0.0;
                        u[i][j]   = -u[i][j-1];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_W: 
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        v[i][j-1] = -v[i-1][j-1];
                        break;
                    case B_NE:
                        v[i][j]   = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j-1] = -v[i+1][j-1];
                        u[i-1][j] = -u[i-1][j+1];
                        break;
                    case B_SE:
                        v[i][j-1] = 0.0;
                        u[i][j]   = 0.0;
                        v[i][j]   = -v[i+1][j];
                        u[i-1][j] = -u[i-1][j-1];
                        break;
                    case B_SW:
                        v[i][j-1] = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j]   = -v[i-1][j];
                        u[i][j]   = -u[i][j-1];
                        break;
                    case B_NW:
                        v[i][j]   = 0.0;
                        u[i-1][j] = 0.0;
                        v[i][j-1] = -v[i-1][j-1];
                        u[i][j]   = -u[i][j+1];
                        break;
                }
            }
	    }
	    }  */

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v[0][0] = 2*vi-v[1][0];
    for (j=1;j<=jmax;j++) {
        u[0][j] = ui;
        v[0][j] = 2*vi-v[1][j];
    }
}

/* 

for(i=0;i<=(imax+1);i++){
  for(j=0;j<=(jmax+1);j++){
    if (i==0 && j>=1 && j <=jmax) {
         u[i][j] = ui;
    }
  }
}

for(i=0;i<=(imax+1);i++){
  for(j=0;j<=(jmax+1);j++){
     if (i==0 && j<=jmax) {
       v[i][j] = 2*vi-v[1][j];
     }
  }
}

 */
