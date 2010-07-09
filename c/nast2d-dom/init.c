#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "datadef.h"


/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void init_flag(char **flag, int imax, int jmax, double delx, double dely,
    int *ibound)
{
    int i, j;
    double mx, my, x, y, rad1;

    /* Mask of a circular obstacle */
    mx = 20.0/41.0*jmax*dely;
    my = mx;
    rad1 = 5.0/41.0*jmax*dely;
    for (i=1;i<=imax;i++) {
        for (j=1;j<=jmax;j++) {
            x = (i-0.5)*delx - mx;
            y = (j-0.5)*dely - my;
            flag[i][j] = (x*x + y*y <= rad1*rad1)?C_B:C_F;
        }
    }
    
    /* Mark the north & south boundary cells */
    for (i=0; i<=imax+1; i++) {
        flag[i][0]      = C_B;
        flag[i][jmax+1] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (j=1; j<=jmax; j++) {
        flag[0][j]      = C_B;
        flag[imax+1][j] = C_B;
    }

    /* flags for boundary cells */
    *ibound = 0;
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (!(flag[i][j] & C_F)) {
                (*ibound)++;
                if (flag[i-1][j] & C_F) flag[i][j] |= B_W;
                if (flag[i+1][j] & C_F) flag[i][j] |= B_E;
                if (flag[i][j-1] & C_F) flag[i][j] |= B_S;
                if (flag[i][j+1] & C_F) flag[i][j] |= B_N;
            }
        }
    }
}
