#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <sys/stat.h>
#include "alloc.h"
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))


void calc_psi_zeta(double **u, double **v, double **psi, double **zeta,
		   char **flag, int imax, int jmax, double delx, double dely);

void write_ppm(double** u, double** v, double** p, char** flag,
	      int imax, int jmax, double xlength, double ylength, 
	       char* outname, int iters, int freq) {
  int i, j;
  double zmax = -1e10, zmin = 1e10;
  double pmax = -1e10, pmin = 1e10;
  double umax = -1e10, umin = 1e10;
  double vmax = -1e10, vmin = 1e10;

  char outfile[64];
  sprintf(outfile, "%s/%06d.ppm", outname, (iters/freq));
  
  mode_t process_mask = umask(0);
  mkdir(outname, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);

  FILE *fout = fopen(outfile, "wb");
  if (!fout) {
    fprintf (stderr, "Could not open '%s'\n", outfile);
    return;
  }
  
  double delx = xlength/imax;
  double dely = ylength/jmax;

  int outmode = 0;
  double** psi;
  double** zeta;

  if (outmode == 0 || outmode == 1) {
    if (outmode == 1) {
      psi  = alloc_doublematrix(imax+2, jmax+2);
    } else {
      zeta = alloc_doublematrix(imax+2, jmax+2);
    }
    calc_psi_zeta(u, v, psi, zeta, flag, imax, jmax, delx, dely);
  } else {
     for (j = 1; j < jmax+1 ; j++) {
         for (i = 1; i < imax+1 ; i++) {
            if (flag[i][j] & C_F) {
	      //zmax = max(zmax, zeta[i][j]);
	      //zmin = min(zmin, zeta[i][j]);
	      pmax = max(pmax, p[i][j]);
	      pmin = min(pmin, p[i][j]);
	      //umax = max(umax, u[i][j]);
	      //umin = min(umin, u[i][j]);
	      //vmax = max(vmax, v[i][j]);
	      //vmin = min(vmin, v[i][j]);
	    }
	 }
     }
  }

  fprintf(fout, "P6 %d %d 255\n", imax, jmax);

  for (j = 1; j < jmax+1 ; j++) {
      for (i = 1; i < imax+1 ; i++) {
          int r, g, b;
          if (!(flag[i][j] & C_F)) {
                r = 0; b = 0; g = 255;
          } else {
	    //zmax = max(zmax, zeta[i][j]);
	    // zmin = min(zmin, zeta[i][j]);
	    //  pmax = max(pmax, psi[i][j]);
	    //  pmin = min(pmin, psi[i][j]);
	    //  umax = max(umax, u[i][j]);
	    //  umin = min(umin, u[i][j]);
	    //  vmax = max(vmax, v[i][j]);
	    //  vmin = min(vmin, v[i][j]);
                //// set rgb based on zeta
		if (outmode == 0) {
                    double z = (i < imax && j < jmax)?zeta[i][j]:0.0;
                    r = g = b = pow(fabs(z/12.6),.4) * 255;
                } else if (outmode == 1) {
                    double p = (i < imax && j < jmax)?psi[i][j]:0.0;
                    r = g = b = (p+3.0)/7.5 * 255; 
                } else if (outmode == 2) {
		  r = g = b = (p[i][j]-pmin)/(pmax-pmin) * 255;
		}
                //double z = (i < imax && j < jmax)?zeta[i][j]:0.0;
                //r = g = b = pow(fabs(z/12.6),.4) * 255;
            }
            fprintf(fout, "%c%c%c", r, g, b);
     }
 }

 fclose(fout);
 
 if (outmode==1) {
   free_matrix(psi);
 }
 if (outmode==0) {
   free_matrix(zeta);
 } 
}

/* Computation of stream function and vorticity */
void calc_psi_zeta(double **u, double **v, double **psi, double **zeta,
    char **flag, int imax, int jmax, double delx, double dely)
{
    int i, j;

    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
    for (i=1;i<=imax-1;i++) {
        for (j=1;j<=jmax-1;j++) {
            if ( (flag[i][j] & C_F) && (flag[i+1][j] & C_F) &&
                (flag[i][j+1] & C_F) && (flag[i+1][j+1] & C_F)) {
                zeta[i][j] = (u[i][j+1]-u[i][j])/dely
                            -(v[i+1][j]-v[i][j])/delx;
            } else {
                zeta[i][j] = 0.0;
            }
        }
    }

    /* Computation of the stream function at the upper right corner    */
    /* of cell (i,j) (only if bother lower cells are fluid cells)      */
    //for (i=0;i<=imax;i++) {
    //    psi[i][0] = 0.0;
    ///    for (j=1;j<=jmax;j++) {
    //        psi[i][j] = psi[i][j-1];
    //        if ((flag[i][j] & C_F) || (flag[i+1][j] & C_F)) {
    //           psi[i][j] += u[i][j]*dely;
    //      }
    //  }
    //}
}


