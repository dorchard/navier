#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include "alloc.h"
#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

int main(int argc, char *argv[])
{
    int verbose = 1;          /* Verbosity level */
    double xlength = 22.0;     /* Width of simulated domain */
    double ylength = 4.1;      /* Height of simulated domain */
    int imax = 660;           /* Number of cells horizontally */
    int jmax = 120;           /* Number of cells vertically */

    char *outname;
    int output = 0;
    int output_frequency = 0;

    double t_end = 40; //2.1       /* Simulation runtime */
    double del_t = 0.003;      /* Duration of each timestep */
    double tau = 0.5;          /* Safety factor for timestep control */

    int itermax = 100;        /* Maximum number of iterations in SOR */
    double eps = 0.001;        /* Stopping error threshold for SOR */
    double omega = 1.7;        /* Relaxation parameter for SOR */
    double gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    double Re = 150.0;         /* Reynolds number */
    double ui = 1.0;           /* Initial X velocity */
    double vi = 0.0;           /* Initial Y velocity */

    double t, delx, dely;
    int  i, j, itersor = 0, ifluid = 0, ibound = 0;
    double res;
    double **u, **v, **p, **rhs, **f, **g;
    char  **flag;
    int init_case, iters = 0;
    int show_help = 0, show_usage = 0, show_version = 0;

    if (argc > 1) {
      output = 1;
      outname = argv[1];
      output_frequency = 1;
    }

    if (argc > 2) {
      output_frequency = atoi(argv[2]);
    }

    delx = xlength/imax;
    dely = ylength/jmax;

    /* Allocate arrays */
    u    = alloc_doublematrix(imax+2, jmax+2);
    v    = alloc_doublematrix(imax+2, jmax+2);
    f    = alloc_doublematrix(imax+2, jmax+2);
    g    = alloc_doublematrix(imax+2, jmax+2);
    p    = alloc_doublematrix(imax+2, jmax+2);
    rhs  = alloc_doublematrix(imax+2, jmax+2); 
    flag = alloc_charmatrix(imax+2, jmax+2);                    

    if (!u || !v || !f || !g || !p || !rhs || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
        return 1;
    }

    unsigned long checker = 0;
    double checker1 = 0.0;

    // Set up initial values
    for (i=0;i<=imax+1;i++) {
         for (j=0;j<=jmax+1;j++) {
   	     checker += (i*jmax)+ j + 1;
	     checker1 += (i*jmax) + j + 1.0;
             u[i][j] = ui;
             v[i][j] = vi;
             p[i][j] = 0.0;
         }
     }

    init_flag(flag, imax, jmax, delx, dely, &ibound);
    printf("check flag = %ld\n", simplest_checksum_char(flag, imax, jmax));
    apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
    
    // Main loop

    for (t = 0.0; t < t_end; t += del_t, iters++) {
        printf("check u = %f\n", simplest_checksum(u, imax, jmax));
        printf("check v = %f\n", simplest_checksum(v, imax, jmax));
        set_timestep_interval(&del_t, imax, jmax, delx, dely, u, v, Re, tau);

        ifluid = (imax * jmax) - ibound;

	printf("check f = %f\n", simplest_checksum(f, imax, jmax));
        printf("check g = %f\n", simplest_checksum(g, imax, jmax));
        compute_tentative_velocity(u, v, f, g, flag, imax, jmax,
            del_t, delx, dely, gamma, Re);
	printf("check f' = %f\n", simplest_checksum(f, imax, jmax));
        printf("check g' = %f\n", simplest_checksum(g, imax, jmax));

	printf("check rhs = %f\n", simplest_checksum(rhs, imax, jmax));
        compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);
	printf("check rhs' = %f\n", simplest_checksum(rhs, imax, jmax));

        if (ifluid > 0) {
	    printf("check p = %f\n", simplest_checksum(p, imax, jmax));
            itersor = poisson(p, rhs, flag, imax, jmax, delx, dely,
                        eps, itermax, omega, &res, ifluid);
	    printf("check p' = %f\n", simplest_checksum(p, imax, jmax));
        } else {
            itersor = 0;
        }

         printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
                iters, t+del_t, del_t, itersor, res, ibound);

	
        update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);
	printf("check u' = %f\n", simplest_checksum(u, imax, jmax));
	printf("check v' = %f\n", simplest_checksum(v, imax, jmax));

        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
	printf("check u'' = %f\n", simplest_checksum(u, imax, jmax));
	printf("check v'' = %f\n", simplest_checksum(v, imax, jmax));

	if (output && (iters % output_frequency == 0)) {
	  write_ppm(u, v, p, flag, imax, jmax, xlength, ylength, outname,
		    iters, output_frequency);
	}
    }

    free_matrix(u);
    free_matrix(v);
    free_matrix(f);
    free_matrix(g);
    free_matrix(p);
    free_matrix(rhs);
    free_matrix(flag);

    return 0;
}

unsigned int simplest_checksum_char(char** in, int imax, int jmax) {
  unsigned int checksum = 0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*(i);
    }
  }
  return checksum;
}

double simplest_checksum(double** in, int imax, int jmax) {
  double checksum = 0.0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*((double)(i*jmax)+j);
    }
  }
  return checksum;
}
