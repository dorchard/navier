void compute_tentative_velocity(double **u, double **v, double **f, double **g,
    char **flag, int imax, int jmax, double del_t, double delx, double dely,
    double gamma, double Re);

void compute_rhs(double **f, double **g, double **rhs, char **flag, int imax,
    int jmax, double del_t, double delx, double dely);

int poisson(double **p, double **rhs, char **flag, int imax, int jmax,
    double delx, double dely, double eps, int itermax, double omega,
    double *res, int ifull);

void update_velocity(double **u, double **v, double **f, double **g, double **p,
    char **flag, int imax, int jmax, double del_t, double delx, double dely);

void set_timestep_interval(double *del_t, int imax, int jmax, double delx,
    double dely, double **u, double **v, double Re, double tau);
