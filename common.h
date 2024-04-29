#define NUMROWS 100
#define NUMCOLS 100
#define pi 3.14159265358979323

#define delx 0.01                     // Cell Size
#define delt (delx / (3e8 * sqrt(2))) // Time Step
#define eps0 8.8e-12
#define mu0 (4 * pi * 1e-7)

#define courant (sqrt(2) / 2)

#define nsteps 500

void init_simulation(double *dx, double *ex, double *hy, double *hz);

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step);