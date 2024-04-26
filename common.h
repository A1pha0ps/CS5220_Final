#define NUMROWS 10000
#define NUMCOLS 10000
#define pi 3.14159265358979323

#define delx 0.0015       // Cell Size
#define delt (delx / 6e8) // Time Step
#define eps0 8.8e-12;
#define mu0 (4 * pi * 1e-7)

#define nsteps 1000

void init_simulation(double *dx, double *ex, double *hy, double *hz);

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step);