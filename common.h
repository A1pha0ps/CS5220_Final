extern int fdtd_steps;
extern int topop_steps;
extern int NUMROWS;
extern int NUMCOLS;

#define pi 3.14159265358979323

#define delx 0.01                           // Cell Size
#define delt (imp0 * courant * eps0 / delx) // Time Step
#define eps0 8.8e-12
#define mu0 (4 * pi * 1e-7)
#define imp0 377.0

#define courant (sqrt(2) / 2)
#define t_1 1.0/courant + 2.0 + courant



#define	abc_c0 -(1.0 / courant - 2.0 + courant) / t_1
#define	abc_c1 -2.0 * (courant - 1.0 / courant) / t_1
#define	abc_c2 4.0 * (courant + 1.0 / courant) / t_1

#define abs(x) x < 0 ? -x : x

#define dx(i, j) dx[(i) * NUMROWS + (j)]
#define ex(i, j) ex[(i) * NUMROWS + (j)]
#define hy(i, j) hy[(i) * NUMROWS + (j)]
#define hz(i, j) hz[(i) * NUMROWS + (j)]
#define ix(i, j) ix[(i) * NUMROWS + (j)]
#define relative_eps(i, j) relative_eps[(i) * NUMROWS + j]
#define sigma(i, j) sigma[(i) * NUMROWS + j]
#define MAXLOSS 0.35

#define dxL_abc(i, j, k) dxL_abc[(k) * 6 + (j) * 3 + i]
#define dxR_abc(i, j, k) dxR_abc[(k) * 6 + (j) * 3 + i]
#define dxT_abc(i, j, k) dxT_abc[(k) * 6 + (j) * 3 + i]
#define dxB_abc(i, j, k) dxB_abc[(k) * 6 + (j) * 3 + i]

#define deltx (delt / delx)

// offset starting point for sinusoidal wave
#define IC (NUMROWS / 2)
#define JC (NUMCOLS / 2)

// cuda

#define NUMTHREADS 256
#define NUMBLKS 256

// cuda

#define NUMTHREADS 256
#define NUMBLKS 256

void reset_simulation(double *dx, double *ex, double *hy, double *hz);
void reset_bounds(double * dxL_abc, double * dxR_abc, double * dxT_ab, double * dxB_abc);


void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step,
		double *relative_eps, double *sigma, double *dxL_abc, double *dxR_abc, double *dxT_abc,
		double *dxB_abc,
		int * s_x, int * s_y, double * s_amp, double * s_off, int * s_type, int s_count);

