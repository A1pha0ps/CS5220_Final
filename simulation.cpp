#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <omp.h>

#include <iostream>
using namespace std;

#include "common.h"

double abc_c0, abc_c1, abc_c2;

//define design and target space

int * target_x;
int * target_y;

int * design_x;
int * design_y;
//Silicon for Device
double t_perm = 19.3;
double t_cond = 5.21;
//Copper for Shield
double d_perm = 1;
double d_cond = 5.8e7;


void init_abc(double *dxL_abc, double *dxR_abc, double *dxT_abc, double *dxB_abc)
{
  for (int i = 0; i < NUMROWS * 6; i++)
  {
    dxL_abc[i] = 0;
    dxR_abc[i] = 0;
  }

  for (int i = 0; i < NUMCOLS * 6; i++)
  {
    dxT_abc[i] = 0;
    dxB_abc[i] = 0;
  }

  double t0 = courant;
  double t1 = 1.0 / t0 + 2.0 + t0;
  abc_c0 = -(1.0 / t0 - 2.0 + t0) / t1;
  abc_c1 = -2.0 * (t0 - 1.0 / t0) / t1;
  abc_c2 = 4.0 * (t0 + 1.0 / t0) / t1;

  cout << "finished init_abc" << endl;
}

void init_mat_properties(double *relative_eps, double *sigma)
{
  /*
    Code to specify properties of cell

    Reference:
    Material: Lossy Dielectric (Silicon @ 20 GHz)
    Relative Permitivity: 19.3
    Conductivity: 5.21 S/m
  */

  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {

      bool cond1 = (i >= 100 && i <= 102 && j <= 60);
      bool cond2 = (i >= 100 && i <= 102 && j >= 90 && j <= 110);
      bool cond3 = (i >= 100 && i <= 102 && j >= 140);

      if (cond1 || cond2 || cond3)
      {
        relative_eps(i, j) = 2.5;
        sigma(i, j) = 8e6;
      }
      else
      {
        relative_eps(i, j) = 1;
        sigma(i, j) = 0;
      }
    }
  }

  cout << "finished init mat properties" << endl;
}

void free_all(double *dx, double *ex, double *hy, double *hz, double *dxL_abc,
              double *dxR_abc, double *dxT_abc, double *dxB_abc,
              double *relative_eps, double *sigma)
{
  free(dx);
  free(ex);
  free(hy);
  free(hz);
  free(dxL_abc);
  free(dxR_abc);
  free(dxT_abc);
  free(dxB_abc);
  free(relative_eps);
  free(sigma);
}


int main()
{

  // fundamental fields
  double *dx, *ex, *hy, *hz;

  // abcs
  double *dxL_abc, *dxR_abc, *dxT_abc, *dxB_abc;

  // describes properties of grid material
  double *relative_eps, *sigma;

  // double *ix;

  dx = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  ex = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  hy = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  hz = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  // ix = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  relative_eps = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  sigma = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  dxL_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
  dxR_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
  dxT_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));
  dxB_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));

  init_simulation(dx, ex, hy, hz);
  init_abc(dxL_abc, dxR_abc, dxT_abc, dxB_abc);
  init_mat_properties(relative_eps, sigma);

  FILE *fp;

  fp = fopen("ex_field.txt", "w");

  for (int step = 0; step < nsteps; step++)
  {
    // cout << step << endl;
    simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma,
                       abc_c0, abc_c1, abc_c2, dxL_abc, dxR_abc, dxT_abc, dxB_abc);

    /* Write the E field out to a file "Ez" */

    // fprintf(fp, "%d\n", step);

    // if (step % 1 == 0)
    // {
    for (int i = 0; i < NUMROWS; i++)
    {
      for (int j = 0; j < NUMCOLS; j++)
      {
        fprintf(fp, "%6.3f ", ex[i * NUMROWS + j]);
      }
      fprintf(fp, " \n");
    }

    fprintf(fp, "\n");
    // }
  }

  fclose(fp);

  free_all(dx, ex, hy, hz, dxL_abc, dxR_abc, dxT_abc, dxB_abc, relative_eps, sigma);

  return 0;
}
