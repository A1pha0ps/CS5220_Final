#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include <iostream>
using namespace std;

#include "common.h"

/*

  (row i, col j)
  arr[i,j] = arr[i * numRows + j]

*/

#define dx(i, j) dx[(i) * NUMROWS + (j)]
#define ex(i, j) ex[(i) * NUMROWS + (j)]
#define hy(i, j) hy[(i) * NUMROWS + (j)]
#define hz(i, j) hz[(i) * NUMROWS + (j)]
#define ix(i, j) ix[(i) * NUMROWS + (j)]

#define deltx (delt / delx)

// offset starting point for sinusoidal wave
#define IC (NUMROWS / 2)
#define JC (NUMCOLS / 2)

double *ix;

// test values
/*

  Material: Lossy Dielectric (Silicon @ 20 GHz)
  Relative Permitivity: 19.3
  Conductivity: 5.21 S/m
*/
double eps = 1;
double sigma = 0;

double eps_eff;

double pulse = 0.0;

double start_time, end_time;
double dtime;
double etime;
double hytime;
double hztime;

double t0 = 40.0;
int spread = 12;

void init_simulation(double *dx, double *ex, double *hy, double *hz)
{

  start_time = omp_get_wtime();

  ix = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      dx(i, j) = 0.0;
      ex(i, j) = 0.0;
      hy(i, j) = 0.0;
      hz(i, j) = 0.0;
      ix(i, j) = 0.0;
    }
  }

  eps_eff = eps + (sigma * delt) / eps0;

  end_time = omp_get_wtime();

  printf("Init Time: %f seconds\n", end_time - start_time);
}

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step)
{

  // cout << "time step" << cur_step << "\n"
  //      << endl;

  // calculate D field
  start_time = omp_get_wtime();

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {
      int iidx = (i - 1) % NUMROWS;
      int jidx = (j - 1) % NUMCOLS;
      dx(i, j) = dx(i, j) + courant * (hz(i, j) - hz(iidx, j) - hy(i, j) + hy(i, jidx));
    }
  }

  end_time = omp_get_wtime();
  dtime += (end_time - start_time);

  // Sinusoidal Source
  // 1 GHz
  pulse = 5 * exp(-.1 * (pow((t0 - cur_step) / spread, 2.0)));
  dx(IC, JC) = pulse;
  // cout << dx(IC, JC) << endl;

  start_time = omp_get_wtime();
  // Compute E-field due to permittivity
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      ex(i, j) = (1 / eps_eff) * (dx(i, j) - ix(i, j));
      ix(i, j) = ix(i, j) + ((sigma * delt) / eps0) * ex(i, j);
      // if (ex(i, j) != 0.0)
      // {
      //   cout << 1 / eps_eff << endl;
      //   cout << dx(i, j) << endl;
      //   cout << ex(i, j) << endl;
      //   cout << hy(i, j) << endl;
      //   cout << hz(i, j) << endl;
      //   cout << "\n"
      //        << endl;
      // }
    }
  }
  end_time = omp_get_wtime();
  etime += (end_time - start_time);

  start_time = omp_get_wtime();
  // Calculate Hy
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS - 1; j++)
    {
      int idx = (j + 1) % NUMCOLS;
      hy(i, j) = hy(i, j) - courant * (ex(i, idx) - ex(i, j));
      // cout << hy(i, j) << endl;
    }
  }
  end_time = omp_get_wtime();
  hytime += (end_time - start_time);

  start_time = omp_get_wtime();
  // Calculate Hz
  for (int i = 0; i < NUMROWS - 1; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      int idx = (i + 1) % NUMROWS;
      hz(i, j) = hz(i, j) + courant * (ex(idx, j) - ex(i, j));
    }
  }
  end_time = omp_get_wtime();
  hztime += (end_time - start_time);

  if (cur_step == nsteps - 1)
  {
    printf("D-field time: %f seconds\n", dtime);
    printf("E-field time: %f seconds\n", etime);
    printf("Hy-field time: %f seconds\n", hytime);
    printf("Hz-field time: %f seconds\n", hztime);
  }
}