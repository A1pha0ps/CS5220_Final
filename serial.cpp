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
#define IC (NUMROWS / 2 - 20)
#define JC (NUMCOLS / 2 - 20)

double *ix;

// test values
/*

  Material: Lossy Dielectric (Silicon @ 20 GHz)
  Relative Permitivity: 19.3
  Conductivity: 5.21 S/m
*/
double eps = 19.3 * eps0;
double sigma = 5.21;

double eps_eff;

double pulse = 0.0;

double start_time, end_time;
double dtime;
double etime;
double hytime;
double hztime;

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

  eps_eff = eps + sigma * delt;

  end_time = omp_get_wtime();

  printf("Init Time: %f seconds\n", end_time - start_time);
}

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step)
{
  // calculate D field
  start_time = omp_get_wtime();

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {
      // cout << i << "," << j << endl;
      // cout << dx[i * NUMROWS + j] << endl;
      // cout << hz(i, j) << endl;
      dx(i, j) = dx(i, j) + deltx * (hz(i, j) - hz(i - 1, j) - hy(i, j) + hy(i, j - 1));
    }
  }

  end_time = omp_get_wtime();
  dtime += (end_time - start_time);

  // Sinusoidal Source
  // 20 GHz
  pulse = sin(2 * pi * 2 * 1e10 * delt * cur_step);
  dx(IC, JC) = pulse;

  start_time = omp_get_wtime();
  // Compute E-field due to permittivity
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      ex(i, j) = eps_eff * (dx(i, j) - ix[i, j]);
      ix(i, j) = ix(i, j) + (sigma * delt) * ex(i, j);
    }
  }
  end_time = omp_get_wtime();
  etime += (end_time - start_time);

  start_time = omp_get_wtime();
  // Calculate Hy
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      hy(i, j) = hy(i, j) - deltx * (1 / mu0) * (ex[i, j + 1] - ex[i, j]);
    }
  }
  end_time = omp_get_wtime();
  hytime += (end_time - start_time);

  start_time = omp_get_wtime();
  // Calculate Hz
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      hz(i, j) = hz(i, j) + deltx * (1 / mu0) * (ex[i + 1, j] - ex[i, j]);
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