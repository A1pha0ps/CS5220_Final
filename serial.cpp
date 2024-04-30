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

int ra = 5;
int rb = NUMROWS - ra - 1;
int ca = 5;
int cb = NUMCOLS - rb - 1;

double *ex_incident;
double *hz_incident;
double ex_incident_low0, ex_incident_low1, ex_incident_high0, ex_incident_high1;

void init_simulation(double *dx, double *ex, double *hy, double *hz)
{

  start_time = omp_get_wtime();

  ix = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  ex_incident = (double *)malloc(NUMROWS * sizeof(double));
  hz_incident = (double *)malloc(NUMROWS * sizeof(double));

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

  for (int i = 0; i < NUMROWS; i++)
  {
    ex_incident[i] = 0;
    hz_incident[i] = 0;
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

  // for (int j = 1; j < NUMROWS; j++)
  // {
  //   ex_incident[j] = ex_incident[j] + courant * (hz_incident[j - 1] - hz_incident[j]);
  // }

  // ex_incident[0] = ex_incident_low1;
  // ex_incident_low1 = ex_incident_low0;
  // ex_incident_low0 = ex_incident[1];

  // ex_incident[NUMROWS - 1] = ex_incident_high1;
  // ex_incident_high1 = ex_incident_high0;
  // ex_incident_high0 = ex_incident[NUMROWS - 2];

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {
      int iidx = (i - 1);
      int jidx = (j - 1);
      dx(i, j) = dx(i, j) + courant * (hz(i, j) - hz(iidx, j) - hy(i, j) + hy(i, jidx));
    }
  }

  end_time = omp_get_wtime();
  dtime += (end_time - start_time);

  // Sinusoidal Source
  // 1 GHz
  pulse = 500000 * exp(-.1 * (pow((t0 - cur_step) / spread, 2.0)));
  dx(IC, JC) = pulse;
  // ex_incident[20] = pulse;
  // cout << dx(IC, JC) << endl;

  /* Incident Dx values for left and right */
  // for (int i = ra; i <= rb; i++)
  // {
  //   dx(ca, i) = dx(ca, i) - courant * hz_incident[ca - 1];
  //   dx(cb, i) = dx(cb, i) + courant * hz_incident[cb];
  // }

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

  // for (int j = 0; j < NUMROWS - 1; j++)
  // {
  //   hz_incident[j] = hz_incident[j] + courant * (ex_incident[j] - ex_incident[j + 1]);
  // }

  // Calculate Hy
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS - 1; j++)
    {
      int idx = (j + 1);
      hy(i, j) = hy(i, j) - courant * (ex(i, idx) - ex(i, j));
      // cout << hy(i, j) << endl;
    }
  }
  end_time = omp_get_wtime();
  hytime += (end_time - start_time);

  /* Incident Hy values for top and bottom */
  // for (int i = ca; i <= cb; i++)
  // {
  //   hy(i, ra - 1) = hy(i, ra - 1) + courant * ex_incident[ra - 1];
  //   hy(i, rb) = hy(i, rb) - courant * ex_incident[rb];
  // }

  start_time = omp_get_wtime();
  // Calculate Hz
  for (int i = 0; i < NUMROWS - 1; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      int idx = (i + 1);
      hz(i, j) = hz(i, j) + courant * (ex(idx, j) - ex(i, j));
    }
  }
  end_time = omp_get_wtime();
  hztime += (end_time - start_time);

  /* Incident Hz values for left and right */
  // for (int j = ra; j <= rb; j++)
  // {
  //   hz(ca - 1, j) = hz(ca - 1, j) - courant * ex_incident[ca];
  //   hz(cb, j) = hz(cb, j) + courant * ex_incident[cb];
  // }

  if (cur_step == nsteps - 1)
  {
    printf("D-field time: %f seconds\n", dtime);
    printf("E-field time: %f seconds\n", etime);
    printf("Hy-field time: %f seconds\n", hytime);
    printf("Hz-field time: %f seconds\n", hztime);
  }
}
