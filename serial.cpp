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
#define relative_eps(i, j) relative_eps[(i) * NUMROWS + j]
#define sigma(i, j) sigma[(i) * NUMROWS + j]
#define Idx(i, j) Idx[(i) * NUMROWS + (j)]
#define Ice_y(i, j) Ice_y[(i) * NUMROWS + (j)]
#define Ice_z(i, j) Ice_z[(i) * NUMROWS + (j)]

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
// double eps = 4;
// double sigma = 0.04;

double eps_eff;

double pulse;

double start_time, end_time;
double dtime;
double etime;
double hytime;
double hztime;

double t0 = 40.0;
int spread = 12;

int ra = 10;
int rb = NUMROWS - ra - 1;
int ca = 10;
int cb = NUMCOLS - ca - 1;

// variables for TF/SF
double *ex_incident;
double *hy_incident;
double ex_incident_low0, ex_incident_low1, ex_incident_high0, ex_incident_high1;

// variables for PML

int pmlLen = 10;

double *sigy;
double *sigz;

double *Idx;
double *Ice_y;
double *Ice_z;

// describes properties of grid material
double *relative_eps, *sigma;

void init_simulation(double *dx, double *ex, double *hy, double *hz)
{

  start_time = omp_get_wtime();

  ix = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  ex_incident = (double *)malloc(NUMROWS * sizeof(double));
  hy_incident = (double *)malloc(NUMROWS * sizeof(double));

  relative_eps = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  sigma = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  Idx = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  Ice_y = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  Ice_z = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  sigy = (double *)malloc(NUMCOLS * sizeof(double));
  sigz = (double *)malloc(NUMROWS * sizeof(double));

  for (int i = 0; i < NUMROWS; i++)
  {
    sigz[i] = 0.0;
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    sigy[i] = 0.0;
  }

  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      dx(i, j) = 0.0;
      ex(i, j) = 0.0;
      hy(i, j) = 0.0;
      hz(i, j) = 0.0;
      ix(i, j) = 0.0;
      Idx(i, j) = 0.0;
      Ice_y(i, j) = 0.0;
      Ice_z(i, j) = 0.0;

      // sigy(i, j) = (eps0 / (2.0 * delt)) * pow((j / yLenPML), 3);
      // sigz(i, j) = (eps0 / (2.0 * delt)) * pow((i / zLenPML), 3);

      /*
        Code to specify properties of cell
      */
      if (i > (NUMROWS / 4) && i < (3 * NUMROWS / 4) && j > NUMCOLS / 4 && j < (3 * NUMCOLS / 4))
      {
        relative_eps(i, j) = 4;
        sigma(i, j) = 0.1;
      }
      else
      {
        relative_eps(i, j) = 1;
        sigma(i, j) = 0;
      }

      // printf("%1.2f ", relative_eps(i, j));
    }
    // printf("\n");
  }

  for (int i = 0; i < pmlLen; i++)
  {
    sigy[i] = (eps0 / (2.0 * delt)) * pow((i / pmlLen), 3);
    sigy[NUMCOLS - i - 1] = (eps0 / (2.0 * delt)) * pow((i / pmlLen), 3);

    sigz[i] = (eps0 / (2.0 * delt)) * pow((i / pmlLen), 3);
    sigy[NUMROWS - i - 1] = (eps0 / (2.0 * delt)) * pow((i / pmlLen), 3);
  }

  for (int i = 0; i < NUMROWS; i++)
  {
    ex_incident[i] = 0;
    hy_incident[i] = 0;
  }

  // eps_eff = eps + (sigma * delt) / eps0;

  end_time = omp_get_wtime();

  printf("Init Time: %f seconds\n", end_time - start_time);
}

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step)
{

  // cout << "time step" << cur_step << "\n"
  //      << endl;

  // calculate D field
  start_time = omp_get_wtime();

  for (int j = 1; j < NUMROWS; j++)
  {
    ex_incident[j] = ex_incident[j] + courant * (hy_incident[j - 1] - hy_incident[j]);
  }

  ex_incident[0] = ex_incident_low1;
  ex_incident_low1 = ex_incident_low0;
  ex_incident_low0 = ex_incident[1];

  ex_incident[NUMROWS - 1] = ex_incident_high1;
  ex_incident_high1 = ex_incident_high0;
  ex_incident_high0 = ex_incident[NUMROWS - 2];

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {
      // int iidx = (i - 1);
      // int jidx = (j - 1);

      double curl_h = hz(i, j) - hz(i, j - 1) - hy(i, j) + hy(i - 1, j);
      double c1 = (2 / delt * (sigy[j] + sigz[i]));
      double c2 = (sigy[j] * sigz[i] * delt) / pow(eps0, 2);
      double c3 = c2 * delt;

      Idx(i, j) = Idx(i, j) + dx(i, j);

      dx(i, j) = (1 / (4 + c1 + c3)) * (dx(i, j) * (4 - c1 - c3) + 4 * courant * curl_h - (c2 / 4) * Idx(i, j));
    }
  }

  end_time = omp_get_wtime();
  dtime += (end_time - start_time);

  // Gaussian Pulse
  pulse = 5 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

  // Sinusoidal Source
  // 20 GHz
  // pulse = 10 * sin(2 * pi * 7 * 1e8 * delt * cur_step);

  // dx(IC - 20, JC) += pulse;
  ex_incident[5] = pulse;
  // cout << dx(IC, JC) << endl;

  /* Incident Dx values for left and right */
  for (int i = ca; i <= cb; i++)
  {
    dx(ra, i) = dx(ra, i) + courant * hy_incident[ra - 1];
    dx(rb, i) = dx(rb, i) - courant * hy_incident[rb];
  }

  start_time = omp_get_wtime();
  // Compute E-field due to permittivity
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      double eps_eff = relative_eps(i, j) + (sigma(i, j) * delt) / eps0;
      ex(i, j) = (1 / eps_eff) * (dx(i, j) - ix(i, j));
      ix(i, j) = ix(i, j) + (sigma(i, j) * delt) / eps0 * ex(i, j);
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

  for (int j = 0; j < NUMROWS - 1; j++)
  {
    hy_incident[j] = hy_incident[j] + courant * (ex_incident[j] - ex_incident[j + 1]);
  }

  // Calculate Hy
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS - 1; j++)
    {
      double curl_e = ex(i + 1, j) - ex(i, j);

      Ice_z(i, j) += curl_e;

      double c1 = (1 - (sigz[i] * delt) / eps0);

      hy(i, j) = hy(i, j) * c1 - (sigy[j] / eps0) * courant * Ice_z(i, j) - courant * curl_e;
      // cout << hy(i, j) << endl;
    }
  }
  end_time = omp_get_wtime();
  hytime += (end_time - start_time);

  /* Incident Hy values for top and bottom */
  for (int i = ca; i <= cb; i++)
  {
    hy(ra - 1, i) = hy(ra - 1, i) + courant * ex_incident[ra];
    hy(rb, i) = hy(rb, i) - courant * ex_incident[rb];
  }

  start_time = omp_get_wtime();
  // Calculate Hz
  for (int i = 0; i < NUMROWS - 1; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      // int idx = (i + 1);

      double curl_e = ex(i, j + 1) - ex(i, j);
      Ice_y(i, j) = Ice_y(i, j) + curl_e;

      double c1 = (1 - (sigy[j] * delt) / eps0);

      hz(i, j) = courant * curl_e + (sigz[i] / eps0) * courant * Ice_y(i, j) + hz(i, j) * c1;
    }
  }
  end_time = omp_get_wtime();
  hztime += (end_time - start_time);

  /* Incident Hz values for left and right */
  for (int j = ra; j <= rb; j++)
  {
    hz(j, ca - 1) = hz(j, ca - 1) - courant * ex_incident[j];
    hz(j, cb) = hz(j, cb) + courant * ex_incident[j];
  }

  if (cur_step == nsteps - 1)
  {
    printf("D-field time: %f seconds\n", dtime);
    printf("E-field time: %f seconds\n", etime);
    printf("Hy-field time: %f seconds\n", hytime);
    printf("Hz-field time: %f seconds\n", hztime);
  }
}