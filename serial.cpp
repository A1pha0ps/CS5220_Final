#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include <iostream>
using namespace std;

#include "common.h"

/*

  (row i, col j)
  (y, z)
  arr[i,j] = arr[i * numRows + j]

*/

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

// variables for TF/SF
int ra = 10;
int rb = NUMCOLS - ra - 1;
int ca = 10;
int cb = NUMROWS - ca - 1;

int NLOSSY1D = 20;
double *ex_incident;
double *hz_incident;
double *hz_inc_c1, *hz_inc_c2, *ex_inc_c1, *ex_inc_c2;

// second order ABC

double abc_c0, abc_c1, abc_c2;
double *dxL_abc, *dxR_abc, *dxT_abc, *dxB_abc;

// describes properties of grid material
double *relative_eps, *sigma;

void init_abc()
{
  dxL_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
  dxR_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
  dxT_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));
  dxB_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));

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
}

void init_tfsf()
{
  ex_incident = (double *)malloc(NUMCOLS * sizeof(double));
  hz_incident = (double *)malloc(NUMCOLS * sizeof(double));
  hz_inc_c1 = (double *)malloc(NUMCOLS * sizeof(double));
  hz_inc_c2 = (double *)malloc(NUMCOLS * sizeof(double));
  ex_inc_c1 = (double *)malloc(NUMCOLS * sizeof(double));
  ex_inc_c2 = (double *)malloc(NUMCOLS * sizeof(double));

  for (int i = 0; i < NUMCOLS; i++)
  {
    ex_incident[i] = 0.0;
    hz_incident[i] = 0.0;
    hz_inc_c1[i] = 0.0;
    hz_inc_c2[i] = 0.0;
    ex_inc_c1[i] = 0.0;
    ex_inc_c2[i] = 0.0;
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    if (i < NUMCOLS - NLOSSY1D - 1)
    {
      hz_inc_c1[i] = 1.0;
      hz_inc_c2[i] = courant / imp0;
      ex_inc_c1[i] = 1.0;
      ex_inc_c2[i] = courant * imp0;
    }
    else
    {
      double temp = i - (NUMCOLS - 1 - NLOSSY1D) + 0.5;
      double lossFactor = MAXLOSS * pow(temp / NLOSSY1D, 2);
      ex_inc_c1[i] = (1.0 - lossFactor) / (1.0 + lossFactor);
      ex_inc_c2[i] = courant * imp0 / (1.0 + lossFactor);
      temp += 0.5;
      lossFactor = MAXLOSS * pow(temp / NLOSSY1D, 2);
      hz_inc_c1[i] = (1.0 - lossFactor) / (1.0 + lossFactor);
      hz_inc_c2[i] = courant / imp0 / (1.0 + lossFactor);
    }
  }
}

void init_simulation(double *dx, double *ex, double *hy, double *hz)
{

  start_time = omp_get_wtime();

  ix = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  relative_eps = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  sigma = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  // Idx = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  // Ice_y = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  // Ice_z = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      dx(i, j) = 0.0;
      ex(i, j) = 0.0;
      hy(i, j) = 0.0;
      hz(i, j) = 0.0;
      ix(i, j) = 0.0;

      /*
        Code to specify properties of cell
      */
      if (i > (NUMROWS / 2 - 20) && i < (NUMROWS / 2 + 20) && j > (NUMCOLS / 2 - 20) && j < (NUMCOLS / 2 + 20))
      {
        relative_eps(i, j) = 19.3;
        sigma(i, j) = 0.04;
      }
      else
      {
        relative_eps(i, j) = 1;
        sigma(i, j) = 0;
      }
    }
  }

  // init_tfsf();

  init_abc();

  end_time = omp_get_wtime();

  printf("Init Time: %f seconds\n", end_time - start_time);
}

void updateHFields(double *dx, double *ex, double *hy, double *hz, int cur_step)
{
  // Update H-fields
  // Calculate Hy
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS - 1; j++)
    {
      double curl_e = ex(i, j + 1) - ex(i, j);
      hy(i, j) = hy(i, j) - courant / imp0 * curl_e;
    }
  }

  // Calculate Hz
  for (int i = 0; i < NUMROWS - 1; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      double curl_e = ex(i + 1, j) - ex(i, j);
      hz(i, j) = hz(i, j) + courant / imp0 * curl_e;
    }
  }
}

void updateEFields(double *dx, double *ex, double *hy, double *hz, int cur_step)
{
  // calculate D field

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {

      double curl_h = hz(i, j) - hz(i - 1, j) - hy(i, j) + hy(i, j - 1);

      double loss_factor_cond = sigma(i, j) * delt / (2 * relative_eps(i, j));

      ex(i, j) = ((1.0 - loss_factor_cond) / (1.0 + loss_factor_cond)) * ex(i, j) +
                 courant * imp0 / relative_eps(i, j) / (1.0 + loss_factor_cond) * curl_h;
    }
  }

  // Compute E-field due to permittivity
  // for (int i = 0; i < NUMROWS; i++)
  // {
  //   for (int j = 0; j < NUMCOLS; j++)
  //   {
  //     double eps_eff = relative_eps(i, j) + (sigma(i, j) * delt) / eps0;
  //     ex(i, j) = (1 / eps_eff) * (dx(i, j) - ix(i, j));
  //     ix(i, j) = ix(i, j) + (sigma(i, j) * delt) / eps0 * ex(i, j);
  //   }
  // }
}

void updateTfSf(double *dx, double *ex, double *hy, double *hz, int cur_step)
{

  // update Hz left and right
  for (int i = ca; i <= cb; i++)
  {
    hz(ra - 1, i) -= courant / imp0 * ex_incident[ra];
    hz(rb, i) += courant / imp0 * ex_incident[rb];
  }

  // update Hy
  for (int i = ra; i <= rb; i++)
  {
    hy(i, ca - 1) += courant / imp0 * ex_incident[i];
    hy(i, cb) -= courant / imp0 * ex_incident[i];
  }

  // update hz incident
  for (int i = 0; i < NUMCOLS - 1; i++)
  {
    hz_incident[i] = hz_inc_c1[i] * hz_incident[i] + hz_inc_c2[i] * (ex_incident[i + 1] - ex_incident[i]);
  }

  // update ex incident
  for (int i = 1; i < NUMCOLS; i++)
  {
    ex_incident[i] = ex_inc_c1[i] * ex_incident[i] + ex_inc_c2[i] * (hz_incident[i] - hz_incident[i - 1]);
  }

  // set source node

  // Gaussian Pulse
  // pulse = 5 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

  // Sinusoidal Source
  // 20 GHz
  // pulse = 10 * sin(2 * pi * 7 * 1e8 * delt * cur_step);

  // dx(IC - 20, JC) += pulse;

  ex_incident[0] = 5 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

  for (int i = ca; i <= cb; i++)
  {
    dx(ra, i) -= courant * imp0 * hz_incident[ra - 1];
    dx(rb, i) += courant * imp0 * hz_incident[rb];
  }
}

void apply_abc(double *dx)
{
  for (int i = 0; i < NUMROWS; i++)
  {
    dx(0, i) = abc_c0 * (dx(2, i) + dxL_abc(0, 1, i)) +
               abc_c1 * (dxL_abc(0, 0, i) + dxL_abc(2, 0, i) - dx(1, i) -
                         dxL_abc(1, 1, i)) +
               abc_c2 * dxL_abc(1, 0, i) - dxL_abc(2, 1, i);

    // dx(i, 0) = abc_c0 * (dx(i, 2) + dxL_abc(0, 1, i)) +
    //            abc_c1 * (dxL_abc(0, 0, i) + dxL_abc(2, 0, i) - dx(i, 1) -
    //                      dxL_abc(1, 1, i)) +
    //            abc_c2 * dxL_abc(1, 0, i) - dxL_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxL_abc(j, 1, i) = dxL_abc(j, 0, i);
      dxL_abc(j, 0, i) = dx(j, i);
      // dxL_abc(j, 0, i) = dx(i, j);
    }
  }

  for (int i = 0; i < NUMROWS; i++)
  {
    dx(NUMCOLS - 1, i) = abc_c0 * (dx(NUMCOLS - 3, i) + dxR_abc(0, 1, i)) +
                         abc_c1 * (dxR_abc(0, 0, i) + dxR_abc(2, 0, i) - dx(NUMCOLS - 2, i) -
                                   dxR_abc(1, 1, i)) +
                         abc_c2 * dxR_abc(1, 0, i) - dxR_abc(2, 1, i);
    // dx(i, NUMCOLS - 1) = abc_c0 * (dx(i, NUMCOLS - 3) + dxR_abc(0, 1, i)) +
    //                      abc_c1 * (dxR_abc(0, 0, i) + dxR_abc(2, 0, i) - dx(i, NUMCOLS - 2) -
    //                                dxR_abc(1, 1, i)) +
    //                      abc_c2 * dxR_abc(1, 0, i) - dxR_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxR_abc(j, 1, i) = dxR_abc(j, 0, i);
      dxR_abc(j, 0, i) = dx(NUMCOLS - 1 - j, i);
      //  dxR_abc(j, 0, i) = dx(i, NUMCOLS - 1 - j);
    }
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    dx(i, 0) = abc_c0 * (dx(i, 2) + dxB_abc(0, 1, i)) +
               abc_c1 * (dxB_abc(0, 0, i) + dxB_abc(2, 0, i) - dx(i, 1) -
                         dxB_abc(1, 1, i)) +
               abc_c2 * dxB_abc(1, 0, i) - dxB_abc(2, 1, i);
    // dx(0, i) = abc_c0 * (dx(2, i) + dxB_abc(0, 1, i)) +
    //            abc_c1 * (dxB_abc(0, 0, i) + dxB_abc(2, 0, i) - dx(1, i) -
    //                      dxB_abc(1, 1, i)) +
    //            abc_c2 * dxB_abc(1, 0, i) - dxB_abc(2, 1, i);
    for (int j = 0; j < 3; j++)
    {
      dxB_abc(j, 1, i) = dxB_abc(j, 0, i);
      dxB_abc(j, 0, i) = dx(i, j);
      // dxB_abc(j, 0, i) = dx(j, i);
    }
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    dx(i, NUMROWS - 1) = abc_c0 * (dx(i, NUMROWS - 3) + dxT_abc(0, 1, i)) +
                         abc_c1 * (dxT_abc(0, 0, i) + dxT_abc(2, 0, i) - dx(i, NUMROWS - 2) -
                                   dxT_abc(1, 1, i)) +
                         abc_c2 * dxT_abc(1, 0, i) - dxT_abc(2, 1, i);
    // dx(NUMROWS - 1, i) = abc_c0 * (dx(NUMROWS - 3, i) + dxT_abc(0, 1, i)) +
    //                      abc_c1 * (dxT_abc(0, 0, i) + dxT_abc(2, 0, i) - dx(NUMROWS - 2, i) -
    //                                dxT_abc(1, 1, i)) +
    //                      abc_c2 * dxT_abc(1, 0, i) - dxT_abc(2, 1, i);
    for (int j = 0; j < 3; j++)
    {
      dxT_abc(j, 1, i) = dxT_abc(j, 0, i);
      dxT_abc(j, 0, i) = dx(i, NUMROWS - 1 - j);
      // dxT_abc(j, 0, i) = dx(NUMROWS - 1 - j, i);
    }
  }
}

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step)
{
  updateHFields(dx, ex, hy, hz, cur_step);

  // updateTfSf(dx, ex, hy, hz, cur_step);

  updateEFields(dx, ex, hy, hz, cur_step);

  pulse = 10 * sin(2 * pi * 7 * 1e8 * delt * cur_step);

  ex(NUMROWS / 2 - 40, NUMCOLS / 2) = pulse;

  // ex(NUMROWS / 2 + 40, NUMCOLS / 2) = pulse;

  apply_abc(ex);

  // ex(NUMROWS / 2 - 30, NUMCOLS / 2) = pulse = 5 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

  // if (cur_step == nsteps - 1)
  // {
  //   printf("D-field time: %f seconds\n", dtime);
  //   printf("E-field time: %f seconds\n", etime);
  //   printf("Hy-field time: %f seconds\n", hytime);
  //   printf("Hz-field time: %f seconds\n", hztime);
  // }
}