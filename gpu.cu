#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <cuda.h>
#include <iostream>
using namespace std;

#include "common.h"

extern int nsteps, NUMROWS, NUMCOLS;

/*

  (row i, col j)
  (y, z)
  arr[i,j] = arr[i * numRows + j]

*/
double pulse;

// double start_time, end_time;
double dtime;
double etime;
double hytime;
double hztime;

double t0 = 40.0;
int spread = 15;

// variables for TF/SF
int ra = 10;
int rb = NUMCOLS - ra - 1;
int ca = 10;
int cb = NUMROWS - ca - 1;

int NLOSSY1D = 20;
double *ex_incident;
double *hz_incident;
double *hz_inc_c1, *hz_inc_c2, *ex_inc_c1, *ex_inc_c2;

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

  // start_time = omp_get_wtime();
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
    }
  }

  // init_tfsf();

  // init_abc();

  // end_time = omp_get_wtime();

  // printf("Init Time: %f seconds\n", end_time - start_time);
}

__global__ void updateHyGPU(double *ex, double *hy, double *hz)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMROWS * NUMCOLS; i += stride)
  {
    int rowIdx = i / NUMCOLS;
    int colIdx = i / NUMROWS;

    if (colIdx == NUMCOLS - 1)
    {
      continue;
    }

    double curl_e = ex(rowIdx, colIdx + 1) - ex(rowIdx, colIdx);
    hy(rowIdx, colIdx) = hy(rowIdx, colIdx) - courant / imp0 * curl_e;
  }
}

__global__ void updateHzGPU(double *ex, double *hy, double *hz)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMROWS * NUMCOLS; i += stride)
  {
    int rowIdx = i / NUMCOLS;
    int colIdx = i / NUMROWS;

    if (colIdx == NUMROWS - 1)
    {
      continue;
    }

    double curl_e = ex(rowIdx + 1, colIdx) - ex(rowIdx, colIdx);
    hz(rowIdx, colIdx) = hz(rowIdx, colIdx) + courant / imp0 * curl_e;
  }
}

__global__ void updateEzGPU(double *ex, double *hy, double *hz, double *relative_eps, double *sigma)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMROWS * NUMCOLS; i += stride)
  {
    int rowIdx = i / NUMCOLS;
    int colIdx = i / NUMROWS;

    if (colIdx == 0 || rowIdx == 0)
    {
      continue;
    }

    double curl_h = hz(rowIdx, colIdx) - hz(rowIdx - 1, colIdx) - hy(rowIdx, colIdx) + hy(rowIdx, colIdx - 1);

    double loss_factor_cond = sigma(rowIdx, colIdx) * delt / (2 * relative_eps(rowIdx, colIdx));

    ex(rowIdx, colIdx) = ((1.0 - loss_factor_cond) / (1.0 + loss_factor_cond)) * ex(rowIdx, colIdx) +
                         courant * imp0 / relative_eps(rowIdx, colIdx) / (1.0 + loss_factor_cond) * curl_h;
  }
}

/*

void updateHFields(double *dx, double *ex, double *hy, double *hz)
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

*/

/*
void updateEFields(double *dx, double *ex, double *hy, double *hz, int cur_step, double *relative_eps, double *sigma)
{
  // calculate D field

  for (int i = 1; i < NUMROWS; i++)
  {
    for (int j = 1; j < NUMCOLS; j++)
    {

      double curl_h = hz(i, j) - hz(i - 1, j) - hy(i, j) + hy(i, j - 1);

      double loss_factor_cond = sigma(i, j) * delt / (2 * relative_eps(i, j));

      // cout << loss_factor_cond << endl;

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
*/

/*
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

  ex_incident[0] = 50000 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

  for (int i = ca; i <= cb; i++)
  {
    dx(ra, i) -= courant * imp0 * hz_incident[ra - 1];
    dx(rb, i) += courant * imp0 * hz_incident[rb];
  }
}
*/

__global__ void apply_abc_left_col(double *ex, double abc_c0, double abc_c1, double abc_c2, double *dxL_abc, int NUMROWS, int NUMCOLS)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMROWS; i += stride)
  {
    ex(0, i) = abc_c0 * (ex(2, i) + dxL_abc(0, 1, i)) +
               abc_c1 * (dxL_abc(0, 0, i) + dxL_abc(2, 0, i) - ex(1, i) -
                         dxL_abc(1, 1, i)) +
               abc_c2 * dxL_abc(1, 0, i) - dxL_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxL_abc(j, 1, i) = dxL_abc(j, 0, i);
      dxL_abc(j, 0, i) = ex(j, i);
    }
  }
}

__global__ void apply_abc_right_col(double *ex, double abc_c0, double abc_c1, double abc_c2, double *dxR_abc, int NUMROWS, int NUMCOLS)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMROWS; i += stride)
  {
    ex(NUMCOLS - 1, i) = abc_c0 * (ex(NUMCOLS - 3, i) + dxR_abc(0, 1, i)) +
                         abc_c1 * (dxR_abc(0, 0, i) + dxR_abc(2, 0, i) - ex(NUMCOLS - 2, i) -
                                   dxR_abc(1, 1, i)) +
                         abc_c2 * dxR_abc(1, 0, i) - dxR_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxR_abc(j, 1, i) = dxR_abc(j, 0, i);
      dxR_abc(j, 0, i) = ex(NUMCOLS - 1 - j, i);
    }
  }
}

__global__ void apply_abc_top_row(double *ex, double abc_c0, double abc_c1, double abc_c2, double *dxB_abc, int NUMROWS, int NUMCOLS)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMCOLS; i += stride)
  {
    ex(i, 0) = abc_c0 * (ex(i, 2) + dxB_abc(0, 1, i)) +
               abc_c1 * (dxB_abc(0, 0, i) + dxB_abc(2, 0, i) - ex(i, 1) -
                         dxB_abc(1, 1, i)) +
               abc_c2 * dxB_abc(1, 0, i) - dxB_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxB_abc(j, 1, i) = dxB_abc(j, 0, i);
      dxB_abc(j, 0, i) = ex(i, j);
    }
  }
}

__global__ void apply_abc_bot_row(double *ex, double abc_c0, double abc_c1, double abc_c2, double *dxT_abc, int NUMROWS, int NUMCOLS)
{
  // Get thread ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = tid; i < NUMCOLS; i += stride)
  {
    ex(i, NUMROWS - 1) = abc_c0 * (ex(i, NUMROWS - 3) + dxT_abc(0, 1, i)) +
                         abc_c1 * (dxT_abc(0, 0, i) + dxT_abc(2, 0, i) - ex(i, NUMROWS - 2) -
                                   dxT_abc(1, 1, i)) +
                         abc_c2 * dxT_abc(1, 0, i) - dxT_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxT_abc(j, 1, i) = dxT_abc(j, 0, i);
      dxT_abc(j, 0, i) = ex(i, NUMROWS - 1 - j);
    }
  }
}

/*

void apply_abc(double *ex, double abc_c0, double abc_c1, double abc_c2,
               double *dxL_abc, double *dxR_abc, double *dxT_abc, double *dxB_abc)
{

  for (int i = 0; i < NUMROWS; i++)
  {
    ex(0, i) = abc_c0 * (ex(2, i) + dxL_abc(0, 1, i)) +
               abc_c1 * (dxL_abc(0, 0, i) + dxL_abc(2, 0, i) - ex(1, i) -
                         dxL_abc(1, 1, i)) +
               abc_c2 * dxL_abc(1, 0, i) - dxL_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxL_abc(j, 1, i) = dxL_abc(j, 0, i);
      dxL_abc(j, 0, i) = ex(j, i);
    }
  }

  for (int i = 0; i < NUMROWS; i++)
  {
    ex(NUMCOLS - 1, i) = abc_c0 * (ex(NUMCOLS - 3, i) + dxR_abc(0, 1, i)) +
                         abc_c1 * (dxR_abc(0, 0, i) + dxR_abc(2, 0, i) - ex(NUMCOLS - 2, i) -
                                   dxR_abc(1, 1, i)) +
                         abc_c2 * dxR_abc(1, 0, i) - dxR_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxR_abc(j, 1, i) = dxR_abc(j, 0, i);
      dxR_abc(j, 0, i) = ex(NUMCOLS - 1 - j, i);
    }
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    ex(i, 0) = abc_c0 * (ex(i, 2) + dxB_abc(0, 1, i)) +
               abc_c1 * (dxB_abc(0, 0, i) + dxB_abc(2, 0, i) - ex(i, 1) -
                         dxB_abc(1, 1, i)) +
               abc_c2 * dxB_abc(1, 0, i) - dxB_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxB_abc(j, 1, i) = dxB_abc(j, 0, i);
      dxB_abc(j, 0, i) = ex(i, j);
    }
  }

  for (int i = 0; i < NUMCOLS; i++)
  {
    ex(i, NUMROWS - 1) = abc_c0 * (ex(i, NUMROWS - 3) + dxT_abc(0, 1, i)) +
                         abc_c1 * (dxT_abc(0, 0, i) + dxT_abc(2, 0, i) - ex(i, NUMROWS - 2) -
                                   dxT_abc(1, 1, i)) +
                         abc_c2 * dxT_abc(1, 0, i) - dxT_abc(2, 1, i);

    for (int j = 0; j < 3; j++)
    {
      dxT_abc(j, 1, i) = dxT_abc(j, 0, i);
      dxT_abc(j, 0, i) = ex(i, NUMROWS - 1 - j);
    }
  }
}
*/

__global__ void compute_srcs(int cur_step, double *ex)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid != 0)
  {
    return;
  }

  // Sinusoidal source
  // 20 GHz
  pulse = 10 * sin(2 * pi * 7 * 1e8 * delx / 3e8 * cur_step);

  // pulse = 5 * exp(-.2 * (pow((t0 - cur_step) / spread, 2.0)));

  ex(10, NUMCOLS / 2) = pulse;
}

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step,
                        double *relative_eps, double *sigma, double abc_c0, double abc_c1, double abc_c2,
                        double *dxL_abc, double *dxR_abc, double *dxT_abc,
                        double *dxB_abc)
{
  updateHyGPU<<<NUMBLKS, NUMTHREADS>>>(ex, hy, hz);
  updateHzGPU<<<NUMBLKS, NUMTHREADS>>>(ex, hy, hz);
  updateEzGPU<<<NUMBLKS, NUMTHREADS>>>(ex, hy, hz, relative_eps, sigma);

  compute_srcs<<<NUMBLKS, NUMTHREADS>>>(cur_step, ex);

  apply_abc_bot_row<<<NUMBLKS, NUMTHREADS>>>(ex, abc_c0, abc_c1, abc_c2, dxT_abc, NUMROWS, NUMCOLS);
  apply_abc_top_row<<<NUMBLKS, NUMTHREADS>>>(ex, abc_c0, abc_c1, abc_c2, dxB_abc, NUMROWS, NUMCOLS);
  apply_abc_left_col<<<NUMBLKS, NUMTHREADS>>>(ex, abc_c0, abc_c1, abc_c2, dxL_abc, NUMROWS, NUMCOLS);
  apply_abc_right_col<<<NUMBLKS, NUMTHREADS>>>(ex, abc_c0, abc_c1, abc_c2, dxR_abc, NUMROWS, NUMCOLS);
}