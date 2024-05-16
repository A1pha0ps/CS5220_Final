#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

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

void reset_simulation(double *dx, double *ex, double *hy, double *hz)
{

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
}
void reset_bounds(double * dxL_abc, double * dxR_abc, double * dxT_abc, double * dxB_abc){

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

void apply_abc(double *ex, double *dxL_abc, double *dxR_abc, double *dxT_abc, double *dxB_abc)
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

void simulate_time_step(double *dx, double *ex, double *hy, double *hz, int cur_step,
		double *relative_eps, double *sigma, double *dxL_abc, double *dxR_abc, double *dxT_abc,
		double *dxB_abc,
		int * s_x, int * s_y, double * s_amp, double * s_off, int * s_type, int s_count)
{
	updateHFields(dx, ex, hy, hz, cur_step);

	// updateTfSf(dx, ex, hy, hz, cur_step);

	updateEFields(dx, ex, hy, hz, cur_step, relative_eps, sigma);

	//Apply loads

	for (int i = 0; i < s_count; ++i){
		double pulse;
		if(s_type[i] == 0){
			//Sin pulse
			pulse = s_amp[i] * sin(2 * pi * 7 * 1e8 * delx/3e8 * (cur_step*delt - s_off[i]));
		}
		else if(s_type[i] == 1){
			//Gaussian pulse
			pulse = s_amp[i] * exp(-.5 * (pow((t0 - (cur_step*delt - s_off[i])) / spread, 2.0)));
		}
		ex(s_x[i], s_y[i]) = pulse;
	}


	apply_abc(ex, dxL_abc, dxR_abc, dxT_abc, dxB_abc);

	// ex(NUMROWS / 2 - 30, NUMCOLS / 2) = pulse = 5 * exp(-.5 * (pow((t0 - cur_step) / spread, 2.0)));

	// if (cur_step == nsteps - 1)
	// {
	//   printf("D-field time: %f seconds\n", dtime);
	//   printf("E-field time: %f seconds\n", etime);
	//   printf("Hy-field time: %f seconds\n", hytime);
	//   printf("Hz-field time: %f seconds\n", hztime);
	// }
}
