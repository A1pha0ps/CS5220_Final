#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <omp.h>
#include <random>
#include <iostream>
#include <cuda.h>
using namespace std;

#include "common.h"
#define NUM_SOURCE NUMCOLS / 2
#define NUM_TARGET (NUMROWS / 8) * (NUMCOLS / 8) * 2
#define NUM_DESIGN (NUMROWS * 3 / 4) * (NUMCOLS / 4)
#define MIN_DEN 0.001
#define grad_w 0.0005
#define OP_ITERS 20
double abc_c0, abc_c1, abc_c2;

// define design and target space

int *target_x;
int *target_y;

int *design_x;
int *design_y;
double *design_density;
double *design_grad;
double d_max;

int *source_x;
int *source_y;
double *source_amp;
double *source_offset;
double *source_freq;

// Silicon for Device
double t_perm = 19.3;
double t_cond = 5.21;
// Copper for Shield
double d_perm = 1.0;
double d_cond = 5.8e7;

// fundamental fields
double *dx, *ex, *hy, *hz;

// abcs
double *dxL_abc, *dxR_abc, *dxT_abc, *dxB_abc;

// describes properties of grid material
double *relative_eps, *sigma;

int nsteps, NUMROWS, NUMCOLS;

void init_abc()
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
}

void init_mat_properties()
{
  /*
     Code to specify properties of cell

    Reference:
    Material: Lossy Dielectric (Silicon @ 20 GHz)
    Relative Permitivity: 19.3
    Conductivity: 5.21 S/m
   */
  // initialize free space
  for (int i = 0; i < NUMROWS; i++)
  {
    for (int j = 0; j < NUMCOLS; j++)
    {
      relative_eps(i, j) = 1;
      sigma(i, j) = 0;
    }
  }

  // initialize targets
  // for (int i = 0; i < NUM_TARGET; ++i)
  // {
  //   relative_eps(target_x[i], target_y[i]) = t_perm;
  //   sigma(target_x[i], target_y[i]) = t_cond;
  // }
}

void update_d_mat()
{
  for (int i = 0; i < NUM_DESIGN; ++i)
  {
    relative_eps(design_x[i], design_y[i]) = design_density[i] * (d_perm - 1) + 1;
    sigma(design_x[i], design_y[i]) = design_density[i] * d_cond;
  }
}

void set_d_mat(int dof, double dens)
{
  relative_eps(design_x[dof], design_y[dof]) = dens * (d_perm - 1) + 1;
  sigma(design_x[dof], design_y[dof]) = dens * d_cond;
}

void free_all()
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

void all_init()
{
  init_simulation(dx, ex, hy, hz);

  init_abc();

  init_mat_properties();
}

double simulate()
{
  init_simulation(dx, ex, hy, hz);
  init_abc();

  double damage = 0.0;

  for (int step = 0; step < nsteps; step++)
  {
    // cout << step << endl;
    simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma, abc_c0, abc_c1, abc_c2, dxL_abc, dxR_abc, dxT_abc, dxB_abc);

    for (int i = 0; i < NUM_TARGET; ++i)
    {
      damage += abs(ex(target_x[i], target_y[i])) / (nsteps * NUMCOLS * NUMROWS);
    }
  }

  return damage;
}

void update_density()
{
  double sum = 0.0;
  for (int d = 0; d < NUM_DESIGN; ++d)
  {
    cout << "Doing " << d << "/" << NUM_DESIGN << endl;
    double old_val = design_density[d];
    // we need to estimate the gradient of the damage - central difference with each design variable
    set_d_mat(d, old_val - grad_w);
    double back = simulate();
    set_d_mat(d, old_val + grad_w);
    double front = simulate();
    design_grad[d] = (front - back) / (2 * grad_w);
    sum += design_grad[d];
    set_d_mat(d, old_val);
  }

  // take a step backward
  for (int d = 0; d < NUM_DESIGN; ++d)
  {
    design_density[d] += design_grad[d] / sum;
    if (design_density[d] < MIN_DEN)
    {
      design_density[d] = MIN_DEN;
    }
    if (design_density[d] >= 1.0)
    {
      design_density[d] = 1.0;
    }
  }
  update_d_mat();
}

void simulate_with_output(FILE *fp)
{

  init_simulation(dx, ex, hy, hz);
  init_abc();

  for (int step = 0; step < nsteps; step++)
  {
    // cout << step << endl;
    simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma, abc_c0, abc_c1, abc_c2, dxL_abc, dxR_abc, dxT_abc, dxB_abc);
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
}

void do_top_op()
{
  design_grad = (double *)malloc(NUM_DESIGN * sizeof(double));
  // set up a dummy problem: 1 source, boxlike design space and target

  vector<int> sx;
  vector<int> sy;
  vector<double> samp;
  vector<double> off;
  vector<double> freq;

  for (int i = -NUMCOLS / 4; i < NUMCOLS / 4; ++i)
  {
    sx.push_back(NUMCOLS * 1 / 8);
    sy.push_back(NUMROWS / 2 + i);
    samp.push_back(320000);
    off.push_back(0);
    freq.push_back(7);
  }
  source_x = sx.data();
  source_y = sy.data();
  source_amp = samp.data();
  source_offset = off.data();
  source_freq = freq.data();

  int t_ind = 0;
  int d_ind = 0;
  vector<int> tx;
  vector<int> ty;
  vector<int> ddx;
  vector<int> ddy;
  vector<double> ddden;
  for (int i = NUMROWS / 8; i < NUMROWS / 4; ++i)
  {
    for (int j = NUMCOLS * 3 / 4; j < 7 * NUMCOLS / 8; ++j)
    {
      tx.push_back(j);
      ty.push_back(i);
    }
  }

  for (int i = 3 * NUMROWS / 4; i < 7 * NUMROWS / 8; ++i)
  {
    for (int j = NUMCOLS * 3 / 4; j < 7 * NUMCOLS / 8; ++j)
    {
      tx.push_back(j);
      ty.push_back(i);
    }
  }

  for (int i = NUMROWS / 8; i < NUMROWS * 7 / 8; ++i)
  {
    for (int j = 3 * NUMCOLS / 8; j < 5 * NUMCOLS / 8; ++j)
    {
      ddx.push_back(j);
      ddy.push_back(i);
      ddden.push_back(0.1);
    }
  }

  //	for(int i = 0; i < OP_ITERS; ++i){
  //		update_density();
  //	}

  target_x = tx.data();
  target_y = ty.data();

  design_x = ddx.data();
  design_y = ddy.data();
  design_density = ddden.data();
  init_mat_properties();
  FILE *DATA;
  DATA = fopen("ML_DATA.txt", "w");
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(0.1, 1.0);

  for (int i = 0; i < 6000; ++i)
  {
    cout << i << endl;
    // populate design density with random stuff
    for (int d = 0; d < NUM_DESIGN; ++d)
    {
      design_density[d] = dis(gen);
      fprintf(DATA, "%.10f,", design_density[d]);
    }
    update_d_mat();
    fprintf(DATA, "%.10f\n", simulate());
  }
  fclose(DATA);

  FILE *F_DAMAGE;
  F_DAMAGE = fopen("damage_plot", "w");
  FILE *FPB;
  FPB = fopen("e_field_before.txt", "w");
  FILE *FPA;
  FPA = fopen("e_field_after.txt", "w");
  update_d_mat();
  simulate_with_output(FPB);
  fprintf(F_DAMAGE, "%6.3f\n", simulate());
  FILE *fp;
  fp = fopen("conductance.txt", "w");
  for (int i = 0; i < 30; ++i)
  {
    update_density();
    fprintf(F_DAMAGE, "%6.3f\n", simulate());

    for (int i = 0; i < NUMCOLS; ++i)
    {
      for (int j = 0; j < NUMROWS; ++j)
      {
        if (i < 3 * NUMCOLS / 4 && i >= NUMCOLS / 2)
        {
          fprintf(fp, "%6.3f ", sigma[i * NUMROWS + j] / d_cond);
        }
        else
        {
          fprintf(fp, "%6.3f ", 0.0);
        }
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  fclose(F_DAMAGE);
  simulate_with_output(FPA);
  free(design_grad);
}

int main(int argc, char **argv)
{

  if (argc != 3)
  {
    cout << "Wrong number of arguments provided: ./serial <NUM_TIME_STEPS> <SIDE_LENGTH_OF_GRID>" << endl;
    return 0;
  }

  nsteps = stoi(argv[1]);
  NUMROWS = stoi(argv[2]);
  NUMCOLS = stoi(argv[2]);

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

  // do_top_op();

  all_init();

  FILE *fp;

  fp = fopen("ex_field.txt", "w");

  simulate_with_output(fp);
  free_all();

  return 0;
}