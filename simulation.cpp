#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <omp.h>

#include <iostream>
using namespace std;

#include "common.h"
#define NUM_SOURCE  1
#define NUM_TARGET  (NUMROWS/4)*(NUMCOLS/8)
#define NUM_DESIGN  (NUMROWS)*(NUMCOLS/4)
#define MIN_DEN 0.001
#define grad_w 0.0005
#define OP_ITERS 20
double abc_c0, abc_c1, abc_c2;

//define design and target space

int * target_x;
int * target_y;

int * design_x;
int * design_y;
double * design_density;
double * design_grad;
double d_max;

int * source_x;
int * source_y;
double * source_amp;
double * source_offset;
double * source_freq;

//Silicon for Device
double t_perm = 19.3;
double t_cond = 5.21;
//Copper for Shield
double d_perm = 1;
double d_cond = 5.8e7;



// fundamental fields
double *dx, *ex, *hy, *hz;

// abcs
double *dxL_abc, *dxR_abc, *dxT_abc, *dxB_abc;

// describes properties of grid material
double *relative_eps, *sigma;


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
//initialize free space
	for (int i = 0; i < NUMROWS; i++)
	{
		for (int j = 0; j < NUMCOLS; j++)
		{
			relative_eps(i, j) = 1;
			sigma(i, j) = 0;
		}
	}

//initialize targets
	for(int i = 0; i < NUM_TARGET; ++i){
		relative_eps(target_x[i], target_y[i]) = t_perm;
		sigma(target_x[i], target_y[i]) = t_cond;
	}

}

void update_d_mat(){
	for(int i = 0; i < NUM_DESIGN; ++i){
		relative_eps(design_x[i], design_y[i]) = design_density[i] * (d_perm - 1) + 1;
		sigma(design_x[i], design_y[i]) = design_density[i] * d_cond; 
	}

}
void set_d_mat(int dof, double dens){
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
void all_init(){

	init_simulation(dx, ex, hy, hz);
	init_abc();
	init_mat_properties();

}
double simulate(){
	init_simulation(dx, ex, hy, hz);
	init_abc();	

	double damage = 0.0;

	for (int step = 0; step < nsteps; step++)
	{
		// cout << step << endl;
		simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma);

		apply_source(ex, NUM_SOURCE, step, source_x, source_y, source_amp, source_offset, source_freq);

		apply_abc(ex, abc_c0, abc_c1, abc_c2, dxL_abc, dxR_abc, dxT_abc, dxB_abc);

		for(int i = 0; i < NUM_TARGET; ++i){
			damage += abs(ex(target_x[i], target_y[i]))/(nsteps * NUMCOLS * NUMROWS);
		}
	}
	
	return damage;



}


void update_density(){
	double sum = 0.0;
	for(int d = 0; d < NUM_DESIGN; ++d){
		cout << "Doing " << d << endl;
		double old_val = design_density[d];
		//we need to estimate the gradient of the damage - central difference with each design variable
		set_d_mat(d, old_val - grad_w);
		double back = simulate();	
		set_d_mat(d, old_val + grad_w);
		double front = simulate();
		design_grad[d] = (front - back)/(2*grad_w);
		cout << design_grad[d] << endl;
		sum += design_grad[d];
		set_d_mat(d, old_val);
	}

	


	//take a step backward
	for(int d = 0; d < NUM_DESIGN; ++d){
		design_density[d] -= design_grad[d];
		if(design_density[d] < MIN_DEN){
			design_density[d] = MIN_DEN;
		}
		if(design_density[d] >= 1.0){
			design_density[d] = 1.0;
		}
	}
	update_d_mat();
}

void simulate_with_output(){
	init_mat_properties();
	init_simulation(dx, ex, hy, hz);
	init_abc();	


	FILE *fp;

	fp = fopen("ex_field.txt", "w");

	for (int step = 0; step < nsteps; step++)
	{
		// cout << step << endl;
		simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma);

		apply_source(ex, NUM_SOURCE, step, source_x, source_y, source_amp, source_offset, source_freq);

		apply_abc(ex, abc_c0, abc_c1, abc_c2, dxL_abc, dxR_abc, dxT_abc, dxB_abc);


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
void do_top_op(){
	design_grad = (double *)malloc(NUM_DESIGN *sizeof(double));


	source_x = (int *)malloc(NUM_SOURCE *sizeof(int));
	source_y = (int *)malloc(NUM_SOURCE *sizeof(int));
	source_amp = (double *)malloc(NUM_SOURCE *sizeof(double));
	source_offset = (double *)malloc(NUM_SOURCE *sizeof(double));
	source_freq = (double *)malloc(NUM_SOURCE *sizeof(double));

	//set up a dummy problem: 1 source, boxlike design space and target
	source_x[0] = NUMROWS/8;
	source_y[0] = NUMCOLS/8;
	source_amp[0] = 16000.0;
	source_offset[0] = 0.0;
	source_freq[0] = 7.0;
	int t_ind = 0;
	int d_ind = 0;
	vector<int> tx;
	vector<int> ty;
	vector<int> ddx;
	vector<int> ddy;
	vector<double> ddden;	
	for(int i = 0; i < NUMROWS; ++i){
		for(int j = NUMCOLS/2; j < NUMCOLS; ++j){
			if(j >= NUMCOLS * 3/4 && j < NUMCOLS*7/8 && i >= 3*NUMROWS/8 && i < 5*NUMROWS/8){

				tx.push_back(j);
				ty.push_back(i);	

			}
			else if(j < 3*NUMCOLS/4){	
				ddx.push_back(j);
				ddy.push_back(i);
				ddden.push_back(0.5);
			}
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
	//update_density();
	simulate_with_output();
	free(design_grad);
	free(source_x);
	free(source_y);
	free(source_amp);
	free(source_offset);
	free(source_freq);
}
int main()
{

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
	//all_init();
	do_top_op();

	//simulate_with_output();
	free_all();

	return 0;
}
