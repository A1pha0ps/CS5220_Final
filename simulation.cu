#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <omp.h>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cuda.h>
using namespace std;

#include "common.h"



#define MIN_DEN 0.001
#define grad_w 0.0001
#define t_step_size 0.1

// define design and target space


int n_source;
int n_targets;
int n_design;
int n_mats;

int *target_x;
int *target_y;

int *design_x;
int *design_y;
double *design_density;
double *design_grad;


int *source_x;
int *source_y;
double *source_amp;
double *source_offset;
int * source_type;

int * mat_x;
int * mat_y;
double * mat_eps;
double * mat_sig;

double d_eps;
double d_sig;

// fundamental fields
double *dx, *ex, *hy, *hz;

// abcs
double *dxL_abc, *dxR_abc, *dxT_abc, *dxB_abc;

// describes properties of grid material
double *relative_eps, *sigma;
double damage_max;
int fdtd_steps, topop_steps, NUMROWS, NUMCOLS;

double mat_cost(){
	double cost = 0.0;
	for(int i = 0; i < n_design; ++i){
		cost += design_density[i];
	}
	cost/=n_design;	
	return cost * 0.1;
}
void init_mat_properties()
{
	// initialize free space
	for (int i = 0; i < NUMROWS; i++)
	{
		for (int j = 0; j < NUMCOLS; j++)
		{
			relative_eps(i, j) = 1;
			sigma(i, j) = 0;
		}
	}
	cout << "Free space done " << endl;
	//initialize materials other than design ps
	for (int i = 0; i < n_mats; ++i)
	{

		relative_eps(mat_x[i], mat_y[i]) = mat_eps[i];
		sigma(mat_x[i], mat_y[i]) = mat_sig[i];
	}
	cout << "Custom materials done " << endl;
}

void update_d_mat()
{
	//update material space with design params
	for (int i = 0; i < n_design; ++i)
	{
		
			cout << design_density[i] << " ";
		
		relative_eps(design_x[i], design_y[i]) = design_density[i] * (d_eps - 1) + 1;
		sigma(design_x[i], design_y[i]) = design_density[i] * d_sig;
	}
	cout << endl;
}

void set_d_mat(int dof, double dens)
{
	//set a single design material square to some density
	design_density[dof] = dens;
	relative_eps(design_x[dof], design_y[dof]) = dens * (d_eps - 1) + 1;
	sigma(design_x[dof], design_y[dof]) = dens * d_sig;
}


double simulate()
{
	reset_simulation(dx, ex, hy, hz);
	reset_bounds(dxL_abc, dxR_abc, dxT_abc, dxB_abc);

	double damage = 0.0;

	for (int step = 0; step < fdtd_steps; step++)
	{
		simulate_time_step(dx, ex, hy, hz, step, relative_eps, sigma, dxL_abc, dxR_abc, dxT_abc, dxB_abc, source_x, source_y, source_amp, source_offset, source_type, n_source);

		for (int i = 0; i < n_targets; ++i)
		{
			damage += pow(ex(target_x[i], target_y[i]),2)*delt*delx*delx;
		}
	}

	return damage;
}

void update_density()
{
	double sum = 0;
	for (int d = 0; d < n_design; ++d)
	{

		double old_val = design_density[d];
		// we need to estimate the gradient of the score - central difference with each design variable
		set_d_mat(d, old_val - grad_w);
		double back = simulate()/damage_max + mat_cost();
		set_d_mat(d, old_val + grad_w);
		double front = simulate()/damage_max + mat_cost();
		
		design_grad[d] = (front - back) / (2 * grad_w);
		sum += design_grad[d] * design_grad[d];
		//set it back
		set_d_mat(d, old_val);
	}

	sum = sqrt(sum);
	
	// take a step backward, penalize, and cap
	for (int d = 0; d < n_design; ++d)
	{
		design_density[d] -= t_step_size * design_grad[d]/sum;

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


vector<int> sx;
vector<int> sy;
vector<double> s_amp;
vector<int> s_type;
vector<double> s_off;

vector<int> mtx;
vector<int> mty;
vector<double> mEps;
vector<double> mSig;


vector<int> d_x;
vector<int> d_y;

vector<int> tx;
vector<int> ty;

void read_in_file(string f_name){

	ifstream file(f_name);
	string line;

	getline(file, line);

	//Read Field Size
	stringstream lineStream;
	lineStream.str(line);

	lineStream >> NUMROWS;
	lineStream >> NUMCOLS;
	cout << "field size read: " << NUMROWS << " x " << NUMCOLS << endl;
	//Read Time Steps, Top op steps
	getline(file, line);
	lineStream.clear();
	lineStream.str(line);
	lineStream >> fdtd_steps;
	lineStream >> topop_steps;
	cout << "step counts read: " << fdtd_steps << " fdtd, " << topop_steps << " topop" << endl;
	//Read in Source Points
	getline(file, line);
	lineStream.clear();
	lineStream.str(line);
	lineStream >> n_source;
	for(int i = 0; i < n_source; ++i){
		getline(file, line);
		lineStream.clear();
		lineStream.str(line);
		int x, y, type;
		double amp, off;
		lineStream >> x;
		lineStream >> y;
		lineStream >> amp;
		lineStream >> type;
		lineStream >> off;

		sx.push_back(x);
		sy.push_back(y);
		s_amp.push_back(amp);
		s_type.push_back(type);
		s_off.push_back(off);	
	}

	cout << "source info read " << n_source << " sources" << endl;
	//Read in Material Properties
	getline(file, line);
	lineStream.clear();
	lineStream.str(line);
	int mat_rects;
	lineStream >> mat_rects;

	for(int i = 0; i < mat_rects; ++i){
		getline(file, line);
		lineStream.clear();
		lineStream.str(line);
		int TL_X, TL_Y, L, W;
		double eps;
		double sig;
		lineStream >> TL_X;
		lineStream >> TL_Y;
		lineStream >> L;
		lineStream >> W;
		lineStream >> eps;
		lineStream >> sig;
		for(int x = 0; x < L; ++x){
			for(int y = 0; y < W; ++y){
				mtx.push_back(TL_X + x);
				mty.push_back(TL_Y + y);
				mEps.push_back(eps);
				mSig.push_back(sig);
			}
		}
	}
	n_mats = mtx.size();
	cout << "material info read " << n_mats << " mats"  << endl;
	//Read in targets
	getline(file, line);
	lineStream.clear();
	lineStream.str(line);
	int t_rects;
	lineStream >> t_rects;

	for(int i = 0; i < t_rects; ++i){
		getline(file, line);
		lineStream.clear();
		lineStream.str(line);
		int TLX, TLY, L, W;
		lineStream >> TLX;
		lineStream >> TLY;
		lineStream >> L;
		lineStream >> W;
		for(int x = 0; x < L; ++x){
			for(int y = 0; y < W; ++y){
				tx.push_back(TLX + x);
				ty.push_back(TLY + y);
			}
		}
	}
	n_targets = tx.size();
	cout << "target info read " << n_targets <<" targets" << endl;
	//Read in design parameter information
	getline(file, line);
	lineStream.clear();
	lineStream.str(line);
	int dp_rects;
	lineStream >> dp_rects;
	lineStream >> d_eps;
	lineStream >> d_sig;
	for(int i  = 0; i < dp_rects; ++i){
		getline(file, line);
		lineStream.clear();
		lineStream.str(line);
		int TLX, TLY, L, W;
		lineStream >> TLX;
		lineStream >> TLY;
		lineStream >> L;
		lineStream >> W;
		for(int x = 0; x < L; ++x){
			for(int y = 0; y < W; ++y){
				d_x.push_back(x + TLX);
				d_y.push_back(y + TLY);
			}
		}
	}

	n_design = d_x.size();
	cout << "design params done " << n_design << " params" << endl;
	//allocate necessary pointers

	dx = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
	ex = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
	hy = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
	hz = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

	relative_eps = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
	sigma = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

	dxL_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
	dxR_abc = (double *)malloc(NUMROWS * 6 * sizeof(double));
	dxT_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));
	dxB_abc = (double *)malloc(NUMCOLS * 6 * sizeof(double));


	design_grad = (double *)malloc(n_design * sizeof(double));
	design_density = (double *) malloc(n_design * sizeof(double));
	for(int i = 0; i < n_design; ++i){
		design_density[i] = 0.5;
	}

	target_x = tx.data();
	target_y = ty.data();

	design_x = d_x.data();
	design_y = d_y.data();




	source_x = sx.data();
	source_y = sy.data();
	source_amp = s_amp.data();
	source_offset = s_off.data();
	source_type = s_type.data();

	mat_x = mtx.data();
	mat_y = mty.data();
	mat_eps = mEps.data();
	mat_sig = mSig.data();



	cout << "storing done" << endl;
}

void do_top_op()
{

	cout << "starting topop: " << topop_steps << " total steps " << endl;
	init_mat_properties();
	cout << "initial materials set" << endl;
	update_d_mat();
	damage_max = simulate();
	for (int i = 0; i < topop_steps; ++i)
	{
		cout << "step " << i << "/" << topop_steps << " finished" << endl;
		update_density();
	}
}
void do_top_op_full_out(){
	ofstream damage_file("damage.out");
	ofstream score_file("score.out");
	ofstream cond_file("cond.out");
	ofstream perm_file("perm.out");

	cout << "starting topop: " << topop_steps << " total steps " << endl;
	init_mat_properties();
	cout << "initial materials set" << endl;

	update_d_mat();
	damage_max = simulate();
	damage_file << simulate()/damage_max << endl;
	score_file << simulate()/damage_max + mat_cost()<<endl;
	for(int i = 0; i < NUMROWS; ++i){
		for(int j = 0; j < NUMCOLS; ++j){
			cond_file << sigma(i, j) << " ";
			perm_file << relative_eps(i, j) << " ";
		}
		perm_file << endl;
		cond_file << endl;
	}
	for (int i = 0; i < topop_steps; ++i)
	{
		cout << "step " << i << "/" << topop_steps << " finished" << endl;
		update_density();
		damage_file << simulate()/damage_max << endl;
		score_file << simulate()/damage_max + mat_cost()<<endl;
		for(int i = 0; i < NUMROWS; ++i){
			for(int j = 0; j < NUMCOLS; ++ j){
				cond_file << sigma(i, j) << " ";
				perm_file << relative_eps(i, j) << " ";
			}
			perm_file << endl;
			cond_file << endl;
		}

	}
	damage_file.close();
	perm_file.close();
	cond_file.close();
	score_file.close();
}

int main(int argc, char *argv [])
{
	cout << "Start " << endl;
	int output = 0;
	string file_name;
	if(argc == 3){
		output = stoi(argv[1]);	
		file_name = argv[2];
	}
	else if(argc == 2){
		file_name = argv[1];
	}
	cout << "Got file name " << file_name << endl;
	read_in_file(file_name);
	cout << d_eps << " " << d_sig << endl;
	if(output == 0){
		do_top_op();
	}
	else{
		do_top_op_full_out();
	}


	return 0;
}
