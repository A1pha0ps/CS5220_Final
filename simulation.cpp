#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <mpi.h>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#include "common.h"


#define OMP_THREADS 64
#define MIN_DEN 0.001
#define grad_w 0.0001
#define t_step_size 0.1

// define design and target space

int r_id;
int num_procs;
int d_perr_id;

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
double *lcl_grad;

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

	//initialize materials other than design ps
	for (int i = 0; i < n_mats; ++i)
	{

		relative_eps(mat_x[i], mat_y[i]) = mat_eps[i];
		sigma(mat_x[i], mat_y[i]) = mat_sig[i];
	}

}

void update_d_mat()
{
	//update material space with design params
	for (int i = 0; i < n_design; ++i)
	{	
		relative_eps(design_x[i], design_y[i]) = design_density[i] * (d_eps - 1) + 1;
		sigma(design_x[i], design_y[i]) = design_density[i] * d_sig;
	}

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
	for (int d = 0; d < d_perr_id; ++d)
	{
		int d_id = r_id * d_perr_id + d;
		if(d_id < n_design){
			double old_val = design_density[d_id];
			// we need to estimate the gradient of the score - central difference with each design variable
			set_d_mat(d_id, old_val - grad_w);
			double back = simulate()/damage_max + mat_cost();
			set_d_mat(d_id, old_val + grad_w);
			double front = simulate()/damage_max + mat_cost();

			lcl_grad[d] = (front - back) / (2 * grad_w);
			sum += lcl_grad[d] * lcl_grad[d];
			//set it back
			set_d_mat(d_id, old_val);
		}


	}
	
	//combine magnitudes and gradients
	double total_mag;
	MPI_Allgather(lcl_grad, d_perr_id, MPI_DOUBLE, design_grad, d_perr_id, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
cout << "rank " << r_id << " updating density" << endl;
	MPI_Allreduce(&sum, &total_mag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	total_mag = sqrt(total_mag);

	// take a step backward, penalize, and cap
	for (int d = 0; d < n_design; ++d)
	{
		design_density[d] -= t_step_size * design_grad[d]/total_mag;

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
	
	if(r_id == 0){
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
	}
	int * counts = new int[8];
	if(r_id == 0){
		counts[0] = NUMROWS;
		counts[1] = NUMCOLS;
		counts[2] = n_source;
		counts[3] = n_mats;
		counts[4] = n_targets;
		counts[5] = n_design;
		counts[6] = fdtd_steps;
		counts[7] = topop_steps;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(counts, 8, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	NUMROWS = counts[0];
	NUMCOLS = counts[1];
	n_source = counts[2];
	n_mats = counts[3];
	n_targets = counts[4];
	n_design = counts[5];
	fdtd_steps = counts[6];
	topop_steps = counts[7];


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

	d_perr_id = (n_design - 1)/num_procs + 1;
	lcl_grad = (double *)malloc(d_perr_id*sizeof(double));	
	design_grad = (double *)malloc(num_procs * d_perr_id * sizeof(double));
	design_density = (double *) malloc(n_design * sizeof(double));

	for(int i = 0; i < n_design; ++i){
		design_density[i] = 0.5;
	}
	for(int i = 0; i < num_procs * d_perr_id; ++i){
		design_grad[i] = 0.0;
	}
	for(int i = 0; i < d_perr_id; ++i){
		lcl_grad[i] = 0.0;
	}
	if(r_id == 0){
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
	}
	else{
		//malloc pointers that are not on r_id 0	
		target_x = (int *)malloc(n_targets * sizeof(int));
		target_y = (int *)malloc(n_targets * sizeof(int));

		design_x = (int *)malloc(n_design * sizeof(int));
		design_y = (int *)malloc(n_design * sizeof(int));

		source_x = (int *)malloc(n_source * sizeof(int));
		source_y = (int *)malloc(n_source * sizeof(int));
		source_amp = (double *)malloc(n_source * sizeof(double));
		source_offset = (double *) malloc(n_source * sizeof(double));
		source_type = (int *) malloc(n_source * sizeof(int));

		mat_x = (int *) malloc(n_mats * sizeof(int));
		mat_y = (int *) malloc(n_mats * sizeof(int));
		mat_eps = (double *) malloc(n_mats * sizeof(double));
		mat_sig = (double *) malloc(n_mats * sizeof(double));


	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(target_x, n_targets, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(target_y, n_targets, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(design_x, n_targets, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(design_y, n_targets, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(source_x, n_source, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(source_y, n_source, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(source_amp, n_source, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(source_offset, n_source, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(source_type, n_source, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(mat_x, n_mats, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(mat_y, n_mats, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(mat_eps, n_mats, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(mat_sig, n_mats, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(r_id == 0){
		cout << "storing done" << endl;
	}
}

void do_top_op()
{


	init_mat_properties();

	update_d_mat();
	damage_max = simulate();
	for (int i = 0; i < topop_steps; ++i)
	{	
		if(r_id == 0){
			cout << "doing " << i << " / " << topop_steps << endl;
		}
		update_density();
	}
}
void do_top_op_full_out(){


	ofstream cond_file("cond.out");


	cout << "starting topop: " << topop_steps << " total steps " << endl;
	init_mat_properties();
	cout << "initial materials set" << endl;

	update_d_mat();



	for(int i = 0; i < NUMROWS; ++i){
		for(int j = 0; j < NUMCOLS; ++j){
			cond_file << sigma(i, j) << " ";

		}

		cond_file << endl;
	}

	cond_file << endl;

	for (int i = 0; i < topop_steps; ++i)
	{
		cout << "step " << i << "/" << topop_steps << " finished" << endl;
		update_density();


		for(int i = 0; i < NUMROWS; ++i){
			for(int j = 0; j < NUMCOLS; ++ j){
				cond_file << sigma(i, j) << " ";

			}

			cond_file << endl;
		}
		cond_file << endl;


	}


	cond_file.close();

}

int main(int argc, char *argv [])
{	

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &r_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if(r_id == 0){
		cout << "Start " << endl;
	}
	int output = 0;
	string file_name;
	
	if(argc == 3){
		output = stoi(argv[1]);	
		file_name = argv[2];
	}
	else if(argc == 2){
		file_name = argv[1];
	}

	read_in_file(file_name);
	if(r_id == 0){
		cout << d_eps << " " << d_sig << endl;
	}
	if(output == 0){
		do_top_op();
	}
	else if(num_procs == 1){
		do_top_op_full_out();
	}
	else{
		if(r_id == 0){
			cout << "bad input" << endl;
		}
	}

	MPI_Finalize();
	return 0;
}
