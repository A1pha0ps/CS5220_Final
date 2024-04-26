#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <omp.h>

#include <iostream>
using namespace std;

#include "common.h"

int main()
{

  double *dx, *ex, *hy, *hz;

  dx = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  ex = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  hy = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));
  hz = (double *)malloc(NUMROWS * NUMCOLS * sizeof(double));

  init_simulation(dx, ex, hy, hz);

  for (int step = 0; step < nsteps; step++)
  {
    // cout << step << endl;
    simulate_time_step(dx, ex, hy, hz, step);
  }

  return 0;
}
