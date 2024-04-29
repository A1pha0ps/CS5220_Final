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

  FILE *fp;

  fp = fopen("ex_field.txt", "w");

  for (int step = 1; step <= nsteps; step++)
  {
    // cout << step << endl;
    simulate_time_step(dx, ex, hy, hz, step);

    /* Write the E field out to a file "Ez" */

    // fprintf(fp, "%d\n", step);

    for (int i = 0; i < NUMROWS; i++)
    {
      for (int j = 0; j < NUMCOLS; j++)
      {
        fprintf(fp, "%6.3f ", ex[i * NUMROWS + j]);
      }
      fprintf(fp, " \n");
    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  return 0;
}
