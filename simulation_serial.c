#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define KE 200

/* FD1D_1.3.c.
/* Simulation of a pulse hitting a dielectric medium */

main()
{
  float ex[KE], hy[KE];
  float cb[KE];
  int n, k, kc, ke, kstart, nsteps;
  float ddx, dt, T, epsz, epsilon, sigma, eaf;
  float t0, spread, pi, pulse;
  FILE *fp, *fopen();
  float ex_low_m1, ex_low_m2, ex_high_m1, ex_high_m2;

  for (k = 1; k <= KE; k++)
  { /* Initialize to free space */
    cb[k] = .5;
  }

  printf("Dielectric starts at --> ");
  scanf("%d", &kstart);
  printf("Epsilon --> ");
  scanf("%f", &epsilon);
  printf("%d %6.2f \n", kstart, epsilon);

  for (k = kstart; k <= KE; k++)
  {
    cb[k] = .5 / epsilon;
  }

  for (k = 1; k <= KE; k++)
  {
    printf("%2d %4.2f\n", k, cb[k]);
  } 
  
  /* Main part of the program */

  while (nsteps > 0)
  {
    printf("nsteps --> ");
    scanf("%d", &nsteps);
    printf("%d \n", nsteps);
    for (n = 1; n <= nsteps; n++)
    {
      T = T + 1;
      /* Calculate the Ex field */
      for (k = 0; k < KE; k++)
      {
        ex[k] = ex[k] + cb[k] * (hy[k - 1] - hy[k]);
      }
      /* Put a Gaussian pulse at the low end */
      pulse = exp(-.5 * (pow((t0 - T) / spread, 2.0)));
      ex[5] = ex[5] + pulse;
      printf("%5.1f %6.2f %6.2f\n", T, pulse, ex[5]);
    }
  }
}
