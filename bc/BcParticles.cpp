#include <mpi.h>
#include "BcParticles.h"

/** boundary conditions for particles:
  <ul>
  <li>bcFace = 0 : loose particles</li>
  <li>bcFace = 1 : perfect mirror</li>
  <li>bcFace = 2 : riemission</li>
  </ul> */

/** apply left boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_left(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt, int bcFaceXleft) {
  switch (bcFaceXleft) {
    case 1:
      // perfect mirror
      *x = -*x;
      *u = -*u;
      break;
    case 2:
      // riemmission
      *x = -*x;
      double harvest, prob, theta;
      // u
      harvest = rand() / (double) RAND_MAX;
      prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
      harvest = rand() / (double) RAND_MAX;
      theta = 2.0 * M_PI * harvest;
      *u = fabs(ut * prob * cos(theta));
      // v
      *v = vt * prob * sin(theta);
      // w
      harvest = rand() / (double) RAND_MAX;
      prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
      harvest = rand() / (double) RAND_MAX;
      theta = 2.0 * M_PI * harvest;
      *w = wt * prob * cos(theta);
      break;
  }
}

/** apply right boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_right(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt, int bcFaceXright) {
  switch (bcFaceXright) {
    case 1:
      // perfect mirror
      *x = 2 * Lx - *x;
      *u = -*u;
      break;
    case 2:
      // riemmission
      double harvest, prob, theta;
      *x = 2 * Lx - *x;
      // u
      harvest = rand() / (double) RAND_MAX;
      prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
      harvest = rand() / (double) RAND_MAX;
      theta = 2.0 * M_PI * harvest;
      *u = -fabs(ut * prob * cos(theta));
      // v
      *v = vt * prob * sin(theta);
      // w
      harvest = rand() / (double) RAND_MAX;
      prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
      harvest = rand() / (double) RAND_MAX;
      theta = 2.0 * M_PI * harvest;
      *w = wt * prob * cos(theta);
      break;
  }
}

