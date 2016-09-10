#include <mpi.h>
#include "BcParticles.h"

/** apply left mirror boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_left_mirror(double *x, double *u, double Lx) {
  // perfect mirror
  *x = -*x;
  *u = -*u;
}

/** apply right mirror boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_right_mirror(double *x, double *u, double Lx) {
  // perfect mirror
  *x = 2 * Lx - *x;
  *u = -*u;
}

/** apply left riemission boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_left_riemission(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt) {
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
}

/** apply right riemission boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_right_riemission(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt) {
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
}

