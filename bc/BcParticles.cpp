#include <mpi.h>
#include "BcParticles.h"

/** available boundary conditions:
  <ul>
  <li>bcFace = 0 : loose particles</li>
  <li>bcFace = 1 : perfect mirror</li>
  <li>bcFace = 2 : re-emission</li>
  </ul> */

/** apply left degenerated boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_left_degenerated(double *x, double Lx) {
  // degenerated case
  *x = *x + Lx;
}

/** apply right degenerated boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_right_degenerated(double *x, double Lx) {
  // degenerated case
  *x = *x - Lx;
}

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

/** apply left re-emission boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_left_reemission(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt) {
  // re-emission
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

/** apply right re-emission boundary condition for particles along one coordinate direction (here named 'x') */
void BCpart_right_reemission(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt) {
  // re-emission
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

