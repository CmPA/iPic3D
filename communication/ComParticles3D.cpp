
#include "ComParticles3D.h"
#include "ipicdefs.h"

/** comunicate particles and receive particles to and from 6 processors */
void communicateParticles(int buffer_size, double *b_Xleft, double *b_Xright, double *b_Yleft, double *b_Yright, double *b_Zleft, double *b_Zright, VirtualTopology3D * vct) {
  // DIR X
  communicateParticlesDIR(buffer_size, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Xright, b_Xleft);
  // DIR Y
  communicateParticlesDIR(buffer_size, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Yright, b_Yleft);
  // DIR Z
  communicateParticlesDIR(buffer_size, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Zright, b_Zleft);
}

/** communicate the number of particles are not in the right domain*/
int reduceNumberParticles(int rightDomain) {
  int result = 0;
  former_MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&rightDomain, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}

/** communicate the maximum number of particles from a domain */
int reduceMaxNpExiting(int npExitingMax) {
  int result = 0;
  former_MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&npExitingMax, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return (result);
}

