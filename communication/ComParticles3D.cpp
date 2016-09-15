
#include <climits>
#include "ComParticles3D.h"

/** comunicate particles and receive particles to and from 6 processors */
void communicateParticles(long long buffer_size_x, long long buffer_size_y, long long buffer_size_z, double *b_Xleft, double *b_Xright, double *b_Yleft, double *b_Yright, double *b_Zleft, double *b_Zright, VirtualTopology3D * vct) {

  // DIR X
  if (buffer_size_x > INT_MAX) {
    cout << "ERROR: X-buffer too large: " << buffer_size_x << " > " << INT_MAX << endl;
  }
  communicateParticlesDIR((int) buffer_size_x, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Xright, b_Xleft);

  // DIR Y
  if (buffer_size_y > INT_MAX) {
    cout << "ERROR: Y-buffer too large: " << buffer_size_y << " > " << INT_MAX << endl;
  }
  communicateParticlesDIR((int) buffer_size_y, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Yright, b_Yleft);

  // DIR Z
  if (buffer_size_z > INT_MAX) {
    cout << "ERROR: Z-buffer too large: " << buffer_size_z << " > " << INT_MAX << endl;
  }
  communicateParticlesDIR((int) buffer_size_z, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), b_Zright, b_Zleft);

}
