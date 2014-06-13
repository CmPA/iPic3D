/***************************************************************************
  BcFields3D.h  -  Library to manage boundary conditions for fields 
  -------------------
begin                : Fri Jan 2009
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef BcFields_H
#define BcFields_H

#include "VirtualTopology3D.h"

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, double ***vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);

/** set the boundary condition on boundaries */
void BCface_P(int nx, int ny, int nz, double ***vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);

/** set the boundary condition on boundaries Particles*/
void BCface_P(int nx, int ny, int nz, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);

#endif
