/*******************************************************************************************************
  Grid3DCUCU.h  -  uniform cartesian 3D local grid for each process, including che guard cells
  -------------------

 *******************************************************************************************************/

#ifndef GRID3DCU_H
#define GRID3DCU_H

#include <iostream>

#include "Grid.h"
#include "CollectiveIO.h"
#include "ComInterpNodes3D.h"
#include "ComNodes3D.h"
#include "VirtualTopology3D.h"
#include "Alloc.h"

using std::cout;
using std::endl;

/**
 * Uniform cartesian local grid 3D
 *
 * @date Fri Jun 4 2008
 * @par Copyright:
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 3.0
 *
 */
class Grid3DCU                  // :public Grid
{
public:
  /** constructor */
  Grid3DCU(CollectiveIO * col, VirtualTopology3D * vct);
  /** destructor */
  ~Grid3DCU();
  /** allocate grid arrays for this domain */
  void allocate(CollectiveIO * ptC, VirtualTopology3D * ptVCT);
  /** deallocate grid arrays for this domain */
  void deallocate();
  /** print grid info */
  void print(VirtualTopology3D * ptVCT);
  /** calculate a derivative along a direction on nodes */
  void derivN(array_ref3_double& derN,
    const_arr4_double& scFieldC, int ns, int dir);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradC2N(array_ref3_double& gradXN, array_ref3_double& gradYN, array_ref3_double& gradZN,
    const_arr3_double& scFieldC);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradN2C(array_ref3_double& gradXC, array_ref3_double& gradYC, array_ref3_double& gradZC,
    const_arr3_double& scFieldN);
  /** calculate divergence on central points, given a vector field defined on nodes  */
  void divN2C(array_ref3_double& divC,
    const_arr3_double& vecFieldXN,
    const_arr3_double& vecFieldYN,
    const_arr3_double& vecFieldZN);
  /** calculate divergence on nodes, given a vector field defined on central points  */
  void divC2N(array_ref3_double& divN,
    const_arr3_double& vecFieldXC,
    const_arr3_double& vecFieldYC,
    const_arr3_double& vecFieldZC);
  /** calculate curl on nodes, given a vector field defined on central points  */
  void curlC2N(array_ref3_double& curlXN, array_ref3_double& curlYN,
    array_ref3_double& curlZN,
    const_arr3_double& vecFieldXC,
    const_arr3_double& vecFieldYC,
    const_arr3_double& vecFieldZC);
  /** calculate curl on central points, given a vector field defined on nodes  */
  void curlN2C(array_ref3_double& curlXC, array_ref3_double& curlYC, array_ref3_double& curlZC,
    const_arr3_double& vecFieldXN,
    const_arr3_double& vecFieldYN,
    const_arr3_double& vecFieldZN);

  /** calculate divergence on central points, given a Tensor field defined on nodes  */
  void divSymmTensorN2C(array_ref3_double& divCX, array_ref3_double& divCY, array_ref3_double& divCZ,
    const_arr4_double& pXX,
    const_arr4_double& pXY,
    const_arr4_double& pXZ,
    const_arr4_double& pYY,
    const_arr4_double& pYZ,
    const_arr4_double& pZZ, int ns);

  /** calculate laplacian on nodes, given a scalar field defined on nodes */
  void lapN2N(array_ref3_double& lapN,
    const_arr3_double& scFieldN, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
  void lapC2Cpoisson(array_ref3_double& lapC,
    array_ref3_double& scFieldC, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points */
  void lapC2C(array_ref3_double& lapC,
    const_arr3_double& scFieldC, VirtualTopology3D * vct);

  /** calculate divergence on boundaries */
  void divBCleft(array_ref3_double& divBC,
    const_arr3_double& vectorX,
    const_arr3_double& vectorY,
    const_arr3_double& vectorZ, int leftActiveNode, int dirDER);
  /** calculate divergence on boundaries */
  void divBCright(array_ref3_double& divBC,
    const_arr3_double& vectorX,
    const_arr3_double& vectorY,
    const_arr3_double& vectorZ, int rightActiveNode, int dirDER);
  /** calculate derivative on boundaries */
  void derBC(array_ref3_double& derBC,
    const_arr3_double& vector, int leftActiveNode, int dirDER);


  /** interpolate on nodes from central points */
  void interpC2N(array_ref3_double& vecFieldN, const_arr3_double& vecFieldC);
  /** interpolate on central points from nodes */
  void interpN2C(array_ref3_double& vecFieldC, const_arr3_double& vecFieldN);
  /** interpolate on central points from nodes */
  void interpN2C(array_ref4_double& vecFieldC, int ns, const_arr4_double& vecFieldN);

  // /////////// PRIVATE VARIABLES //////////////
private:
  /** number of cells - X direction, including + 2 (guard cells) */
  int nxc;
  /** number of nodes - X direction, including + 2 extra nodes for guard cells */
  int nxn;
  /** number of cell - Y direction, including + 2 (guard cells) */
  int nyc;
  /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
  int nyn;
  /** number of cell - Z direction, including + 2 (guard cells) */
  int nzc;
  /** number of nodes - Z direction, including + 2 extra nodes for guard cells */
  int nzn;
  /** dx = space step - X direction */
  double dx;
  /** dy = space step - Y direction */
  double dy;
  /** dz = space step - Z direction */
  double dz;
  /** invdx = 1/dx */
  double invdx;
  /** invdy = 1/dy */
  double invdy;
  /** invdz = 1/dz */
  double invdz;
  /** invol = inverse of volume*/
  double invVOL;
  /** node coordinate */
  double *node_xcoord;
  double *node_ycoord;
  double *node_zcoord;
  /** center coordinate */
  double *center_xcoord;
  double *center_ycoord;
  double *center_zcoord;
  /** local grid boundaries coordinate  */
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;

public: // accessors (inline)
  int getNXC() { return (nxc); }
  int getNXN() { return (nxn); }
  int getNYC() { return (nyc); }
  int getNYN() { return (nyn); }
  int getNZC() { return (nzc); }
  int getNZN() { return (nzn); }
  double getDX() { return (dx); }
  double getDY() { return (dy); }
  double getDZ() { return (dz); }
  //
  // coordinate accessors
  //
  // calculated equivalents (preferred for accelerator?):
  //const double &calcXN(int X) { return xStart+(X-1)*dx;}
  //const double &calcYN(int Y) { return yStart+(Y-1)*dy;}
  //const double &calcZN(int Z) { return zStart+(Z-1)*dz;}
  const double &getXN(int X) { return node_xcoord[X];}
  const double &getYN(int Y) { return node_ycoord[Y];}
  const double &getZN(int Z) { return node_zcoord[Z];}
  const double &getXC(int X) { return center_xcoord[X];}
  const double &getYC(int Y) { return center_ycoord[Y];}
  const double &getZC(int Z) { return center_zcoord[Z];}
  //
  // The following could be eliminated in favor of the previous
  // unless we truly anticipate generalizing to a deformed
  // logically cartesian mesh.  See issue #40.
  //
  const double &getXN(int X, int Y, int Z) { return node_xcoord[X];}
  const double &getYN(int X, int Y, int Z) { return node_ycoord[Y];}
  const double &getZN(int X, int Y, int Z) { return node_zcoord[Z];}
  const double &getXC(int X, int Y, int Z) { return center_xcoord[X];}
  const double &getYC(int X, int Y, int Z) { return center_ycoord[Y];}
  const double &getZC(int X, int Y, int Z) { return center_zcoord[Z];}
  //
  double getXstart() { return (xStart); }
  double getXend() { return (xEnd); }
  double getYstart() { return (yStart); }
  double getYend() { return (yEnd); } 
  double getZstart() { return (zStart); }
  double getZend() { return (zEnd); }
  double getInvVOL() { return (invVOL); }
};

typedef Grid3DCU Grid;

#endif
