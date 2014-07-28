/*******************************************************************************************************
  Grid3DCUCU.h  -  uniform cartesian 3D local grid for each process, including che guard cells
  -------------------

 *******************************************************************************************************/

#ifndef GRID3DCU_H
#define GRID3DCU_H

#include <iostream>

#include "Grid.h"
#include "Collective.h"
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
  Grid3DCU(Collective * col, VirtualTopology3D * vct);
  /** destructor */
  ~Grid3DCU();
  /** allocate grid arrays for this domain */
  void allocate(Collective * ptC, VirtualTopology3D * ptVCT);
  /** deallocate grid arrays for this domain */
  void deallocate();
  /** print grid info */
  void print(VirtualTopology3D * ptVCT);
  /** calculate a derivative along a direction on nodes */
  void derivN(double ***derN, double ****scFieldC, int ns, int dir);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradN2C(double ***gradXC, double ***gradYC, double ***gradZC, double ***scFieldN);
  /** calculate divergence on central points, given a vector field defined on nodes  */
  void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN);
  /** calculate divergence on nodes, given a vector field defined on central points  */
  void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC);
  /** calculate curl on nodes, given a vector field defined on central points  */
  void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC);
  /** calculate curl on central points, given a vector field defined on nodes  */
  void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN);

  /** calculate divergence on central points, given a Tensor field defined on nodes  */
  void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns);

  /** calculate laplacian on nodes, given a scalar field defined on nodes */
  void lapN2N(double ***lapN, double ***scFieldN, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
  void lapC2Cpoisson(double ***lapC, double ***scFieldC, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points */
  void lapC2C(double ***lapC, double ***scFieldC, VirtualTopology3D * vct);

  /** calculate divergence on boundaries */
  void divBCleft(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER);
  /** calculate divergence on boundaries */
  void divBCright(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int rightActiveNode, int dirDER);
  /** calculate derivative on boundaries */
  void derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER);


  /** interpolate on nodes from central points */
  void interpC2N(double ****vecFieldN, int ns, double ****vecFieldC);
  /** interpolate on nodes from central points */
  void interpC2N(double ***vecFieldN, double ***vecFieldC);
  /** interpolate on central points from nodes */
  void interpN2C(double ***vecFieldC, double ***vecFieldN);
  /** interpolate on central points from nodes */
  void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN);

  /** return nxc */
  int getNXC();
  /** return nxn */
  int getNXN();
  /** return nyc */
  int getNYC();
  /** return nyn */
  int getNYN();
  /** return nzc */
  int getNZC();
  /** return nzn */
  int getNZN();
  /** return dx */
  double getDX();
  /** return dy */
  double getDY();
  /** return dz */
  double getDZ();
  /** get xn(X,Y,Z) */
  double &getXN(int indexX, int indexY, int indexZ);
  /** get yn(X,Y,Z) */
  double &getYN(int indexX, int indexY, int indexZ);
  /** get zn(X,Y,Z) */
  double &getZN(int indexX, int indexY, int indexZ);
  /** get the whole vector of nodes*/
  double ****getN();
  /** get xc(X,Y,Z) */
  double &getXC(int indexX, int indexY, int indexZ);
  /** get yc(X,Y,Z) */
  double &getYC(int indexX, int indexY, int indexZ);
  /** get zc(X,Y,Z) */
  double &getZC(int indexX, int indexY, int indexZ);
  /** get Xstart */
  double getXstart();
  /** get Xend */
  double getXend();
  /** get Ystart */
  double getYstart();
  /** get Yend */
  double getYend();
  /** get Zstart */
  double getZstart();
  /** get Zend */
  double getZend();
  /** get the inverse of volume */
  double getInvVOL();

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
  double ****node_coordinate;
  /** center coordinate */
  double ****center_coordinate;
  /** local grid boundaries coordinate  */
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;

};

typedef Grid3DCU Grid;

#endif
