/***************************************************************************
  Grid.h  -  Abstract Base class for grids
  -------------------
begin                : Wed Jun 2 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef Grid_H
#define Grid_H
/**
 *  Abstract base class for different kinds of grid.
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
//class Grid {
//public:
//  /** print grid info */
//  virtual void print(VirtualTopology3D * ptVCT) = 0;
//  /** calculate gradient on nodes, given a scalar field defined on central points  */
//  virtual void gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC) = 0;
//  /** calculate gradient on nodes, given a scalar field defined on central points  */
//  virtual void gradN2C(double ***gradXC, double ***gradYC, double ***gradZC, double ***scFieldN) = 0;
//  /** calculate divergence on central points, given a vector field defined on nodes  */
//  virtual void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN) = 0;
//  /** calculate divergence on nodes, given a vector field defined on central points  */
//  virtual void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC) = 0;
//  /** calculate curl on nodes, given a vector field defined on central points  */
//  virtual void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC) = 0;
//  /** calculate curl on central points, given a vector field defined on nodes  */
//  virtual void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN) = 0;
//  /** calculate laplacian on nodes, given a scalar field defined on nodes */
//  virtual void lapN2N(double ***lapN, double ***scFieldN, VirtualTopology3D * vct) = 0;
//  /** calculate laplacian on central points, given a scalar field defined on central points */
//  virtual void lapC2C(double ***lapC, double ***scFieldC, VirtualTopology3D * vct) = 0;
//  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
//  virtual void lapC2Cpoisson(double ***lapC, double ***scFieldC, VirtualTopology3D * vct) = 0;
//  /** calculate divergence on central points, given a Tensor field defined on nodes  */
//  virtual void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns) = 0;
//  /** calculate divergence on boudaries */
//  virtual void divBCleft(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dir) = 0;
//  /** calculate divergence on boudaries */
//  virtual void divBCright(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int rightActiveNode, int dir) = 0;
//  /** calculate derivative on boundaries */
//  virtual void derBC(double ***derBC, double ***vector, int leftActiveNode, int dir) = 0;
//  /** interpolate on nodes from central points */
//  virtual void interpC2N(double ***vecFieldN, double ***vecFieldC) = 0;
//  /** interpolate on central points from nodes */
//  virtual void interpN2C(double ***vecFieldC, double ***vecFieldN) = 0;
//  /** interpolate on central points from nodes */
//  virtual void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN) = 0;
//  /** return nxc */
//  virtual int getNXC() = 0;
//  /** return nxn */
//  virtual int getNXN() = 0;
//  /** return nyc */
//  virtual int getNYC() = 0;
//  /** return nyn */
//  virtual int getNYN() = 0;
//  /** return nzc */
//  virtual int getNZC() = 0;
//  /** return nzn */
//  virtual int getNZN() = 0;
//  /** return dx */
//  virtual double getDX() = 0;
//  /** return dy */
//  virtual double getDY() = 0;
//  /** return dz */
//  virtual double getDZ() = 0;
//  /** get xn(X,Y,Z) */
//  virtual double &getXN(int indexX, int indexY, int indexZ) = 0;
//  /** get yn(X,Y,Z) */
//  virtual double &getYN(int indexX, int indexY, int indexZ) = 0;
//  /** get zn(X,Y,Z) */
//  virtual double &getZN(int indexX, int indexY, int indexZ) = 0;
//  /** get xc(X,Y,Z) */
//  virtual double &getXC(int indexX, int indexY, int indexZ) = 0;
//  /** get yc(X,Y,Z) */
//  virtual double &getYC(int indexX, int indexY, int indexZ) = 0;
//  /** get the whole vector xc*/
//  virtual double ***getXC() = 0;
//  /** get the whole vector yc*/
//  virtual double ***getYC() = 0;
//  /** get the whole vector zc*/
//  virtual double ***getZC() = 0;
//  /** get zc(X,Y,Z) */
//  virtual double &getZC(int indexX, int indexY, int indexZ) = 0;
//  /** get Xstart */
//  virtual double getXstart() = 0;
//  /** get Xend */
//  virtual double getXend() = 0;
//  /** get Ystart */
//  virtual double getYstart() = 0;
//  /** get Yend */
//  virtual double getYend() = 0;
//  /** get Zstart */
//  virtual double getZstart() = 0;
//  /** get Zend */
//  virtual double getZend() = 0;
//  /** get the inverse of Volume */
//  virtual double getInvVOL() = 0;
//
//};
#include "Grid3DCU.h"
#endif
