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
  /** calculate curl on central points, given a vector field defined on nodes
      calculate ghost cells also using ghost node info if at buondary **/
  void curlN2C_Ghost(VirtualTopology3D * vct, double ***curlXC, double ***curlYC, double ***curlZC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN);

  /** calculate divergence on central points, given a Tensor field defined on nodes  */
  void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns);

  /** calculate laplacian on nodes, given a scalar field defined on nodes */
  void lapN2N(double ***lapN, double ***scFieldN, VirtualTopology3D * vct);
  void lapN2N_mlmd(double ***lapN, double ***scFieldN, VirtualTopology3D * vct);
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
  void interpN2C_GC(double ***vecFieldC, double ***vecFieldN);
  /** interpolate on central points from nodes */
  void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN);

  void interpN2C_GC(double ****vecFieldC, int ns, double ****vecFieldN);

  void interpN2C_ActiveCell(double ***vecFieldC, double ***vecFieldN, VirtualTopology3D * vct) ;


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

  const double &getXN(int X, int Y, int Z) { return node_xcoord[X];}
  const double &getYN(int X, int Y, int Z) { return node_ycoord[Y];}
  const double &getZN(int X, int Y, int Z) { return node_zcoord[Z];}
  const double &getXC(int X, int Y, int Z) { return center_xcoord[X];}
  const double &getYC(int X, int Y, int Z) { return center_ycoord[Y];}
  const double &getZC(int X, int Y, int Z) { return center_zcoord[Z];}

  /* mlmd: coordinate on your parent grid 
     NB: i want it to be able to manage also negative indexes or indexes > nxn/ nyn/ nzn 
     (for the phase 1 of particle init BC) */
  double getXN_P(int X, int Y, int Z); 
  double getYN_P(int X, int Y, int Z); 
  double getZN_P(int X, int Y, int Z); 
  /* end mlmd: coordinate on your parent grid */
  /* mlmd: coordinate of centers on the parent grid */
  double getXC_P(int X, int Y, int Z); 
  double getYC_P(int X, int Y, int Z); 
  double getZC_P(int X, int Y, int Z); 
  /* end mlmd: coordinate of centers on the parent grid */

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

  /** get xStart_GC **/
  double getxStart_GC();
  /** get yStart_GC **/
  double getyStart_GC();
  /** get zStart_GC **/
  double getzStart_GC();
  /** get xEnd_GC **/
  double getxEnd_GC();
  /** get yEnd_GC **/
  double getyEnd_GC();
  /** get zEnd_GC **/
  double getzEnd_GC();

  /*! mlmd specific functions */
  int getNumGrid(){return numGrid;}
  /*! return your coordinates of origin on the parent grid */
  double getOx(){return Ox;} double getOy(){return Oy;} double getOz(){return Oz;}
  
  double getDx_mlmd(int NG){return dx_mlmd[NG];}
  double getDy_mlmd(int NG){return dy_mlmd[NG];}
  double getDz_mlmd(int NG){return dz_mlmd[NG];}

  /*! return your coordinates of origin in the system */
  double getOx_SW(){return Ox_SW;} double getOy_SW(){return Oy_SW;} double getOz_SW(){return Oz_SW;}

  /*! return parentLenX, parentLenY, parentLenZ */
  double getparentLenX(){return parentLenX;}
  double getparentLenY(){return parentLenY;}
  double getparentLenZ(){return parentLenZ;}

  /*** pay exceptional attention to this description ***/
  /* nx, ny, nz: index in the current grid, which is a child*/
  /* returns the rank IN THE PARENT-CHILD communicator of the coarse grid core where the point is hosted    
     only the active part of the parent grid is examined*/
  int getParentRankFromGridPoint(VirtualTopology3D * vct, int xn, int yn, int zn);
  int getParentRankFromGridCenter(VirtualTopology3D * vct, int xn, int yn, int zn);
  void RGBCExploreDirection(VirtualTopology3D *vct,string FACE, int DIR, int i0_s, int i0_e, int i1, int i2, double *SPXperC, double *SPYperC, double *SPZperC, int *NPperC, int *rank, int* Ncores, int *IndexFirstPointperC);

  void RGBCExploreDirection_Centers(VirtualTopology3D *vct,string FACE, int DIR, int i0_s, int i0_e, int i1, int i2, double *SPXperC, double *SPYperC, double *SPZperC, int *NPperC, int *rank, int* Ncores, int *IndexFirstPointperC);
  /** given the rank N on the PARENT-CHILD communicator of the core,
      it returns it physical extension of this PARENT core
      -- at the moment, used for debug only --**/
  void getParentLimits(VirtualTopology3D *vct, int N, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax);
  /*! end mlmd specific functions */

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


  double *node_xcoord; /** Node X coordinate */
  double *node_ycoord; /** Node Y coordinate */
  double *node_zcoord; /** Node Z coordinate */

  double *center_xcoord; /** Cell center X coordinate */
  double *center_ycoord; /** Cell center Y coordinate */
  double *center_zcoord; /** Cell center Z coordinate */

  /** local grid boundaries coordinate  */
  /** mlmd: checked: they mark the active part of the grid **/
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;

  /** include ghost cells if boundary core- to be used only in communicateRepopulatedParticles **/
  double xStart_GC, xEnd_GC, yStart_GC, yEnd_GC, zStart_GC, zEnd_GC;

  /*! mlmd specific variables */
  /*! total number of grids in the mlmd hierarchy */
  int Ngrids;
  /*! number of the current grid in the mlmd grid hierarchy */
  int numGrid;
  /* coordinates of the origin on the PARENT grid */
  double Ox, Oy, Oz;

  /* coordinates of the grid of origin System-wide */
  double Ox_SW, Oy_SW, Oz_SW;

  /* portion of ACTIVE grid hosted in each parent core                                           
     -- equivalent of xEnd - xStart on the parent -- used in getParentRankFromGridPoint */
  double parentLenX;
  double parentLenY;
  double parentLenZ;

  // resolution and length of all the grids in the MLMD system
  double * dx_mlmd;
  double * dy_mlmd;
  double * dz_mlmd;

  /*
  double *Lx_mlmd;
  double *Ly_mlmd;
  double *Lz_mlmd;*/

  /*! end mlmd specific variables */
};

typedef Grid3DCU Grid;

#endif
