/*******************************************************************************************************
  Grid3DCUCU.h  -  uniform cartesian 3D local grid for each process, including che guard cells
  -------------------

 *******************************************************************************************************/

#ifndef GRID3DCU_H
#define GRID3DCU_H

#include "arraysfwd.h"
#include "ipicfwd.h"
#include "math.h" // for floor

class VirtualTopology3D;
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
  void derivN(arr3_double derN,
    const_arr4_double scFieldC, int ns, int dir);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradC2N(arr3_double gradXN, arr3_double gradYN, arr3_double gradZN,
    const_arr3_double scFieldC);
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradN2C(arr3_double gradXC, arr3_double gradYC, arr3_double gradZC,
    const_arr3_double scFieldN);
  /** calculate divergence on central points, given a vector field defined on nodes  */
  void divN2C(arr3_double divC,
    const_arr3_double vecFieldXN,
    const_arr3_double vecFieldYN,
    const_arr3_double vecFieldZN);
  /** calculate divergence on nodes, given a vector field defined on central points  */
  void divC2N(arr3_double divN,
    const_arr3_double vecFieldXC,
    const_arr3_double vecFieldYC,
    const_arr3_double vecFieldZC);
  /** calculate curl on nodes, given a vector field defined on central points  */
  void curlC2N(arr3_double curlXN, arr3_double curlYN,
    arr3_double curlZN,
    const_arr3_double vecFieldXC,
    const_arr3_double vecFieldYC,
    const_arr3_double vecFieldZC);
  /** calculate curl on central points, given a vector field defined on nodes  */
  void curlN2C(arr3_double curlXC, arr3_double curlYC, arr3_double curlZC,
    const_arr3_double vecFieldXN,
    const_arr3_double vecFieldYN,
    const_arr3_double vecFieldZN);

  /** calculate divergence on central points, given a Tensor field defined on nodes  */
  void divSymmTensorN2C(arr3_double divCX, arr3_double divCY, arr3_double divCZ,
    const_arr4_double pXX,
    const_arr4_double pXY,
    const_arr4_double pXZ,
    const_arr4_double pYY,
    const_arr4_double pYZ,
    const_arr4_double pZZ, int ns);

  /** calculate laplacian on nodes, given a scalar field defined on nodes */
  void lapN2N(arr3_double lapN,
    const_arr3_double scFieldN, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
  void lapC2Cpoisson(arr3_double lapC,
    arr3_double scFieldC, VirtualTopology3D * vct);
  /** calculate laplacian on central points, given a scalar field defined on central points */
  void lapC2C(arr3_double lapC,
    const_arr3_double scFieldC, VirtualTopology3D * vct);

  /** calculate divergence on boundaries */
  void divBCleft(arr3_double divBC,
    const_arr3_double vectorX,
    const_arr3_double vectorY,
    const_arr3_double vectorZ, int leftActiveNode, int dirDER);
  /** calculate divergence on boundaries */
  void divBCright(arr3_double divBC,
    const_arr3_double vectorX,
    const_arr3_double vectorY,
    const_arr3_double vectorZ, int rightActiveNode, int dirDER);
  /** calculate derivative on boundaries */
  void derBC(arr3_double derBC,
    const_arr3_double vector, int leftActiveNode, int dirDER);


  /** interpolate on nodes from central points */
  void interpC2N(arr3_double vecFieldN, const_arr3_double vecFieldC);
  /** interpolate on central points from nodes */
  void interpN2C(arr3_double vecFieldC, const_arr3_double vecFieldN);
  /** interpolate on central points from nodes */
  void interpN2C(arr4_double vecFieldC, int ns, const_arr4_double vecFieldN);

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
  /** index of last cell including ghost cells */
  // (precomputed for speed)
  int cxlast; // nxc-1;
  int cylast; // nyc-1;
  int czlast; // nzc-1;
  /** node coordinate */
  pfloat *pfloat_node_xcoord;
  pfloat *pfloat_node_ycoord;
  pfloat *pfloat_node_zcoord;
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
  double get_invdx() { return (invdx); }
  double get_invdy() { return (invdy); }
  double get_invdz() { return (invdz); }
  //
  // coordinate accessors
  //
  // calculated equivalents (preferred for accelerator?):
  double calcXN(int X)const{ return xStart+(X-1)*dx;}
  double calcYN(int Y)const{ return yStart+(Y-1)*dy;}
  double calcZN(int Z)const{ return zStart+(Z-1)*dz;}
  const pfloat &get_pfloat_XN(int X)const{ return pfloat_node_xcoord[X];}
  const pfloat &get_pfloat_YN(int Y)const{ return pfloat_node_ycoord[Y];}
  const pfloat &get_pfloat_ZN(int Z)const{ return pfloat_node_zcoord[Z];}
  const double &getXN(int X)const{ return node_xcoord[X];}
  const double &getYN(int Y)const{ return node_ycoord[Y];}
  const double &getZN(int Z)const{ return node_zcoord[Z];}
  const double &getXC(int X)const{ return center_xcoord[X];}
  const double &getYC(int Y)const{ return center_ycoord[Y];}
  const double &getZC(int Z)const{ return center_zcoord[Z];}
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

  // inline methods to calculate mesh cell and weights.
  static void get_weights(double weights[8],
    double w0x, double w0y, double w0z,
    double w1x, double w1y, double w1z)
  {
    // which of the following is faster?
    //
    // this:
    //
    //const double weight00 = w0x*w0y;
    //const double weight01 = w0x*w1y;
    //const double weight10 = w1x*w0y;
    //const double weight11 = w1x*w1y;
    //weights[0] = weight00*w0z; // weight000
    //weights[1] = weight00*w1z; // weight001
    //weights[2] = weight01*w0z; // weight010
    //weights[3] = weight01*w1z; // weight011
    //weights[4] = weight10*w0z; // weight100
    //weights[5] = weight10*w1z; // weight101
    //weights[6] = weight11*w0z; // weight110
    //weights[7] = weight11*w1z; // weight111
    //
    // or this:
    //
    weights[0] = w0x*w0y*w0z; // weight000
    weights[1] = w0x*w0y*w1z; // weight001
    weights[2] = w0x*w1y*w0z; // weight010
    weights[3] = w0x*w1y*w1z; // weight011
    weights[4] = w1x*w0y*w0z; // weight100
    weights[5] = w1x*w0y*w1z; // weight101
    weights[6] = w1x*w1y*w0z; // weight110
    weights[7] = w1x*w1y*w1z; // weight111
  }
  void get_cell_coordinates(
    int& cx, int& cy, int& cz,
    double xpos, double ypos, double zpos)
  {
      // xStart marks start of domain excluding ghosts
      const double rel_xpos = xpos - xStart;
      const double rel_ypos = ypos - yStart;
      const double rel_zpos = zpos - zStart;
      // cell position minus 1 (due to ghost cells)
      const double cxm1_pos = rel_xpos * invdx;
      const double cym1_pos = rel_ypos * invdy;
      const double czm1_pos = rel_zpos * invdz;
      cx = 1 + int(floor(cxm1_pos));
      cy = 1 + int(floor(cym1_pos));
      cz = 1 + int(floor(czm1_pos));
  }
  void make_cell_coordinates_safe(int& cx, int& cy, int& cz)
  {
    // if the cell is outside the domain, then treat it as
    // in the nearest ghost cell.
    //
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    if (cz < 0) cz = 0;
    if (cx > cxlast) cx = cxlast; //nxc-1;
    if (cy > cylast) cy = cylast; //nyc-1;
    if (cz > czlast) cz = czlast; //nzc-1;
  }
  void get_safe_cell_coordinates(
    int& cx, int& cy, int& cz,
    double x, double y, double z)
  {
    get_cell_coordinates(cx,cy,cz,x,y,z);
    make_cell_coordinates_safe(cx,cy,cz);
  }
  void get_safe_cell_and_weights(
    double xpos, double ypos, double zpos,
    int &cx, int& cy, int& cz,
    double weights[8])
  {
    //convert_xpos_to_cxpos(xpos,ypos,zpos,cxpos,cypos,czpos);
    // xStart marks start of domain excluding ghosts
    const double rel_xpos = xpos - xStart;
    const double rel_ypos = ypos - yStart;
    const double rel_zpos = zpos - zStart;
    // cell position minus 1 (due to ghost cells)
    const double cxm1_pos = rel_xpos * invdx;
    const double cym1_pos = rel_ypos * invdy;
    const double czm1_pos = rel_zpos * invdz;
    //
    cx = 1 + int(floor(cxm1_pos));
    cy = 1 + int(floor(cym1_pos));
    cz = 1 + int(floor(czm1_pos));
  
    make_cell_coordinates_safe(cx,cy,cz);
  
    // fraction of the distance from the right of the cell
    const double w1x = cx - cxm1_pos;
    const double w1y = cy - cym1_pos;
    const double w1z = cz - czm1_pos;
    // fraction of distance from the left
    const double w0x = 1.-w1x;
    const double w0y = 1.-w1y;
    const double w0z = 1.-w1z;

    get_weights(weights, w0x, w0y, w0z, w1x, w1y, w1z);
  }
  void get_safe_cell_and_weights(double xpos[3], int cx[3], double weights[8])
  {
    get_safe_cell_and_weights(xpos[0],xpos[1],xpos[2],cx[0],cx[1],cx[2],weights);
  }
};

typedef Grid3DCU Grid;

#endif
