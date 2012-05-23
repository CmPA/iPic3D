 /*******************************************************************************************************
        Grid3DCUCU.h  -  uniform cartesian 3D local grid for each process, including che guard cells
                             -------------------
   
 *******************************************************************************************************/

#ifndef GRID3DCU_H
#define GRID3DCU_H

#include <iostream>

#include "Grid.h"
#include "../communication/ComInterpNodes3D.h"
#include "../communication/ComNodes3D.h"
#include "../utility/Alloc.h"

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
class Grid3DCU : public Grid {
  public:
      /** constructor */
      Grid3DCU(CollectiveIO *col, VirtualTopology3D *vct);
      /** destructor */
      ~Grid3DCU();
      /** allocate grid arrays for this domain */
      void allocate(CollectiveIO *ptC, VirtualTopology3D *ptVCT);
      /** deallocate grid arrays for this domain */
      void deallocate();
      /** print grid info */
      void print(VirtualTopology3D* ptVCT);
      /** calculate a derivative along a direction on nodes */
      void derivN(double ***derN, double ****scFieldC, int ns, int dir);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradN2C(double ***gradXC, double ***gradYC, double ***gradZC, double ***scFieldN);
      /** calculate divergence on central points, given a vector field defined on nodes  */
      void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN);
      /** calculate divergence on nodes, given a vector field defined on central points  */
      void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC);
      /** calculate curl on nodes, given a vector field defined on central points  */
      void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC);
      /** calculate curl on central points, given a vector field defined on nodes  */
      void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN);

      /** calculate divergence on central points, given a Tensor field defined on nodes  */
      void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);
      
      /** calculate laplacian on nodes, given a scalar field defined on nodes */
      void lapN2N(double*** lapN,double ***scFieldN,VirtualTopology3D *vct);
	  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
      void lapC2Cpoisson(double*** lapC,double ***scFieldC,VirtualTopology3D *vct);
      /** calculate laplacian on central points, given a scalar field defined on central points */
      void lapC2C(double*** lapC,double ***scFieldC,VirtualTopology3D *vct);

      /** calculate divergence on boundaries */
      void divBCleft(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER);
      /** calculate divergence on boundaries */
      void divBCright(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int rightActiveNode, int dirDER);
      /** calculate derivative on boundaries */
      void derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER);
      
      
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
      double &getXN(int indexX, int indexY,int indexZ);
      /** get yn(X,Y,Z) */
      double &getYN(int indexX, int indexY,int indexZ);
      /** get zn(X,Y,Z) */
      double &getZN(int indexX, int indexY,int indexZ);
      /** get xc(X,Y,Z) */
      double &getXC(int indexX, int indexY,int indexZ);
      /** get yc(X,Y,Z) */
      double &getYC(int indexX, int indexY,int indexZ);
      /** get zc(X,Y,Z) */
      double &getZC(int indexX, int indexY,int indexZ);
      /** get the whole vector xc*/
      double*** getXC();
      /** get the whole vector yc*/
      double*** getYC();
      /** get the whole vector zc*/
      double*** getZC();
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

  /////////////   PRIVATE VARIABLES //////////////
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
	 double **** node_coordinate;
	 /** center coordinate */
	 double **** center_coordinate;
     /** local grid boundaries coordinate  */
     double xStart, xEnd, yStart, yEnd, zStart, zEnd;
         
};
/** constructor */
inline Grid3DCU::Grid3DCU(CollectiveIO* col, VirtualTopology3D* vct){
   // add 2 for the guard cells
   nxc = (col->getNxc())/(vct->getXLEN()) + 2;
   nyc = (col->getNyc())/(vct->getYLEN()) + 2;
   nzc = (col->getNzc())/(vct->getZLEN()) + 2;
   nxn = nxc + 1;
   nyn = nyc + 1;
   nzn = nzc + 1;
   dx  = col->getLx()/col->getNxc();
   dy  = col->getLy()/col->getNyc();
   dz  = col->getLz()/col->getNzc();
   invVOL = 1.0/(dx*dy*dz);
   invdx = 1.0/dx;
   invdy = 1.0/dy;
   invdz = 1.0/dz;

   // local grid dimensions and boundaries of active nodes
   xStart = vct->getCoordinates(0)*(col->getLx()/(double) vct->getXLEN());
  
   xEnd   = xStart + (col->getLx()/(double) vct->getXLEN());
   
   yStart = vct->getCoordinates(1)*(col->getLy()/(double) vct->getYLEN());
   
   yEnd   = yStart + (col->getLy()/(double) vct->getYLEN());
   
   zStart = vct->getCoordinates(2)*(col->getLz()/(double) vct->getZLEN());
  
   zEnd   = zStart + (col->getLz()/(double) vct->getZLEN());
  
   // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
   node_coordinate = newArr4(double,nxn,nyn,nzn,3); // 0 -> X, 1 -> Y, 2-> Z
   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
        for (int k=0; k < nzn; k++){
          node_coordinate[i][j][k][0] = xStart + (i-1)*dx;
          node_coordinate[i][j][k][1] = yStart + (j-1)*dy;
          node_coordinate[i][j][k][2] = zStart + (k-1)*dz;
		}
     }
   }
  // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
  center_coordinate = newArr4(double,nxc,nyc,nzc,3);
  for (int i=0; i < nxc; i++){
     for (int j=0; j < nyc; j++){
        for (int k=0; k < nzc; k++){
          center_coordinate[i][j][k][0] = .5*(node_coordinate[i][j][k][0] + node_coordinate[i+1][j][k][0]);
          center_coordinate[i][j][k][1] = .5*(node_coordinate[i][j][k][1] + node_coordinate[i][j+1][k][1]);
          center_coordinate[i][j][k][2] = .5*(node_coordinate[i][j][k][2] + node_coordinate[i][j][k+1][2]);
		}
     }
   }
}
/** deallocate the local grid */
inline Grid3DCU::~Grid3DCU(){
   delArr4(node_coordinate,nxn,nyn,nzn);
   delArr4(center_coordinate,nxc,nyc,nzc);
}
/** print the local grid info */
inline void Grid3DCU::print(VirtualTopology3D* ptVCT){
    cout << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ","  << ptVCT->getCoordinates(2) << ")"<< endl;
    cout <<  "Number of cell: -X=" << nxc-2 << " -Y=" << nyc-2 << " -Z=" << nzc-2 << endl;
    cout <<  "Xin = " << node_coordinate[1][1][1][0] << "; Xfin = " << node_coordinate[nxn-2][1][1][0] << endl;
    cout <<  "Yin = " << node_coordinate[1][1][1][1] << "; Yfin = " << node_coordinate[1][nyn-2][1][1] << endl;
    cout <<  "Zin = " << node_coordinate[1][1][1][2] << "; Zfin = " << node_coordinate[1][1][nzn-2][2] << endl;
    cout << endl;
    
}
/** calculate gradient on nodes, given a scalar field defined on central points  */
inline void Grid3DCU::gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC){
   for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
        for (register int k=1;k < nzn-1;k++){
          gradXN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i-1][j][k])*invdx +  .25*(scFieldC[i][j][k-1] - scFieldC[i-1][j][k-1])*invdx + .25*(scFieldC[i][j-1][k] - scFieldC[i-1][j-1][k])*invdx + .25*(scFieldC[i][j-1][k-1] - scFieldC[i-1][j-1][k-1])*invdx;
          gradYN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i][j-1][k])*invdy +  .25*(scFieldC[i][j][k-1] - scFieldC[i][j-1][k-1])*invdy + .25*(scFieldC[i-1][j][k] - scFieldC[i-1][j-1][k])*invdy + .25*(scFieldC[i-1][j][k-1] - scFieldC[i-1][j-1][k-1])*invdy;
          gradZN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i][j][k-1])*invdz +  .25*(scFieldC[i-1][j][k] - scFieldC[i-1][j][k-1])*invdz + .25*(scFieldC[i][j-1][k] - scFieldC[i][j-1][k-1])*invdz + .25*(scFieldC[i-1][j-1][k] - scFieldC[i-1][j-1][k-1])*invdz;
        }
  

}
/** calculate gradient on nodes, given a scalar field defined on central points  */
inline void Grid3DCU::gradN2C(double ***gradXC, double ***gradYC, double ***gradZC, double ***scFieldN){
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++){
          gradXC[i][j][k] = .25*(scFieldN[i+1][j][k] - scFieldN[i][j][k])*invdx +  .25*(scFieldN[i+1][j][k+1] - scFieldN[i][j][k+1])*invdx + .25*(scFieldN[i+1][j+1][k] - scFieldN[i][j+1][k])*invdx + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i][j+1][k+1])*invdx;
          gradYC[i][j][k] = .25*(scFieldN[i][j+1][k] - scFieldN[i][j][k])*invdy +  .25*(scFieldN[i][j+1][k+1] - scFieldN[i][j][k+1])*invdy + .25*(scFieldN[i+1][j+1][k] - scFieldN[i+1][j][k])*invdy + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i+1][j][k+1])*invdy;
          gradZC[i][j][k] = .25*(scFieldN[i][j][k+1] - scFieldN[i][j][k])*invdz +  .25*(scFieldN[i+1][j][k+1] - scFieldN[i+1][j][k])*invdz + .25*(scFieldN[i][j+1][k+1] - scFieldN[i][j+1][k])*invdz + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i+1][j+1][k])*invdz;
       }
 }

/** calculate divergence on central points, given a vector field defined on nodes  */
inline void Grid3DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
  double compX;
  double compY;
  double compZ;
   for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++){
          compX = .25*(vecFieldXN[i+1][j][k] - vecFieldXN[i][j][k])*invdx +  .25*(vecFieldXN[i+1][j][k+1] - vecFieldXN[i][j][k+1])*invdx + .25*(vecFieldXN[i+1][j+1][k] - vecFieldXN[i][j+1][k])*invdx + .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i][j+1][k+1])*invdx;
          compY = .25*(vecFieldYN[i][j+1][k] - vecFieldYN[i][j][k])*invdy +  .25*(vecFieldYN[i][j+1][k+1] - vecFieldYN[i][j][k+1])*invdy + .25*(vecFieldYN[i+1][j+1][k] - vecFieldYN[i+1][j][k])*invdy + .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i+1][j][k+1])*invdy;
          compZ = .25*(vecFieldZN[i][j][k+1] - vecFieldZN[i][j][k])*invdz +  .25*(vecFieldZN[i+1][j][k+1] - vecFieldZN[i+1][j][k])*invdz + .25*(vecFieldZN[i][j+1][k+1] - vecFieldZN[i][j+1][k])*invdz + .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i+1][j+1][k])*invdz;
          divC[i][j][k] = compX + compY + compZ;
       }
}
/** calculate divergence on central points, given a Tensor field defined on nodes  */
inline void Grid3DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
  double comp1Z, comp2Z, comp3Z;
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++){
          comp1X = .25*(pXX[ns][i+1][j][k] - pXX[ns][i][j][k])*invdx +  .25*(pXX[ns][i+1][j][k+1] - pXX[ns][i][j][k+1])*invdx + .25*(pXX[ns][i+1][j+1][k] - pXX[ns][i][j+1][k])*invdx + .25*(pXX[ns][i+1][j+1][k+1] - pXX[ns][i][j+1][k+1])*invdx;
          comp2X = .25*(pXY[ns][i+1][j][k] - pXY[ns][i][j][k])*invdx +  .25*(pXY[ns][i+1][j][k+1] - pXY[ns][i][j][k+1])*invdx + .25*(pXY[ns][i+1][j+1][k] - pXY[ns][i][j+1][k])*invdx + .25*(pXY[ns][i+1][j+1][k+1] - pXY[ns][i][j+1][k+1])*invdx;
          comp3X = .25*(pXZ[ns][i+1][j][k] - pXZ[ns][i][j][k])*invdx +  .25*(pXZ[ns][i+1][j][k+1] - pXZ[ns][i][j][k+1])*invdx + .25*(pXZ[ns][i+1][j+1][k] - pXZ[ns][i][j+1][k])*invdx + .25*(pXZ[ns][i+1][j+1][k+1] - pXZ[ns][i][j+1][k+1])*invdx;
          comp1Y = .25*(pXY[ns][i][j+1][k] - pXY[ns][i][j][k])*invdy +  .25*(pXY[ns][i][j+1][k+1] - pXY[ns][i][j][k+1])*invdy + .25*(pXY[ns][i+1][j+1][k] - pXY[ns][i+1][j][k])*invdy + .25*(pXY[ns][i+1][j+1][k+1] - pXY[ns][i+1][j][k+1])*invdy;
          comp2Y = .25*(pYY[ns][i][j+1][k] - pYY[ns][i][j][k])*invdy +  .25*(pYY[ns][i][j+1][k+1] - pYY[ns][i][j][k+1])*invdy + .25*(pYY[ns][i+1][j+1][k] - pYY[ns][i+1][j][k])*invdy + .25*(pYY[ns][i+1][j+1][k+1] - pYY[ns][i+1][j][k+1])*invdy;
          comp3Y = .25*(pYZ[ns][i][j+1][k] - pYZ[ns][i][j][k])*invdy +  .25*(pYZ[ns][i][j+1][k+1] - pYZ[ns][i][j][k+1])*invdy + .25*(pYZ[ns][i+1][j+1][k] - pYZ[ns][i+1][j][k])*invdy + .25*(pYZ[ns][i+1][j+1][k+1] - pYZ[ns][i+1][j][k+1])*invdy;
          comp1Z = .25*(pXZ[ns][i][j][k+1] - pXZ[ns][i][j][k])*invdz +  .25*(pXZ[ns][i+1][j][k+1] - pXZ[ns][i+1][j][k])*invdz + .25*(pXZ[ns][i][j+1][k+1] - pXZ[ns][i][j+1][k])*invdz + .25*(pXZ[ns][i+1][j+1][k+1] - pXZ[ns][i+1][j+1][k])*invdz;
          comp2Z = .25*(pYZ[ns][i][j][k+1] - pYZ[ns][i][j][k])*invdz +  .25*(pYZ[ns][i+1][j][k+1] - pYZ[ns][i+1][j][k])*invdz + .25*(pYZ[ns][i][j+1][k+1] - pYZ[ns][i][j+1][k])*invdz + .25*(pYZ[ns][i+1][j+1][k+1] - pYZ[ns][i+1][j+1][k])*invdz;
          comp3Z = .25*(pZZ[ns][i][j][k+1] - pZZ[ns][i][j][k])*invdz +  .25*(pZZ[ns][i+1][j][k+1] - pZZ[ns][i+1][j][k])*invdz + .25*(pZZ[ns][i][j+1][k+1] - pZZ[ns][i][j+1][k])*invdz + .25*(pZZ[ns][i+1][j+1][k+1] - pZZ[ns][i+1][j+1][k])*invdz;
          divCX[i][j][k] = comp1X + comp2X + comp3X;
          divCY[i][j][k] = comp1Y + comp2Y + comp3Y;
          divCZ[i][j][k] = comp1Z + comp2Z + comp3Z;
       }  
}
/** calculate divergence on nodes, given a vector field defined on central points  */
inline void Grid3DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC){
  double compX;
  double compY;
  double compZ;
  for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
        for (register int k=1;k < nzn-1;k++){
          compX = .25*(vecFieldXC[i][j][k] - vecFieldXC[i-1][j][k])*invdx +  .25*(vecFieldXC[i][j][k-1] - vecFieldXC[i-1][j][k-1])*invdx + .25*(vecFieldXC[i][j-1][k] - vecFieldXC[i-1][j-1][k])*invdx + .25*(vecFieldXC[i][j-1][k-1] - vecFieldXC[i-1][j-1][k-1])*invdx;
          compY = .25*(vecFieldYC[i][j][k] - vecFieldYC[i][j-1][k])*invdy +  .25*(vecFieldYC[i][j][k-1] - vecFieldYC[i][j-1][k-1])*invdy + .25*(vecFieldYC[i-1][j][k] - vecFieldYC[i-1][j-1][k])*invdy + .25*(vecFieldYC[i-1][j][k-1] - vecFieldYC[i-1][j-1][k-1])*invdy;
          compZ = .25*(vecFieldZC[i][j][k] - vecFieldZC[i][j][k-1])*invdz +  .25*(vecFieldZC[i-1][j][k] - vecFieldZC[i-1][j][k-1])*invdz + .25*(vecFieldZC[i][j-1][k] - vecFieldZC[i][j-1][k-1])*invdz + .25*(vecFieldZC[i-1][j-1][k] - vecFieldZC[i-1][j-1][k-1])*invdz;
          divN[i][j][k] = compX + compY + compZ; 
       }
}
/** calculate curl on nodes, given a vector field defined on central points  */
inline void Grid3DCU::curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
        for (register int k=1;k < nzn-1;k++){
          // curl - X
          compZDY = .25*(vecFieldZC[i][j][k] - vecFieldZC[i][j-1][k])*invdy +  .25*(vecFieldZC[i][j][k-1] - vecFieldZC[i][j-1][k-1])*invdy + .25*(vecFieldZC[i-1][j][k] - vecFieldZC[i-1][j-1][k])*invdy + .25*(vecFieldZC[i-1][j][k-1] - vecFieldZC[i-1][j-1][k-1])*invdy;
          compYDZ = .25*(vecFieldYC[i][j][k] - vecFieldYC[i][j][k-1])*invdz +  .25*(vecFieldYC[i-1][j][k] - vecFieldYC[i-1][j][k-1])*invdz + .25*(vecFieldYC[i][j-1][k] - vecFieldYC[i][j-1][k-1])*invdz + .25*(vecFieldYC[i-1][j-1][k] - vecFieldYC[i-1][j-1][k-1])*invdz;
          // curl - Y
          compXDZ = .25*(vecFieldXC[i][j][k] - vecFieldXC[i][j][k-1])*invdz +  .25*(vecFieldXC[i-1][j][k] - vecFieldXC[i-1][j][k-1])*invdz + .25*(vecFieldXC[i][j-1][k] - vecFieldXC[i][j-1][k-1])*invdz + .25*(vecFieldXC[i-1][j-1][k] - vecFieldXC[i-1][j-1][k-1])*invdz;
          compZDX = .25*(vecFieldZC[i][j][k] - vecFieldZC[i-1][j][k])*invdx +  .25*(vecFieldZC[i][j][k-1] - vecFieldZC[i-1][j][k-1])*invdx + .25*(vecFieldZC[i][j-1][k] - vecFieldZC[i-1][j-1][k])*invdx + .25*(vecFieldZC[i][j-1][k-1] - vecFieldZC[i-1][j-1][k-1])*invdx;
          // curl - Z
          compYDX = .25*(vecFieldYC[i][j][k] - vecFieldYC[i-1][j][k])*invdx +  .25*(vecFieldYC[i][j][k-1] - vecFieldYC[i-1][j][k-1])*invdx + .25*(vecFieldYC[i][j-1][k] - vecFieldYC[i-1][j-1][k])*invdx + .25*(vecFieldYC[i][j-1][k-1] - vecFieldYC[i-1][j-1][k-1])*invdx;
          compXDY = .25*(vecFieldXC[i][j][k] - vecFieldXC[i][j-1][k])*invdy +  .25*(vecFieldXC[i][j][k-1] - vecFieldXC[i][j-1][k-1])*invdy + .25*(vecFieldXC[i-1][j][k] - vecFieldXC[i-1][j-1][k])*invdy + .25*(vecFieldXC[i-1][j][k-1] - vecFieldXC[i-1][j-1][k-1])*invdy;

          curlXN[i][j][k] = compZDY - compYDZ;
          curlYN[i][j][k] = compXDZ - compZDX;
          curlZN[i][j][k] = compYDX - compXDY;
          }
 

}
/** calculate curl on central points, given a vector field defined on nodes  */
inline void Grid3DCU::curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++){
          // curl - X
          compZDY = .25*(vecFieldZN[i][j+1][k] - vecFieldZN[i][j][k])*invdy +  .25*(vecFieldZN[i][j+1][k+1] - vecFieldZN[i][j][k+1])*invdy + .25*(vecFieldZN[i+1][j+1][k] - vecFieldZN[i+1][j][k])*invdy + .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i+1][j][k+1])*invdy;
          compYDZ = .25*(vecFieldYN[i][j][k+1] - vecFieldYN[i][j][k])*invdz +  .25*(vecFieldYN[i+1][j][k+1] - vecFieldYN[i+1][j][k])*invdz + .25*(vecFieldYN[i][j+1][k+1] - vecFieldYN[i][j+1][k])*invdz + .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i+1][j+1][k])*invdz;
          // curl - Y
          compXDZ = .25*(vecFieldXN[i][j][k+1] - vecFieldXN[i][j][k])*invdz +  .25*(vecFieldXN[i+1][j][k+1] - vecFieldXN[i+1][j][k])*invdz + .25*(vecFieldXN[i][j+1][k+1] - vecFieldXN[i][j+1][k])*invdz + .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i+1][j+1][k])*invdz;
          compZDX = .25*(vecFieldZN[i+1][j][k] - vecFieldZN[i][j][k])*invdx +  .25*(vecFieldZN[i+1][j][k+1] - vecFieldZN[i][j][k+1])*invdx + .25*(vecFieldZN[i+1][j+1][k] - vecFieldZN[i][j+1][k])*invdx + .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i][j+1][k+1])*invdx;
          // curl - Z
          compYDX = .25*(vecFieldYN[i+1][j][k] - vecFieldYN[i][j][k])*invdx +  .25*(vecFieldYN[i+1][j][k+1] - vecFieldYN[i][j][k+1])*invdx + .25*(vecFieldYN[i+1][j+1][k] - vecFieldYN[i][j+1][k])*invdx + .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i][j+1][k+1])*invdx;
          compXDY = .25*(vecFieldXN[i][j+1][k] - vecFieldXN[i][j][k])*invdy +  .25*(vecFieldXN[i][j+1][k+1] - vecFieldXN[i][j][k+1])*invdy + .25*(vecFieldXN[i+1][j+1][k] - vecFieldXN[i+1][j][k])*invdy + .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i+1][j][k+1])*invdy;


          curlXC[i][j][k] = compZDY - compYDZ;
          curlYC[i][j][k] = compXDZ - compZDX;
          curlZC[i][j][k] = compYDX - compXDY;
      }
 


}
/** calculate laplacian on nodes, given a scalar field defined on nodes */
inline void Grid3DCU::lapN2N(double*** lapN,double ***scFieldN,VirtualTopology3D *vct){
   // calculate laplacian as divercence of gradient
   // allocate 3 gradients: defined on central points
   double*** gradXC = newArr3(double,nxc,nyc,nzc);
   double*** gradYC = newArr3(double,nxc,nyc,nzc);
   double*** gradZC = newArr3(double,nxc,nyc,nzc);
 
   gradN2C(gradXC,gradYC,gradZC,scFieldN);
   // communicate with BC
   communicateCenterBC(nxc,nyc,nzc,gradXC,1,1,1,1,1,1,vct);
   communicateCenterBC(nxc,nyc,nzc,gradYC,1,1,1,1,1,1,vct);
   communicateCenterBC(nxc,nyc,nzc,gradZC,1,1,1,1,1,1,vct);
   divC2N(lapN,gradXC,gradYC,gradZC);
   // deallocate
   delArr3(gradXC,nxc,nyc);
   delArr3(gradYC,nxc,nyc);
   delArr3(gradZC,nxc,nyc);
}
/** calculate laplacian on central points, given a scalar field defined on central points */
inline void Grid3DCU::lapC2C(double*** lapC,double ***scFieldC,VirtualTopology3D *vct){
    //  calculate laplacian as divercence of gradient
    // allocate 3 gradients: defined on nodes
    double*** gradXN = newArr3(double,nxn,nyn,nzn);
    double*** gradYN = newArr3(double,nxn,nyn,nzn);
    double*** gradZN = newArr3(double,nxn,nyn,nzn);
    
    gradC2N(gradXN,gradYN,gradZN,scFieldC);
	 if (vct->getYleft_neighbor()==MPI_PROC_NULL){
	   for (int ii=0; ii < nxn;ii++)
	     for (int kk=0; kk < nzn;kk++){
		                   gradXN[ii][0][kk] = 0.0;
                                   gradXN[ii][1][kk] = 0.0;
                                   gradXN[ii][2][kk] = 0.0;
				   gradZN[ii][0][kk] = 0.0;
                                   gradZN[ii][1][kk] = 0.0;
                                   gradZN[ii][2][kk] = 0.0;
				  // gradYN[ii][1][kk] = gradYN[ii][2][kk];
		 }
	}
	if (vct->getYright_neighbor()==MPI_PROC_NULL){
	   for (int ii=0; ii < nxn;ii++)
	     for (int kk=0; kk < nzn;kk++){
                                  gradXN[ii][nyn-1][kk] = 0.0;
				  gradXN[ii][nyn-2][kk] = 0.0;
                                  gradXN[ii][nyn-3][kk] = 0.0;
				  gradZN[ii][nyn-1][kk] = 0.0;
                                  gradZN[ii][nyn-2][kk] = 0.0;
                                  gradZN[ii][nyn-3][kk] = 0.0;  
				  // gradYN[ii][nyn-2][kk]  = gradYN[ii][nyc-3][kk];
		 }
	}
    divN2C(lapC,gradXN,gradYN,gradZN);
    
	delArr3(gradXN,nxn,nyn);
    delArr3(gradYN,nxn,nyn);
    delArr3(gradZN,nxn,nyn);

}    
/** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
inline void Grid3DCU::lapC2Cpoisson(double*** lapC,double ***scFieldC,VirtualTopology3D *vct){
    // communicate first the scFieldC
	communicateCenterBoxStencilBC(nxc,nyc,nzc,scFieldC,1,1,1,1,1,1,vct);
	for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++)
         lapC[i][j][k] = (scFieldC[i-1][j][k] -2*scFieldC[i][j][k] + scFieldC[i+1][j][k])*invdx*invdx + (scFieldC[i][j-1][k] -2*scFieldC[i][j][k] + scFieldC[i][j+1][k])*invdy*invdy + (scFieldC[i][j][k-1] -2*scFieldC[i][j][k] + scFieldC[i][j][k+1])*invdz*invdz;
}                          
/** calculate divergence on  boundaries */
inline void Grid3DCU::divBCleft(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER){
   double compX, compY, compZ;
   switch(dirDER){
     case 0:    //  DIVERGENCE DIRECTION  X
      for (register int j=1;j < nyn-2;j++)
        for (register int k=1;k < nzn-2;k++){
          compX = .25*(vectorX[leftActiveNode+1][j][k] - vectorX[leftActiveNode][j][k])*invdx +  .25*(vectorX[leftActiveNode+1][j][k+1] - vectorX[leftActiveNode][j][k+1])*invdx + .25*(vectorX[leftActiveNode+1][j+1][k] - vectorX[leftActiveNode][j+1][k])*invdx + .25*(vectorX[leftActiveNode+1][j+1][k+1] - vectorX[leftActiveNode][j+1][k+1])*invdx;
          compY = .25*(vectorY[leftActiveNode][j+1][k] - vectorY[leftActiveNode][j][k])*invdy +  .25*(vectorY[leftActiveNode][j+1][k+1] - vectorY[leftActiveNode][j][k+1])*invdy + .25*(vectorY[leftActiveNode+1][j+1][k] - vectorY[leftActiveNode+1][j][k])*invdy + .25*(vectorY[leftActiveNode+1][j+1][k+1] - vectorY[leftActiveNode+1][j][k+1])*invdy;
          compZ = .25*(vectorZ[leftActiveNode][j][k+1] - vectorZ[leftActiveNode][j][k])*invdz +  .25*(vectorZ[leftActiveNode+1][j][k+1] - vectorZ[leftActiveNode+1][j][k])*invdz + .25*(vectorZ[leftActiveNode][j+1][k+1] - vectorZ[leftActiveNode][j+1][k])*invdz + .25*(vectorZ[leftActiveNode+1][j+1][k+1] - vectorZ[leftActiveNode+1][j+1][k])*invdz;
          divBC[leftActiveNode][j][k] = compX + compY + compZ;
      }
     break;
     case 1:    //  DIVERGENCE DIRECTION  Y
     for (register int i=1;i < nxn-2;i++)
        for (register int k=1;k < nzn-2;k++){
          compX = .25*(vectorX[i+1][leftActiveNode][k] - vectorX[i][leftActiveNode][k])*invdx +  .25*(vectorX[i+1][leftActiveNode][k+1] - vectorX[i][leftActiveNode][k+1])*invdx + .25*(vectorX[i+1][leftActiveNode+1][k] - vectorX[i][leftActiveNode+1][k])*invdx + .25*(vectorX[i+1][leftActiveNode+1][k+1] - vectorX[i][leftActiveNode+1][k+1])*invdx;
          compY = .25*(vectorY[i][leftActiveNode+1][k] - vectorY[i][leftActiveNode][k])*invdy +  .25*(vectorY[i][leftActiveNode+1][k+1] - vectorY[i][leftActiveNode][k+1])*invdy + .25*(vectorY[i+1][leftActiveNode+1][k] - vectorY[i+1][leftActiveNode][k])*invdy + .25*(vectorY[i+1][leftActiveNode+1][k+1] - vectorY[i+1][leftActiveNode][k+1])*invdy;
          compZ = .25*(vectorZ[i][leftActiveNode][k+1] - vectorZ[i][leftActiveNode][k])*invdz +  .25*(vectorZ[i+1][leftActiveNode][k+1] - vectorZ[i+1][leftActiveNode][k])*invdz + .25*(vectorZ[i][leftActiveNode+1][k+1] - vectorZ[i][leftActiveNode+1][k])*invdz + .25*(vectorZ[i+1][leftActiveNode+1][k+1] - vectorZ[i+1][leftActiveNode+1][k])*invdz;
          divBC[i][leftActiveNode][k] = compX + compY + compZ;
     }  
     break;
     case 2:   // DIVERGENCE DIRECTION  Z
     for (register int i=1;i < nxn-2;i++)
      for (register int j=1;j < nyn-2;j++){
          compX = .25*(vectorX[i+1][j][leftActiveNode] - vectorX[i][j][leftActiveNode])*invdx +  .25*(vectorX[i+1][j][leftActiveNode+1] - vectorX[i][j][leftActiveNode+1])*invdx + .25*(vectorX[i+1][j+1][leftActiveNode] - vectorX[i][j+1][leftActiveNode])*invdx + .25*(vectorX[i+1][j+1][leftActiveNode+1] - vectorX[i][j+1][leftActiveNode+1])*invdx;
          compY = .25*(vectorY[i][j+1][leftActiveNode] - vectorY[i][j][leftActiveNode])*invdy +  .25*(vectorY[i][j+1][leftActiveNode+1] - vectorY[i][j][leftActiveNode+1])*invdy + .25*(vectorY[i+1][j+1][leftActiveNode] - vectorY[i+1][j][leftActiveNode])*invdy + .25*(vectorY[i+1][j+1][leftActiveNode+1] - vectorY[i+1][j][leftActiveNode+1])*invdy;
          compZ = .25*(vectorZ[i][j][leftActiveNode+1] - vectorZ[i][j][leftActiveNode])*invdz +  .25*(vectorZ[i+1][j][leftActiveNode+1] - vectorZ[i+1][j][leftActiveNode])*invdz + .25*(vectorZ[i][j+1][leftActiveNode+1] - vectorZ[i][j+1][leftActiveNode])*invdz + .25*(vectorZ[i+1][j+1][leftActiveNode+1] - vectorZ[i+1][j+1][leftActiveNode])*invdz;
          divBC[i][j][leftActiveNode] = compX + compY + compZ;
     } 
     break;

    }
   

}
/** calculate divergence on  boundaries */
inline void Grid3DCU::divBCright(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ,int rightActiveNode, int dirDER){
   double compX, compY, compZ;
  

    switch(dirDER){
     case 0:    //  DIVERGENCE DIRECTION  X
      for (register int j=1;j < nxn-2;j++)
        for (register int k=1;k < nzn-2;k++){
          compX = .25*(vectorX[rightActiveNode][j][k] - vectorX[rightActiveNode-1][j][k])*invdx +  .25*(vectorX[rightActiveNode][j][k+1] - vectorX[rightActiveNode-1][j][k+1])*invdx + .25*(vectorX[rightActiveNode][j+1][k] - vectorX[rightActiveNode-1][j+1][k])*invdx + .25*(vectorX[rightActiveNode][j+1][k+1] - vectorX[rightActiveNode-1][j+1][k+1])*invdx;
          compY = .25*(vectorY[rightActiveNode][j+1][k] - vectorY[rightActiveNode][j][k])*invdy +  .25*(vectorY[rightActiveNode][j+1][k+1] - vectorY[rightActiveNode][j][k+1])*invdy + .25*(vectorY[rightActiveNode-1][j+1][k] - vectorY[rightActiveNode-1][j][k])*invdy + .25*(vectorY[rightActiveNode-1][j+1][k+1] - vectorY[rightActiveNode-1][j][k+1])*invdy;
          compZ = .25*(vectorZ[rightActiveNode][j][k+1] - vectorZ[rightActiveNode][j][k])*invdz +  .25*(vectorZ[rightActiveNode-1][j][k+1] - vectorZ[rightActiveNode-1][j][k])*invdz + .25*(vectorZ[rightActiveNode][j+1][k+1] - vectorZ[rightActiveNode][j+1][k])*invdz + .25*(vectorZ[rightActiveNode-1][j+1][k+1] - vectorZ[rightActiveNode-1][j+1][k])*invdz;
          divBC[rightActiveNode][j][k] = compX + compY + compZ;
      }
     break;
     case 1:    //  DIVERGENCE DIRECTION  Y
     for (register int i=1;i < nxn-2;i++)
        for (register int k=1;k < nzn-2;k++){
          compX = .25*(vectorX[i+1][rightActiveNode][k] - vectorX[i][rightActiveNode][k])*invdx +  .25*(vectorX[i+1][rightActiveNode][k+1] - vectorX[i][rightActiveNode][k+1])*invdx + .25*(vectorX[i+1][rightActiveNode-1][k] - vectorX[i][rightActiveNode-1][k])*invdx + .25*(vectorX[i+1][rightActiveNode-1][k+1] - vectorX[i][rightActiveNode-1][k+1])*invdx;
          compY = .25*(vectorY[i][rightActiveNode][k] - vectorY[i][rightActiveNode-1][k])*invdy +  .25*(vectorY[i][rightActiveNode][k+1] - vectorY[i][rightActiveNode-1][k+1])*invdy + .25*(vectorY[i+1][rightActiveNode][k] - vectorY[i+1][rightActiveNode-1][k])*invdy + .25*(vectorY[i+1][rightActiveNode+1][k+1] - vectorY[i+1][rightActiveNode][k+1])*invdy;
          compZ = .25*(vectorZ[i][rightActiveNode][k+1] - vectorZ[i][rightActiveNode][k])*invdz +  .25*(vectorZ[i+1][rightActiveNode][k+1] - vectorZ[i+1][rightActiveNode][k])*invdz + .25*(vectorZ[i][rightActiveNode-1][k+1] - vectorZ[i][rightActiveNode-1][k])*invdz + .25*(vectorZ[i+1][rightActiveNode-1][k+1] - vectorZ[i+1][rightActiveNode-1][k])*invdz;
          divBC[i][rightActiveNode][k] = compX + compY + compZ;
     }
     break;
     case 2:   // DIVERGENCE DIRECTION  Z
     for (register int i=1;i < nxn-2;i++)
      for (register int j=1;j < nyn-2;j++){
          compX = .25*(vectorX[i+1][j][rightActiveNode] - vectorX[i][j][rightActiveNode])*invdx +  .25*(vectorX[i+1][j][rightActiveNode-1] - vectorX[i][j][rightActiveNode-1])*invdx + .25*(vectorX[i+1][j+1][rightActiveNode] - vectorX[i][j+1][rightActiveNode])*invdx + .25*(vectorX[i+1][j+1][rightActiveNode-1] - vectorX[i][j+1][rightActiveNode-1])*invdx;
          compY = .25*(vectorY[i][j+1][rightActiveNode] - vectorY[i][j][rightActiveNode])*invdy +  .25*(vectorY[i][j+1][rightActiveNode-1] - vectorY[i][j][rightActiveNode-1])*invdy + .25*(vectorY[i+1][j+1][rightActiveNode] - vectorY[i+1][j][rightActiveNode])*invdy + .25*(vectorY[i+1][j+1][rightActiveNode-1] - vectorY[i+1][j][rightActiveNode-1])*invdy;
          compZ = .25*(vectorZ[i][j][rightActiveNode] - vectorZ[i][j][rightActiveNode-1])*invdz +  .25*(vectorZ[i+1][j][rightActiveNode] - vectorZ[i+1][j][rightActiveNode-1])*invdz + .25*(vectorZ[i][j+1][rightActiveNode] - vectorZ[i][j+1][rightActiveNode-1])*invdz + .25*(vectorZ[i+1][j+1][rightActiveNode] - vectorZ[i+1][j+1][rightActiveNode-1])*invdz;
          divBC[i][j][rightActiveNode] = compX + compY + compZ;
     }
     break;

    }
  

}
/** calculate derivative on left boundary */
inline void Grid3DCU::derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER){
    switch(dirDER){
     case 0:    //  DERIVATIVE DIRECTION  X
      for (register int j=1;j < nyc-1;j++)
        for (register int k=1;k < nzc-1;k++)
          derBC[leftActiveNode][j][k] = .25*(vector[leftActiveNode+1][j][k] - vector[leftActiveNode][j][k])*invdx +  .25*(vector[leftActiveNode+1][j][k+1] - vector[leftActiveNode][j][k+1])*invdx + .25*(vector[leftActiveNode+1][j+1][k] - vector[leftActiveNode][j+1][k])*invdx + .25*(vector[leftActiveNode+1][j+1][k+1] - vector[leftActiveNode][j+1][k+1])*invdx;;
      
     break;
     case 1:    //  DIVERGENCE DIRECTION  Y
     for (register int i=1;i < nxc-1;i++)
        for (register int k=1;k < nzc-1;k++)
          derBC[i][leftActiveNode][k] = .25*(vector[i][leftActiveNode+1][k] - vector[i][leftActiveNode][k])*invdy +  .25*(vector[i][leftActiveNode+1][k+1] - vector[i][leftActiveNode][k+1])*invdy + .25*(vector[i+1][leftActiveNode+1][k] - vector[i+1][leftActiveNode][k])*invdy + .25*(vector[i+1][leftActiveNode+1][k+1] - vector[i+1][leftActiveNode][k+1])*invdy;
     break;
     case 2:   // DIVERGENCE DIRECTION  Z
     for (register int i=1;i < nxc-1;i++)
      for (register int j=1;j < nyc-1;j++)
          derBC[i][j][leftActiveNode] = .25*(vector[i][j][leftActiveNode+1] - vector[i][j][leftActiveNode])*invdz +  .25*(vector[i+1][j][leftActiveNode+1] - vector[i+1][j][leftActiveNode])*invdz + .25*(vector[i][j+1][leftActiveNode+1] - vector[i][j+1][leftActiveNode])*invdz + .25*(vector[i+1][j+1][leftActiveNode+1] - vector[i+1][j+1][leftActiveNode])*invdz;
     break;

    }
}
/** interpolate on nodes from central points: do this for the magnetic field*/
inline void Grid3DCU::interpC2N(double ***vecFieldN, double ***vecFieldC){
	for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
        for (register int k=1;k < nzn-1;k++) 
            vecFieldN[i][j][k] = .125*(vecFieldC[i][j][k] + vecFieldC[i-1][j][k] + vecFieldC[i][j-1][k] + vecFieldC[i][j][k-1] + vecFieldC[i-1][j-1][k] + vecFieldC[i-1][j][k-1] + vecFieldC[i][j-1][k-1] + vecFieldC[i-1][j-1][k-1]);
}
/** interpolate on central points from nodes */
inline void Grid3DCU::interpN2C(double ***vecFieldC, double ***vecFieldN){
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++)
          vecFieldC[i][j][k] = .125*(vecFieldN[i][j][k] + vecFieldN[i+1][j][k] + vecFieldN[i][j+1][k] + vecFieldN[i][j][k+1] + vecFieldN[i+1][j+1][k] + vecFieldN[i+1][j][k+1] + vecFieldN[i][j+1][k+1] + vecFieldN[i+1][j+1][k+1]);
}
/** interpolate on central points from nodes */
inline void Grid3DCU::interpN2C(double ****vecFieldC, int ns, double ****vecFieldN){
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
       for (register int k=1;k < nzc-1;k++)
          vecFieldC[ns][i][j][k] = .125*(vecFieldN[ns][i][j][k] + vecFieldN[ns][i+1][j][k] + vecFieldN[ns][i][j+1][k] + vecFieldN[ns][i][j][k+1] + vecFieldN[ns][i+1][j+1][k] + vecFieldN[ns][i+1][j][k+1] + vecFieldN[ns][i][j+1][k+1] + vecFieldN[ns][i+1][j+1][k+1]);
 }


/** get nxc */
inline int Grid3DCU::getNXC(){
  return(nxc);
}
/** get nxn */
inline int Grid3DCU::getNXN(){
  return(nxn);
}
/** get nyc */
inline int Grid3DCU::getNYC(){
  return(nyc);
}
/** get nyn */
inline int Grid3DCU::getNYN(){
  return(nyn);
}
/** get nzc */
inline int Grid3DCU::getNZC(){
  return(nzc);
}
/** get nzn */
inline int Grid3DCU::getNZN(){
  return(nzn);
}
/** get dx */
inline double Grid3DCU::getDX(){
  return(dx);
}
/** get dy */
inline double Grid3DCU::getDY(){
  return(dy);
}
/** get dz */
inline double Grid3DCU::getDZ(){
  return(dz);
}

/** get xn[][][] */
inline double &Grid3DCU::getXN(int indexX, int indexY,int indexZ){
  return(node_coordinate[indexX][indexY][indexZ][0]);
}
/** get yn[][][] */
inline double &Grid3DCU::getYN(int indexX, int indexY,int indexZ){
  return(node_coordinate[indexX][indexY][indexZ][1]);
}
/** get zn[][][] */
inline double &Grid3DCU::getZN(int indexX, int indexY,int indexZ){
  return(node_coordinate[indexX][indexY][indexZ][2]);
}

/** get xc[][][] */
inline double &Grid3DCU::getXC(int indexX, int indexY,int indexZ){
  return(center_coordinate[indexX][indexY][indexZ][0]);
}
/** get yc[][][] */
inline double &Grid3DCU::getYC(int indexX, int indexY,int indexZ){
  return(center_coordinate[indexX][indexY][indexZ][1]);
}
/** get zc[][][] */
inline double &Grid3DCU::getZC(int indexX, int indexY,int indexZ){
  return(center_coordinate[indexX][indexY][indexZ][2]);
}

/** get the whole vector xc*/
inline double*** Grid3DCU::getXC(){
  return(center_coordinate[0]);
}
/** get the whole vector yc*/
inline double*** Grid3DCU::getYC(){
  return(center_coordinate[1]);
}
/** get the whole vector zc*/
inline double*** Grid3DCU::getZC(){
  return(center_coordinate[2]);
}
/** get Xstart */
inline double Grid3DCU::getXstart(){
 return(xStart);
}
/** get Xend */
inline double Grid3DCU::getXend(){
 return(xEnd);
}
/** get Ystart */
inline double Grid3DCU::getYstart(){
 return(yStart);
}
/** get Yend */
inline double Grid3DCU::getYend(){
 return(yEnd);
}
/** get Zstart */
inline double Grid3DCU::getZstart(){
 return(zStart);
}
/** get Zend */
inline double Grid3DCU::getZend(){
 return(zEnd);
}
/** get the inverse of volume */
inline double Grid3DCU::getInvVOL(){
 return(invVOL);
}
#endif
