
#include "Grid3DCU.h"

/*! constructor */
Grid3DCU::Grid3DCU(Collective * col, VirtualTopology3D * vct) {
  // FOR TESTS - this must be uncommented next
  // int get_rank();
  // if(!get_rank())
  // {
  // fflush(stdout);
  // bool xerror = false;
  // bool yerror = false;
  // bool zerror = false;
  // if((col->getNxc()) % (vct->getXLEN())) xerror=true;
  // if((col->getNyc()) % (vct->getYLEN())) yerror=true;
  // if((col->getNzc()) % (vct->getZLEN())) zerror=true;
  // if(xerror) printf("!!!ERROR: XLEN=%d does not divide nxc=%d\n", vct->getXLEN(),col->getNxc());
  // if(yerror) printf("!!!ERROR: YLEN=%d does not divide nyc=%d\n", vct->getYLEN(),col->getNyc());
  // if(zerror) printf("!!!ERROR: ZLEN=%d does not divide nzc=%d\n", vct->getZLEN(),col->getNzc());
  // fflush(stdout);
  // bool error = xerror||yerror||zerror;
  // if(error) exit(1);
  // }

  /*! mlmd: know your grid number */
  numGrid = vct->getNumGrid();

  /*! mlmd: know where you start on your PARENT grid */
  Ox= col->getOx_P(numGrid); Oy= col->getOy_P(numGrid); Oz= col->getOz_P(numGrid);

  /*! mlmd: know where you start in the system */
  Ox_SW= col->getOx_SW(numGrid); Oy_SW= col->getOy_SW(numGrid); Oz_SW= col->getOz_SW(numGrid);

  // add 2 for the guard cells
  nxc = (col->getNxc_mlmd(numGrid)) / (vct->getXLEN()) + 2;
  nyc = (col->getNyc_mlmd(numGrid)) / (vct->getYLEN()) + 2;
  nzc = (col->getNzc_mlmd(numGrid)) / (vct->getZLEN()) + 2;
  
  nxn = nxc + 1;
  nyn = nyc + 1;
  nzn = nzc + 1;

  dx = col->getDx_mlmd(numGrid);
  dy = col->getDy_mlmd(numGrid);
  dz = col->getDz_mlmd(numGrid);
  
  invVOL = 1.0 / (dx * dy * dz);
  invdx = 1.0 / dx;
  invdy = 1.0 / dy;
  invdz = 1.0 / dz;

  // local grid dimensions and boundaries of active nodes
  xStart = vct->getCoordinates(0) * (col->getLx_mlmd(numGrid) / (double) vct->getXLEN());
  xEnd = xStart + (col->getLx_mlmd(numGrid) / (double) vct->getXLEN());
  yStart = vct->getCoordinates(1) * (col->getLy_mlmd(numGrid) / (double) vct->getYLEN());
  yEnd = yStart + (col->getLy_mlmd(numGrid) / (double) vct->getYLEN());
  zStart = vct->getCoordinates(2) * (col->getLz_mlmd(numGrid) / (double) vct->getZLEN());
  zEnd = zStart + (col->getLz_mlmd(numGrid) / (double) vct->getZLEN());

  // if @ boudaries with MPI_PROC_NULL, local grid dimensions INCLUDING ghost areas; otherwise, 'usual' values
  // since it is used only by particles, use particle periodicity
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL)
    xStart_GC= xStart-dx;
  else
    xStart_GC= xStart;

  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL)
    yStart_GC= yStart-dy;
  else
    yStart_GC= yStart;

  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL)
    zStart_GC= zStart-dz;
  else
    zStart_GC= zStart;

  if (vct->getXright_neighbor_P() == MPI_PROC_NULL)
    xEnd_GC= xEnd+dx;
  else
    xEnd_GC= xEnd;

  if (vct->getYright_neighbor_P() == MPI_PROC_NULL)
    yEnd_GC= yEnd+dy;
  else
    yEnd_GC= yEnd;

  if (vct->getZright_neighbor_P() == MPI_PROC_NULL)
    zEnd_GC= zEnd+dz;
  else
    zEnd_GC= zEnd;


  /* portion of ACTIVE grid hosted in each parent core                                                                                       
     -- equivalent of xEnd - xStart on the parent */
  if (vct->getCommToParent()!= MPI_COMM_NULL) { // I have a parent               
    int PG= vct->getParentGridNum();
    parentLenX= col->getLx_mlmd(PG) / double(col->getXLEN_mlmd(PG)); 
    parentLenY= col->getLy_mlmd(PG) / double(col->getYLEN_mlmd(PG));
    parentLenZ= col->getLz_mlmd(PG) / double(col->getZLEN_mlmd(PG));
  }
  else{parentLenX=-1.; parentLenY=-1.; parentLenZ=-1.;} // just to set a value

  /* total number of grids */
  Ngrids= col->getNgrids();
  /* dx, dy, dz of all grids */
  dx_mlmd= new double [Ngrids];
  dy_mlmd= new double [Ngrids];
  dz_mlmd= new double [Ngrids];

  for (int g=0; g< Ngrids; g++){
    dx_mlmd[g]= col->getDx_mlmd(g);
    dy_mlmd[g]= col->getDy_mlmd(g);
    dz_mlmd[g]= col->getDz_mlmd(g);
  }

  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if (numGrid==0){
    cout << "I am grid " << numGrid <<endl;
    cout <<"xStart= " << xStart <<", xEnd= " <<xEnd <<", zStart= " << zStart <<", zEnd= " <<zEnd <<endl;
    cout << "Lx= " << col->getLx_mlmd(numGrid) << ", Ly= " <<col->getLy_mlmd(numGrid) <<", Lz= " << col->getLz_mlmd(numGrid) <<endl;
    cout <<"dx= " <<dx <<", dy=" <<dy <<", dz=" <<dz <<endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);*/

  // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
  node_xcoord = new double[nxn];
  node_ycoord = new double[nyn];
  node_zcoord = new double[nzn];
  for (int i=0; i<nxn; i++) node_xcoord[i] = xStart + (i - 1) * dx;
  for (int j=0; j<nyn; j++) node_ycoord[j] = yStart + (j - 1) * dy;
  for (int k=0; k<nzn; k++) node_zcoord[k] = zStart + (k - 1) * dz;
  // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
  center_xcoord = new double[nxc];
  center_ycoord = new double[nyc];
  center_zcoord = new double[nzc];
  for(int i=0; i<nxc; i++) center_xcoord[i] = .5*(node_xcoord[i]+node_xcoord[i+1]);
  for(int j=0; j<nyc; j++) center_ycoord[j] = .5*(node_ycoord[j]+node_ycoord[j+1]);
  for(int k=0; k<nzc; k++) center_zcoord[k] = .5*(node_zcoord[k]+node_zcoord[k+1]);

  // this type used in EMfields3D.cpp and Particles3Dcomm.cpp
  //MPI_RGBC_struct_commit();
  
}

/** deallocate the local grid */
Grid3DCU::~Grid3DCU() {
  delete [] node_xcoord;
  delete [] node_ycoord;
  delete [] node_zcoord;
  delete [] center_xcoord;
  delete [] center_ycoord;
  delete [] center_zcoord;

  delete [] dx_mlmd;
  delete [] dy_mlmd;
  delete [] dz_mlmd;
}

/** print the local grid info */
void Grid3DCU::print(VirtualTopology3D * ptVCT) {
  cout << endl;
  cout << "Grid number: " << numGrid << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << "Number of cell: -X=" << nxc - 2 << " -Y=" << nyc - 2 << " -Z=" << nzc - 2 << endl;
  cout << "Xin = " << node_xcoord[1] << "; Xfin = " << node_xcoord[nxn - 2] << endl;
  cout << "Yin = " << node_ycoord[1] << "; Yfin = " << node_ycoord[nyn - 2] << endl;
  cout << "Zin = " << node_zcoord[1] << "; Zfin = " << node_zcoord[nzn - 2] << endl;
  cout << endl;
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
void Grid3DCU::gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC) {
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        gradXN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i - 1][j][k]) * invdx + .25 * (scFieldC[i][j][k - 1] - scFieldC[i - 1][j][k - 1]) * invdx + .25 * (scFieldC[i][j - 1][k] - scFieldC[i - 1][j - 1][k]) * invdx + .25 * (scFieldC[i][j - 1][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdx;
        gradYN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j - 1][k]) * invdy + .25 * (scFieldC[i][j][k - 1] - scFieldC[i][j - 1][k - 1]) * invdy + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j - 1][k]) * invdy + .25 * (scFieldC[i - 1][j][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdy;
        gradZN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j][k - 1]) * invdz + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j][k - 1]) * invdz + .25 * (scFieldC[i][j - 1][k] - scFieldC[i][j - 1][k - 1]) * invdz + .25 * (scFieldC[i - 1][j - 1][k] - scFieldC[i - 1][j - 1][k - 1]) * invdz;
      }
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
void Grid3DCU::gradC2N(double ***gradXN, double ***gradYN, double ***gradZN, double ***scFieldC, int nxn, int nyn, int nzn) {
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        gradXN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i - 1][j][k]) * invdx + .25 * (scFieldC[i][j][k - 1] - scFieldC[i - 1][j][k - 1]) * invdx + .25 * (scFieldC[i][j - 1][k] - scFieldC[i - 1][j - 1][k]) * invdx + .25 * (scFieldC[i][j - 1][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdx;
        gradYN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j - 1][k]) * invdy + .25 * (scFieldC[i][j][k - 1] - scFieldC[i][j - 1][k - 1]) * invdy + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j - 1][k]) * invdy + .25 * (scFieldC[i - 1][j][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdy;
        gradZN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j][k - 1]) * invdz + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j][k - 1]) * invdz + .25 * (scFieldC[i][j - 1][k] - scFieldC[i][j - 1][k - 1]) * invdz + .25 * (scFieldC[i - 1][j - 1][k] - scFieldC[i - 1][j - 1][k - 1]) * invdz;
      }
}

/** calculate gradient on nodes, given a scalar field defined on central points  -
    2D, for Poisson face */
void Grid3DCU::gradC2N_XSide(double **gradXN, double **gradYN, double **gradZN, double **scFieldC) {

    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        gradXN[j][k] = 0.0;
        gradYN[j][k] = .5 * (scFieldC[j][k] - scFieldC[j - 1][k]) * invdy + .5 * (scFieldC[j][k - 1] - scFieldC[j - 1][k - 1]) * invdy;
	gradZN[j][k] = .25 * (scFieldC[j][k] - scFieldC[j][k - 1]) * invdz + .25 * (scFieldC[j - 1][k] - scFieldC[j - 1][k - 1]) * invdz ;

      }
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
void Grid3DCU::gradN2C(double ***gradXC, double ***gradYC, double ***gradZC, double ***scFieldN) {
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        gradXC[i][j][k] = .25 * (scFieldN[i + 1][j][k] - scFieldN[i][j][k]) * invdx + .25 * (scFieldN[i + 1][j][k + 1] - scFieldN[i][j][k + 1]) * invdx + .25 * (scFieldN[i + 1][j + 1][k] - scFieldN[i][j + 1][k]) * invdx + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i][j + 1][k + 1]) * invdx;
        gradYC[i][j][k] = .25 * (scFieldN[i][j + 1][k] - scFieldN[i][j][k]) * invdy + .25 * (scFieldN[i][j + 1][k + 1] - scFieldN[i][j][k + 1]) * invdy + .25 * (scFieldN[i + 1][j + 1][k] - scFieldN[i + 1][j][k]) * invdy + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i + 1][j][k + 1]) * invdy;
        gradZC[i][j][k] = .25 * (scFieldN[i][j][k + 1] - scFieldN[i][j][k]) * invdz + .25 * (scFieldN[i + 1][j][k + 1] - scFieldN[i + 1][j][k]) * invdz + .25 * (scFieldN[i][j + 1][k + 1] - scFieldN[i][j + 1][k]) * invdz + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i + 1][j + 1][k]) * invdz;
      }
}

/** calculate divergence on central points, given a vector field defined on nodes  */
void Grid3DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN) {
  double compX;
  double compY;
  double compZ;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        compX = .25 * (vecFieldXN[i + 1][j][k] - vecFieldXN[i][j][k]) * invdx + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i][j][k + 1]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i][j + 1][k]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i][j + 1][k + 1]) * invdx;
        compY = .25 * (vecFieldYN[i][j + 1][k] - vecFieldYN[i][j][k]) * invdy + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j][k + 1]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i + 1][j][k]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j][k + 1]) * invdy;
        compZ = .25 * (vecFieldZN[i][j][k + 1] - vecFieldZN[i][j][k]) * invdz + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i + 1][j][k]) * invdz + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j + 1][k]) * invdz + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j + 1][k]) * invdz;
        divC[i][j][k] = compX + compY + compZ;
      }
}

/** calculate divergence on central points, given a vector field defined on nodes  -
    on vectors with side redefined **/
void Grid3DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN, int nxc, int nyc, int nzc) {
  double compX;
  double compY;
  double compZ;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        compX = .25 * (vecFieldXN[i + 1][j][k] - vecFieldXN[i][j][k]) * invdx + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i][j][k + 1]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i][j + 1][k]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i][j + 1][k + 1]) * invdx;
        compY = .25 * (vecFieldYN[i][j + 1][k] - vecFieldYN[i][j][k]) * invdy + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j][k + 1]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i + 1][j][k]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j][k + 1]) * invdy;
        compZ = .25 * (vecFieldZN[i][j][k + 1] - vecFieldZN[i][j][k]) * invdz + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i + 1][j][k]) * invdz + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j + 1][k]) * invdz + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j + 1][k]) * invdz;
        divC[i][j][k] = compX + compY + compZ;
      }
}
/** calculate divergence on central points, given a vector field defined on nodes                     
    - for the face with CELL coordinate= I  */
void Grid3DCU::divN2C_XSide(double **divC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN, int I){ 
  double compX;
  double compY;
  double compZ;
  
  int i=I;

  for (register int j = 1; j < nyc - 1; j++)
    for (register int k = 1; k < nzc - 1; k++) {
      compX = .25 * (vecFieldXN[i + 1][j][k] - vecFieldXN[i][j][k]) * invdx + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i][j][k + 1]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i][j + 1][k]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i][j + 1][k + 1]) * invdx;
      compY = .25 * (vecFieldYN[i][j + 1][k] - vecFieldYN[i][j][k]) * invdy + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j][k + 1]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i + 1][j][k]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j][k + 1]) * invdy;
      compZ = .25 * (vecFieldZN[i][j][k + 1] - vecFieldZN[i][j][k]) * invdz + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i + 1][j][k]) * invdz + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j + 1][k]) * invdz + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j + 1][k]) * invdz;
      //divC[i][j][k] = compX + compY + compZ;
      
      divC[j][k] = compX + compY + compZ; 
    }
}

/** calculate divergence on central points, given a Tensor field defined on nodes  */
void Grid3DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns) {
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
  double comp1Z, comp2Z, comp3Z;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        comp1X = .25 * (pXX[ns][i + 1][j][k] - pXX[ns][i][j][k]) * invdx + .25 * (pXX[ns][i + 1][j][k + 1] - pXX[ns][i][j][k + 1]) * invdx + .25 * (pXX[ns][i + 1][j + 1][k] - pXX[ns][i][j + 1][k]) * invdx + .25 * (pXX[ns][i + 1][j + 1][k + 1] - pXX[ns][i][j + 1][k + 1]) * invdx;
        comp2X = .25 * (pXY[ns][i + 1][j][k] - pXY[ns][i][j][k]) * invdx + .25 * (pXY[ns][i + 1][j][k + 1] - pXY[ns][i][j][k + 1]) * invdx + .25 * (pXY[ns][i + 1][j + 1][k] - pXY[ns][i][j + 1][k]) * invdx + .25 * (pXY[ns][i + 1][j + 1][k + 1] - pXY[ns][i][j + 1][k + 1]) * invdx;
        comp3X = .25 * (pXZ[ns][i + 1][j][k] - pXZ[ns][i][j][k]) * invdx + .25 * (pXZ[ns][i + 1][j][k + 1] - pXZ[ns][i][j][k + 1]) * invdx + .25 * (pXZ[ns][i + 1][j + 1][k] - pXZ[ns][i][j + 1][k]) * invdx + .25 * (pXZ[ns][i + 1][j + 1][k + 1] - pXZ[ns][i][j + 1][k + 1]) * invdx;
        comp1Y = .25 * (pXY[ns][i][j + 1][k] - pXY[ns][i][j][k]) * invdy + .25 * (pXY[ns][i][j + 1][k + 1] - pXY[ns][i][j][k + 1]) * invdy + .25 * (pXY[ns][i + 1][j + 1][k] - pXY[ns][i + 1][j][k]) * invdy + .25 * (pXY[ns][i + 1][j + 1][k + 1] - pXY[ns][i + 1][j][k + 1]) * invdy;
        comp2Y = .25 * (pYY[ns][i][j + 1][k] - pYY[ns][i][j][k]) * invdy + .25 * (pYY[ns][i][j + 1][k + 1] - pYY[ns][i][j][k + 1]) * invdy + .25 * (pYY[ns][i + 1][j + 1][k] - pYY[ns][i + 1][j][k]) * invdy + .25 * (pYY[ns][i + 1][j + 1][k + 1] - pYY[ns][i + 1][j][k + 1]) * invdy;
        comp3Y = .25 * (pYZ[ns][i][j + 1][k] - pYZ[ns][i][j][k]) * invdy + .25 * (pYZ[ns][i][j + 1][k + 1] - pYZ[ns][i][j][k + 1]) * invdy + .25 * (pYZ[ns][i + 1][j + 1][k] - pYZ[ns][i + 1][j][k]) * invdy + .25 * (pYZ[ns][i + 1][j + 1][k + 1] - pYZ[ns][i + 1][j][k + 1]) * invdy;
        comp1Z = .25 * (pXZ[ns][i][j][k + 1] - pXZ[ns][i][j][k]) * invdz + .25 * (pXZ[ns][i + 1][j][k + 1] - pXZ[ns][i + 1][j][k]) * invdz + .25 * (pXZ[ns][i][j + 1][k + 1] - pXZ[ns][i][j + 1][k]) * invdz + .25 * (pXZ[ns][i + 1][j + 1][k + 1] - pXZ[ns][i + 1][j + 1][k]) * invdz;
        comp2Z = .25 * (pYZ[ns][i][j][k + 1] - pYZ[ns][i][j][k]) * invdz + .25 * (pYZ[ns][i + 1][j][k + 1] - pYZ[ns][i + 1][j][k]) * invdz + .25 * (pYZ[ns][i][j + 1][k + 1] - pYZ[ns][i][j + 1][k]) * invdz + .25 * (pYZ[ns][i + 1][j + 1][k + 1] - pYZ[ns][i + 1][j + 1][k]) * invdz;
        comp3Z = .25 * (pZZ[ns][i][j][k + 1] - pZZ[ns][i][j][k]) * invdz + .25 * (pZZ[ns][i + 1][j][k + 1] - pZZ[ns][i + 1][j][k]) * invdz + .25 * (pZZ[ns][i][j + 1][k + 1] - pZZ[ns][i][j + 1][k]) * invdz + .25 * (pZZ[ns][i + 1][j + 1][k + 1] - pZZ[ns][i + 1][j + 1][k]) * invdz;
        
	divCX[i][j][k] = comp1X + comp1Y + comp1Z;
	divCY[i][j][k] = comp2X + comp2Y + comp2Z;
	divCZ[i][j][k] = comp3X + comp3Y + comp3Z;
      }
}

/** calculate divergence on nodes, given a vector field defined on central points  */
void Grid3DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC) {
  double compX;
  double compY;
  double compZ;
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        compX = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i - 1][j][k]) * invdx + .25 * (vecFieldXC[i][j][k - 1] - vecFieldXC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldXC[i][j - 1][k] - vecFieldXC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldXC[i][j - 1][k - 1] - vecFieldXC[i - 1][j - 1][k - 1]) * invdx;
        compY = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i][j - 1][k]) * invdy + .25 * (vecFieldYC[i][j][k - 1] - vecFieldYC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldYC[i - 1][j][k] - vecFieldYC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldYC[i - 1][j][k - 1] - vecFieldYC[i - 1][j - 1][k - 1]) * invdy;
        compZ = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i][j][k - 1]) * invdz + .25 * (vecFieldZC[i - 1][j][k] - vecFieldZC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldZC[i][j - 1][k] - vecFieldZC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldZC[i - 1][j - 1][k] - vecFieldZC[i - 1][j - 1][k - 1]) * invdz;
        divN[i][j][k] = compX + compY + compZ;
      }
}

/** calculate curl on nodes, given a vector field defined on central points  */
void Grid3DCU::curlC2N(double ***curlXN, double ***curlYN, double ***curlZN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC) {
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        // curl - X
        compZDY = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i][j - 1][k]) * invdy + .25 * (vecFieldZC[i][j][k - 1] - vecFieldZC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldZC[i - 1][j][k] - vecFieldZC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldZC[i - 1][j][k - 1] - vecFieldZC[i - 1][j - 1][k - 1]) * invdy;
        compYDZ = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i][j][k - 1]) * invdz + .25 * (vecFieldYC[i - 1][j][k] - vecFieldYC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldYC[i][j - 1][k] - vecFieldYC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldYC[i - 1][j - 1][k] - vecFieldYC[i - 1][j - 1][k - 1]) * invdz;
        // curl - Y
        compXDZ = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i][j][k - 1]) * invdz + .25 * (vecFieldXC[i - 1][j][k] - vecFieldXC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldXC[i][j - 1][k] - vecFieldXC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldXC[i - 1][j - 1][k] - vecFieldXC[i - 1][j - 1][k - 1]) * invdz;
        compZDX = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i - 1][j][k]) * invdx + .25 * (vecFieldZC[i][j][k - 1] - vecFieldZC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldZC[i][j - 1][k] - vecFieldZC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldZC[i][j - 1][k - 1] - vecFieldZC[i - 1][j - 1][k - 1]) * invdx;
        // curl - Z
        compYDX = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i - 1][j][k]) * invdx + .25 * (vecFieldYC[i][j][k - 1] - vecFieldYC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldYC[i][j - 1][k] - vecFieldYC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldYC[i][j - 1][k - 1] - vecFieldYC[i - 1][j - 1][k - 1]) * invdx;
        compXDY = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i][j - 1][k]) * invdy + .25 * (vecFieldXC[i][j][k - 1] - vecFieldXC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldXC[i - 1][j][k] - vecFieldXC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldXC[i - 1][j][k - 1] - vecFieldXC[i - 1][j - 1][k - 1]) * invdy;

        curlXN[i][j][k] = compZDY - compYDZ;
        curlYN[i][j][k] = compXDZ - compZDX;
        curlZN[i][j][k] = compYDX - compXDY;
      }
}

/** calculate curl on central points, given a vector field defined on nodes  */
void Grid3DCU::curlN2C(double ***curlXC, double ***curlYC, double ***curlZC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN) {
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        // curl - X
        compZDY = .25 * (vecFieldZN[i][j + 1][k] - vecFieldZN[i][j][k]) * invdy + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j][k + 1]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i + 1][j][k]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j][k + 1]) * invdy;
        compYDZ = .25 * (vecFieldYN[i][j][k + 1] - vecFieldYN[i][j][k]) * invdz + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i + 1][j][k]) * invdz + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j + 1][k]) * invdz + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j + 1][k]) * invdz;
        // curl - Y
        compXDZ = .25 * (vecFieldXN[i][j][k + 1] - vecFieldXN[i][j][k]) * invdz + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i + 1][j][k]) * invdz + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j + 1][k]) * invdz + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j + 1][k]) * invdz;
        compZDX = .25 * (vecFieldZN[i + 1][j][k] - vecFieldZN[i][j][k]) * invdx + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i][j][k + 1]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i][j + 1][k]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i][j + 1][k + 1]) * invdx;
        // curl - Z
        compYDX = .25 * (vecFieldYN[i + 1][j][k] - vecFieldYN[i][j][k]) * invdx + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i][j][k + 1]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i][j + 1][k]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i][j + 1][k + 1]) * invdx;
        compXDY = .25 * (vecFieldXN[i][j + 1][k] - vecFieldXN[i][j][k]) * invdy + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j][k + 1]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i + 1][j][k]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j][k + 1]) * invdy;


        curlXC[i][j][k] = compZDY - compYDZ;
        curlYC[i][j][k] = compXDZ - compZDX;
        curlZC[i][j][k] = compYDX - compXDY;
      }

}

/** calculate curl on central points, given a vector field defined on nodes                                                                                     
    calculate ghost cells also using ghost node info if at buondary **/
void Grid3DCU::curlN2C_Ghost(VirtualTopology3D * vct, double ***curlXC, double ***curlYC, double ***curlZC, double ***vecFieldXN, double ***vecFieldYN, double ***vecFieldZN){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;

  // I need to do all of them!!! (because i cannot do a communicate node afterwards)
  // I am adding a communicate node before, just to be safe; NOT a communicateNode BC
  
  // very imp: DO NOT DO ANY SORT OF COMMUNICATE HERE

  for (register int i = 0; i < nxc; i++)
    for (register int j = 0; j < nyc; j++)
      for (register int k = 0; k < nzc; k++) {
        // curl - X
        compZDY = .25 * (vecFieldZN[i][j + 1][k] - vecFieldZN[i][j][k]) * invdy + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j][k + 1]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i + 1][j][k]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j][k + 1]) * invdy;
        compYDZ = .25 * (vecFieldYN[i][j][k + 1] - vecFieldYN[i][j][k]) * invdz + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i + 1][j][k]) * invdz + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j + 1][k]) * invdz + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j + 1][k]) * invdz;
        // curl - Y
        compXDZ = .25 * (vecFieldXN[i][j][k + 1] - vecFieldXN[i][j][k]) * invdz + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i + 1][j][k]) * invdz + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j + 1][k]) * invdz + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j + 1][k]) * invdz;
        compZDX = .25 * (vecFieldZN[i + 1][j][k] - vecFieldZN[i][j][k]) * invdx + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i][j][k + 1]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i][j + 1][k]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i][j + 1][k + 1]) * invdx;
        // curl - Z
        compYDX = .25 * (vecFieldYN[i + 1][j][k] - vecFieldYN[i][j][k]) * invdx + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i][j][k + 1]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i][j + 1][k]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i][j + 1][k + 1]) * invdx;
        compXDY = .25 * (vecFieldXN[i][j + 1][k] - vecFieldXN[i][j][k]) * invdy + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j][k + 1]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i + 1][j][k]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j][k + 1]) * invdy;


        curlXC[i][j][k] = compZDY - compYDZ;
        curlYC[i][j][k] = compXDZ - compZDX;
        curlZC[i][j][k] = compYDX - compXDY;
      }



}


/** calculate laplacian on nodes, given a scalar field defined on nodes */
void Grid3DCU::lapN2N(double ***lapN, double ***scFieldN, VirtualTopology3D * vct) {
  // calculate laplacian as divercence of gradient
  // allocate 3 gradients: defined on central points
  double ***gradXC = newArr3(double, nxc, nyc, nzc);
  double ***gradYC = newArr3(double, nxc, nyc, nzc);
  double ***gradZC = newArr3(double, nxc, nyc, nzc);

  // ghost not updated
  gradN2C(gradXC, gradYC, gradZC, scFieldN);
  // communicate with BC - 1: ghost set to 0
  communicateCenterBC(nxc, nyc, nzc, gradXC, 1, 1, 1, 1, 1, 1, vct);
  communicateCenterBC(nxc, nyc, nzc, gradYC, 1, 1, 1, 1, 1, 1, vct);
  communicateCenterBC(nxc, nyc, nzc, gradZC, 1, 1, 1, 1, 1, 1, vct);
  // ghost updated
  divC2N(lapN, gradXC, gradYC, gradZC);
  // deallocate
  delArr3(gradXC, nxc, nyc);
  delArr3(gradYC, nxc, nyc);
  delArr3(gradZC, nxc, nyc);
}

void Grid3DCU::lapN2N_mlmd(double ***lapN, double ***scFieldN, VirtualTopology3D * vct) {
  // calculate laplacian as divercence of gradient
  // allocate 3 gradients: defined on central points
  double ***gradXC = newArr3(double, nxc, nyc, nzc);
  double ***gradYC = newArr3(double, nxc, nyc, nzc);
  double ***gradZC = newArr3(double, nxc, nyc, nzc);

  // if you put ghost node, ghost cell is ok
  gradN2C(gradXC, gradYC, gradZC, scFieldN);
  // communicate with BC
  communicateCenter(nxc, nyc, nzc, gradXC, vct);
  communicateCenter(nxc, nyc, nzc, gradYC, vct);
  communicateCenter(nxc, nyc, nzc, gradZC, vct);
  divC2N(lapN, gradXC, gradYC, gradZC);
  // deallocate
  delArr3(gradXC, nxc, nyc);
  delArr3(gradYC, nxc, nyc);
  delArr3(gradZC, nxc, nyc);
}

/** calculate laplacian on central points, given a scalar field defined on central points */
void Grid3DCU::lapC2C(double ***lapC, double ***scFieldC, VirtualTopology3D * vct) {
  // calculate laplacian as divercence of gradient
  // allocate 3 gradients: defined on nodes
  double ***gradXN = newArr3(double, nxn, nyn, nzn);
  double ***gradYN = newArr3(double, nxn, nyn, nzn);
  double ***gradZN = newArr3(double, nxn, nyn, nzn);

  gradC2N(gradXN, gradYN, gradZN, scFieldC);
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int ii = 0; ii < nxn; ii++)
      for (int kk = 0; kk < nzn; kk++) {
        gradXN[ii][0][kk] = 0.0;
        gradXN[ii][1][kk] = 0.0;
        gradXN[ii][2][kk] = 0.0;
        gradZN[ii][0][kk] = 0.0;
        gradZN[ii][1][kk] = 0.0;
        gradZN[ii][2][kk] = 0.0;
        // gradYN[ii][1][kk] = gradYN[ii][2][kk];
      }
  }
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int ii = 0; ii < nxn; ii++)
      for (int kk = 0; kk < nzn; kk++) {
        gradXN[ii][nyn - 1][kk] = 0.0;
        gradXN[ii][nyn - 2][kk] = 0.0;
        gradXN[ii][nyn - 3][kk] = 0.0;
        gradZN[ii][nyn - 1][kk] = 0.0;
        gradZN[ii][nyn - 2][kk] = 0.0;
        gradZN[ii][nyn - 3][kk] = 0.0;
        // gradYN[ii][nyn-2][kk] = gradYN[ii][nyc-3][kk];
      }
  }
  divN2C(lapC, gradXN, gradYN, gradZN);

  delArr3(gradXN, nxn, nyn);
  delArr3(gradYN, nxn, nyn);
  delArr3(gradZN, nxn, nyn);

}

/** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
void Grid3DCU::lapC2Cpoisson(double ***lapC, double ***scFieldC, VirtualTopology3D * vct) {
  // communicate first the scFieldC
  communicateCenterBoxStencilBC(nxc, nyc, nzc, scFieldC, 1, 1, 1, 1, 1, 1, vct);
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        lapC[i][j][k] = (scFieldC[i - 1][j][k] - 2 * scFieldC[i][j][k] + scFieldC[i + 1][j][k]) * invdx * invdx + (scFieldC[i][j - 1][k] - 2 * scFieldC[i][j][k] + scFieldC[i][j + 1][k]) * invdy * invdy + (scFieldC[i][j][k - 1] - 2 * scFieldC[i][j][k] + scFieldC[i][j][k + 1]) * invdz * invdz;
}

void Grid3DCU::lapC2Cpoisson(double ***lapC, double ***scFieldC, VirtualTopology3D * vct, int nxc, int nyc, int nzc) {
  // communicate first the scFieldC
  //communicateCenterBoxStencilBC(nxc, nyc, nzc, scFieldC, 1, 1, 1, 1, 1, 1, vct);
  cout << "I have to do a communication in lapC2Cpoisson, now nothing there" << endl;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        lapC[i][j][k] = (scFieldC[i - 1][j][k] - 2 * scFieldC[i][j][k] + scFieldC[i + 1][j][k]) * invdx * invdx + (scFieldC[i][j - 1][k] - 2 * scFieldC[i][j][k] + scFieldC[i][j + 1][k]) * invdy * invdy + (scFieldC[i][j][k - 1] - 2 * scFieldC[i][j][k] + scFieldC[i][j][k + 1]) * invdz * invdz;


      }
}

/** calculate divergence on  boundaries */
void Grid3DCU::divBCleft(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER) {
  double compX, compY, compZ;
  switch (dirDER) {
    case 0:                    // DIVERGENCE DIRECTION X
      for (register int j = 1; j < nyn - 2; j++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[leftActiveNode + 1][j][k] - vectorX[leftActiveNode][j][k]) * invdx + .25 * (vectorX[leftActiveNode + 1][j][k + 1] - vectorX[leftActiveNode][j][k + 1]) * invdx + .25 * (vectorX[leftActiveNode + 1][j + 1][k] - vectorX[leftActiveNode][j + 1][k]) * invdx + .25 * (vectorX[leftActiveNode + 1][j + 1][k + 1] - vectorX[leftActiveNode][j + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[leftActiveNode][j + 1][k] - vectorY[leftActiveNode][j][k]) * invdy + .25 * (vectorY[leftActiveNode][j + 1][k + 1] - vectorY[leftActiveNode][j][k + 1]) * invdy + .25 * (vectorY[leftActiveNode + 1][j + 1][k] - vectorY[leftActiveNode + 1][j][k]) * invdy + .25 * (vectorY[leftActiveNode + 1][j + 1][k + 1] - vectorY[leftActiveNode + 1][j][k + 1]) * invdy;
          compZ = .25 * (vectorZ[leftActiveNode][j][k + 1] - vectorZ[leftActiveNode][j][k]) * invdz + .25 * (vectorZ[leftActiveNode + 1][j][k + 1] - vectorZ[leftActiveNode + 1][j][k]) * invdz + .25 * (vectorZ[leftActiveNode][j + 1][k + 1] - vectorZ[leftActiveNode][j + 1][k]) * invdz + .25 * (vectorZ[leftActiveNode + 1][j + 1][k + 1] - vectorZ[leftActiveNode + 1][j + 1][k]) * invdz;
          divBC[leftActiveNode][j][k] = compX + compY + compZ;
        }
      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxn - 2; i++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[i + 1][leftActiveNode][k] - vectorX[i][leftActiveNode][k]) * invdx + .25 * (vectorX[i + 1][leftActiveNode][k + 1] - vectorX[i][leftActiveNode][k + 1]) * invdx + .25 * (vectorX[i + 1][leftActiveNode + 1][k] - vectorX[i][leftActiveNode + 1][k]) * invdx + .25 * (vectorX[i + 1][leftActiveNode + 1][k + 1] - vectorX[i][leftActiveNode + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[i][leftActiveNode + 1][k] - vectorY[i][leftActiveNode][k]) * invdy + .25 * (vectorY[i][leftActiveNode + 1][k + 1] - vectorY[i][leftActiveNode][k + 1]) * invdy + .25 * (vectorY[i + 1][leftActiveNode + 1][k] - vectorY[i + 1][leftActiveNode][k]) * invdy + .25 * (vectorY[i + 1][leftActiveNode + 1][k + 1] - vectorY[i + 1][leftActiveNode][k + 1]) * invdy;
          compZ = .25 * (vectorZ[i][leftActiveNode][k + 1] - vectorZ[i][leftActiveNode][k]) * invdz + .25 * (vectorZ[i + 1][leftActiveNode][k + 1] - vectorZ[i + 1][leftActiveNode][k]) * invdz + .25 * (vectorZ[i][leftActiveNode + 1][k + 1] - vectorZ[i][leftActiveNode + 1][k]) * invdz + .25 * (vectorZ[i + 1][leftActiveNode + 1][k + 1] - vectorZ[i + 1][leftActiveNode + 1][k]) * invdz;
          divBC[i][leftActiveNode][k] = compX + compY + compZ;
        }
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxn - 2; i++)
        for (register int j = 1; j < nyn - 2; j++) {
          compX = .25 * (vectorX[i + 1][j][leftActiveNode] - vectorX[i][j][leftActiveNode]) * invdx + .25 * (vectorX[i + 1][j][leftActiveNode + 1] - vectorX[i][j][leftActiveNode + 1]) * invdx + .25 * (vectorX[i + 1][j + 1][leftActiveNode] - vectorX[i][j + 1][leftActiveNode]) * invdx + .25 * (vectorX[i + 1][j + 1][leftActiveNode + 1] - vectorX[i][j + 1][leftActiveNode + 1]) * invdx;
          compY = .25 * (vectorY[i][j + 1][leftActiveNode] - vectorY[i][j][leftActiveNode]) * invdy + .25 * (vectorY[i][j + 1][leftActiveNode + 1] - vectorY[i][j][leftActiveNode + 1]) * invdy + .25 * (vectorY[i + 1][j + 1][leftActiveNode] - vectorY[i + 1][j][leftActiveNode]) * invdy + .25 * (vectorY[i + 1][j + 1][leftActiveNode + 1] - vectorY[i + 1][j][leftActiveNode + 1]) * invdy;
          compZ = .25 * (vectorZ[i][j][leftActiveNode + 1] - vectorZ[i][j][leftActiveNode]) * invdz + .25 * (vectorZ[i + 1][j][leftActiveNode + 1] - vectorZ[i + 1][j][leftActiveNode]) * invdz + .25 * (vectorZ[i][j + 1][leftActiveNode + 1] - vectorZ[i][j + 1][leftActiveNode]) * invdz + .25 * (vectorZ[i + 1][j + 1][leftActiveNode + 1] - vectorZ[i + 1][j + 1][leftActiveNode]) * invdz;
          divBC[i][j][leftActiveNode] = compX + compY + compZ;
        }
      break;

  }


}

/** calculate divergence on  boundaries */
void Grid3DCU::divBCright(double ***divBC, double ***vectorX, double ***vectorY, double ***vectorZ, int rightActiveNode, int dirDER) {
  double compX, compY, compZ;


  switch (dirDER) {
    case 0:                    // DIVERGENCE DIRECTION X
      for (register int j = 1; j < nxn - 2; j++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[rightActiveNode][j][k] - vectorX[rightActiveNode - 1][j][k]) * invdx + .25 * (vectorX[rightActiveNode][j][k + 1] - vectorX[rightActiveNode - 1][j][k + 1]) * invdx + .25 * (vectorX[rightActiveNode][j + 1][k] - vectorX[rightActiveNode - 1][j + 1][k]) * invdx + .25 * (vectorX[rightActiveNode][j + 1][k + 1] - vectorX[rightActiveNode - 1][j + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[rightActiveNode][j + 1][k] - vectorY[rightActiveNode][j][k]) * invdy + .25 * (vectorY[rightActiveNode][j + 1][k + 1] - vectorY[rightActiveNode][j][k + 1]) * invdy + .25 * (vectorY[rightActiveNode - 1][j + 1][k] - vectorY[rightActiveNode - 1][j][k]) * invdy + .25 * (vectorY[rightActiveNode - 1][j + 1][k + 1] - vectorY[rightActiveNode - 1][j][k + 1]) * invdy;
          compZ = .25 * (vectorZ[rightActiveNode][j][k + 1] - vectorZ[rightActiveNode][j][k]) * invdz + .25 * (vectorZ[rightActiveNode - 1][j][k + 1] - vectorZ[rightActiveNode - 1][j][k]) * invdz + .25 * (vectorZ[rightActiveNode][j + 1][k + 1] - vectorZ[rightActiveNode][j + 1][k]) * invdz + .25 * (vectorZ[rightActiveNode - 1][j + 1][k + 1] - vectorZ[rightActiveNode - 1][j + 1][k]) * invdz;
          divBC[rightActiveNode][j][k] = compX + compY + compZ;
        }
      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxn - 2; i++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[i + 1][rightActiveNode][k] - vectorX[i][rightActiveNode][k]) * invdx + .25 * (vectorX[i + 1][rightActiveNode][k + 1] - vectorX[i][rightActiveNode][k + 1]) * invdx + .25 * (vectorX[i + 1][rightActiveNode - 1][k] - vectorX[i][rightActiveNode - 1][k]) * invdx + .25 * (vectorX[i + 1][rightActiveNode - 1][k + 1] - vectorX[i][rightActiveNode - 1][k + 1]) * invdx;
          compY = .25 * (vectorY[i][rightActiveNode][k] - vectorY[i][rightActiveNode - 1][k]) * invdy + .25 * (vectorY[i][rightActiveNode][k + 1] - vectorY[i][rightActiveNode - 1][k + 1]) * invdy + .25 * (vectorY[i + 1][rightActiveNode][k] - vectorY[i + 1][rightActiveNode - 1][k]) * invdy + .25 * (vectorY[i + 1][rightActiveNode + 1][k + 1] - vectorY[i + 1][rightActiveNode][k + 1]) * invdy;
          compZ = .25 * (vectorZ[i][rightActiveNode][k + 1] - vectorZ[i][rightActiveNode][k]) * invdz + .25 * (vectorZ[i + 1][rightActiveNode][k + 1] - vectorZ[i + 1][rightActiveNode][k]) * invdz + .25 * (vectorZ[i][rightActiveNode - 1][k + 1] - vectorZ[i][rightActiveNode - 1][k]) * invdz + .25 * (vectorZ[i + 1][rightActiveNode - 1][k + 1] - vectorZ[i + 1][rightActiveNode - 1][k]) * invdz;
          divBC[i][rightActiveNode][k] = compX + compY + compZ;
        }
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxn - 2; i++)
        for (register int j = 1; j < nyn - 2; j++) {
          compX = .25 * (vectorX[i + 1][j][rightActiveNode] - vectorX[i][j][rightActiveNode]) * invdx + .25 * (vectorX[i + 1][j][rightActiveNode - 1] - vectorX[i][j][rightActiveNode - 1]) * invdx + .25 * (vectorX[i + 1][j + 1][rightActiveNode] - vectorX[i][j + 1][rightActiveNode]) * invdx + .25 * (vectorX[i + 1][j + 1][rightActiveNode - 1] - vectorX[i][j + 1][rightActiveNode - 1]) * invdx;
          compY = .25 * (vectorY[i][j + 1][rightActiveNode] - vectorY[i][j][rightActiveNode]) * invdy + .25 * (vectorY[i][j + 1][rightActiveNode - 1] - vectorY[i][j][rightActiveNode - 1]) * invdy + .25 * (vectorY[i + 1][j + 1][rightActiveNode] - vectorY[i + 1][j][rightActiveNode]) * invdy + .25 * (vectorY[i + 1][j + 1][rightActiveNode - 1] - vectorY[i + 1][j][rightActiveNode - 1]) * invdy;
          compZ = .25 * (vectorZ[i][j][rightActiveNode] - vectorZ[i][j][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i + 1][j][rightActiveNode] - vectorZ[i + 1][j][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i][j + 1][rightActiveNode] - vectorZ[i][j + 1][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i + 1][j + 1][rightActiveNode] - vectorZ[i + 1][j + 1][rightActiveNode - 1]) * invdz;
          divBC[i][j][rightActiveNode] = compX + compY + compZ;
        }
      break;

  }


}

/** calculate derivative on left boundary */
void Grid3DCU::derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER) {
  switch (dirDER) {
    case 0:                    // DERIVATIVE DIRECTION X
      for (register int j = 1; j < nyc - 1; j++)
        for (register int k = 1; k < nzc - 1; k++)
          derBC[leftActiveNode][j][k] = .25 * (vector[leftActiveNode + 1][j][k] - vector[leftActiveNode][j][k]) * invdx + .25 * (vector[leftActiveNode + 1][j][k + 1] - vector[leftActiveNode][j][k + 1]) * invdx + .25 * (vector[leftActiveNode + 1][j + 1][k] - vector[leftActiveNode][j + 1][k]) * invdx + .25 * (vector[leftActiveNode + 1][j + 1][k + 1] - vector[leftActiveNode][j + 1][k + 1]) * invdx;;

      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxc - 1; i++)
        for (register int k = 1; k < nzc - 1; k++)
          derBC[i][leftActiveNode][k] = .25 * (vector[i][leftActiveNode + 1][k] - vector[i][leftActiveNode][k]) * invdy + .25 * (vector[i][leftActiveNode + 1][k + 1] - vector[i][leftActiveNode][k + 1]) * invdy + .25 * (vector[i + 1][leftActiveNode + 1][k] - vector[i + 1][leftActiveNode][k]) * invdy + .25 * (vector[i + 1][leftActiveNode + 1][k + 1] - vector[i + 1][leftActiveNode][k + 1]) * invdy;
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxc - 1; i++)
        for (register int j = 1; j < nyc - 1; j++)
          derBC[i][j][leftActiveNode] = .25 * (vector[i][j][leftActiveNode + 1] - vector[i][j][leftActiveNode]) * invdz + .25 * (vector[i + 1][j][leftActiveNode + 1] - vector[i + 1][j][leftActiveNode]) * invdz + .25 * (vector[i][j + 1][leftActiveNode + 1] - vector[i][j + 1][leftActiveNode]) * invdz + .25 * (vector[i + 1][j + 1][leftActiveNode + 1] - vector[i + 1][j + 1][leftActiveNode]) * invdz;
      break;

  }
}

/** interpolate on nodes from central points: do this for the magnetic field*/
void Grid3DCU::interpC2N(double ****vecFieldN, int ns, double ****vecFieldC) {
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++)
        vecFieldN[ns][i][j][k] = .125 * (vecFieldC[ns][i][j][k] + vecFieldC[ns][i - 1][j][k] + vecFieldC[ns][i][j - 1][k] + vecFieldC[ns][i][j][k - 1] + vecFieldC[ns][i - 1][j - 1][k] + vecFieldC[ns][i - 1][j][k - 1] + vecFieldC[ns][i][j - 1][k - 1] + vecFieldC[ns][i - 1][j - 1][k - 1]);
}
/** interpolate on nodes from central points: do this for the magnetic field*/
void Grid3DCU::interpC2N(double ***vecFieldN, double ***vecFieldC) {
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++)
        vecFieldN[i][j][k] = .125 * (vecFieldC[i][j][k] + vecFieldC[i - 1][j][k] + vecFieldC[i][j - 1][k] + vecFieldC[i][j][k - 1] + vecFieldC[i - 1][j - 1][k] + vecFieldC[i - 1][j][k - 1] + vecFieldC[i][j - 1][k - 1] + vecFieldC[i - 1][j - 1][k - 1]);
}

/** interpolate on central points from nodes */
void Grid3DCU::interpN2C(double ***vecFieldC, double ***vecFieldN) {
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
	
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);
}

void Grid3DCU::interpN2C_GC(double ***vecFieldC, double ***vecFieldN) {
  for (register int i = 0; i < nxc ; i++)
    for (register int j = 0; j < nyc ; j++)
      for (register int k = 0; k < nzc ; k++)
	
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);
}

void Grid3DCU::interpN2C_ActiveCell(double ***vecFieldC, double ***vecFieldN, VirtualTopology3D * vct) {

  if (vct->getXleft_neighbor() == MPI_PROC_NULL){

    for (register int i = 1; i < 2; i++)
      for (register int j = 1; j < nyc - 1; j++)
	for (register int k = 1; k < nzc - 1; k++)
	  vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);
  }

  if (vct->getXright_neighbor() == MPI_PROC_NULL){

    for (register int i = nxc-2; i < nxc-1; i++)
      for (register int j = 1; j < nyc - 1; j++)
	for (register int k = 1; k < nzc - 1; k++)
	  vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);
  }

  if (vct->getYleft_neighbor() == MPI_PROC_NULL){
    for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < 2; j++)
      for (register int k = 1; k < nzc - 1; k++)
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]); }

  if (vct->getYright_neighbor() == MPI_PROC_NULL){
    for (register int i = 1; i < nxc - 1; i++)
      for (register int j = nyc-2; j < nyc-1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]); }

  if (vct->getZleft_neighbor() == MPI_PROC_NULL){
    for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < 2; k++)
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);}

  if (vct->getZright_neighbor() == MPI_PROC_NULL){

    for (register int i = 1; i < nxc - 1; i++)
      for (register int j = 1; j < nyc - 1; j++)
	for (register int k = nzc-2; k < nzc - 1; k++)
	  vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);}

}


/** interpolate on central points from nodes */
void Grid3DCU::interpN2C(double ****vecFieldC, int ns, double ****vecFieldN) {
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        vecFieldC[ns][i][j][k] = .125 * (vecFieldN[ns][i][j][k] + vecFieldN[ns][i + 1][j][k] + vecFieldN[ns][i][j + 1][k] + vecFieldN[ns][i][j][k + 1] + vecFieldN[ns][i + 1][j + 1][k] + vecFieldN[ns][i + 1][j][k + 1] + vecFieldN[ns][i][j + 1][k + 1] + vecFieldN[ns][i + 1][j + 1][k + 1]);
}

void Grid3DCU::interpN2C_GC(double ****vecFieldC, int ns, double ****vecFieldN) {
  for (register int i = 0; i < nxc ; i++)
    for (register int j = 0; j < nyc ; j++)
      for (register int k = 0; k < nzc ; k++)
        vecFieldC[ns][i][j][k] = .125 * (vecFieldN[ns][i][j][k] + vecFieldN[ns][i + 1][j][k] + vecFieldN[ns][i][j + 1][k] + vecFieldN[ns][i][j][k + 1] + vecFieldN[ns][i + 1][j + 1][k] + vecFieldN[ns][i + 1][j][k + 1] + vecFieldN[ns][i][j + 1][k + 1] + vecFieldN[ns][i + 1][j + 1][k + 1]);
}


/** get nxc */
int Grid3DCU::getNXC() {
  return (nxc);
}

/** get nxn */
int Grid3DCU::getNXN() {
  return (nxn);
}

/** get nyc */
int Grid3DCU::getNYC() {
  return (nyc);
}

/** get nyn */
int Grid3DCU::getNYN() {
  return (nyn);
}

/** get nzc */
int Grid3DCU::getNZC() {
  return (nzc);
}

/** get nzn */
int Grid3DCU::getNZN() {
  return (nzn);
}

/** get dx */
double Grid3DCU::getDX() {
  return (dx);
}

/** get dy */
double Grid3DCU::getDY() {
  return (dy);
}

/** get dz */
double Grid3DCU::getDZ() {
  return (dz);
}

/** get Xstart */
double Grid3DCU::getXstart() {
  return (xStart);
}

/** get Xend */
double Grid3DCU::getXend() {
  return (xEnd);
}

/** get Ystart */
double Grid3DCU::getYstart() {
  return (yStart);
}

/** get Yend */
double Grid3DCU::getYend() {
  return (yEnd);
}

/** get Zstart */
double Grid3DCU::getZstart() {
  return (zStart);
}

/** get Zend */
double Grid3DCU::getZend() {
  return (zEnd);
}

/** get the inverse of volume */
double Grid3DCU::getInvVOL() {
  return (invVOL);
}

double Grid3DCU::getxStart_GC(){return xStart_GC;}
double Grid3DCU::getyStart_GC(){return yStart_GC;}
double Grid3DCU::getzStart_GC(){return zStart_GC;}
double Grid3DCU::getxEnd_GC(){return xEnd_GC;}
double Grid3DCU::getyEnd_GC(){return yEnd_GC;}
double Grid3DCU::getzEnd_GC(){return zEnd_GC;}

/** mlmd specific functions **/
int Grid3DCU::getParentRankFromGridPoint(VirtualTopology3D * vct, int xn, int yn, int zn){

    /*** pay exceptional attention to this description ***/
    /* nx, ny, nz: index in the current grid, which is a child*/
    /* returns the rank IN THE PARENT-CHILD communicator of the coarse grid core where the point is hosted                                   
       only the active part of the parent grid is examined*/

  // March 27: examining also the ghost cells of the parent grid

  // rank on MPI_COMM_WORLD
  int SW_rank=vct->getSystemWide_rank();

  if (vct->getCommToParent()== MPI_COMM_NULL){
    cout << "Grid " <<numGrid  <<" R" << SW_rank << ": Fatal error in getParentRankFromGridPoint:" << endl ;
    MPI_Abort(MPI_COMM_WORLD, -1);
    return -1; // if you are not a child, i return -1 to provoke a segm fault
  }
  // coordinates in the parent grid                                                                 
  double coordX_PG= getXN_P(xn, yn, zn);
  double coordY_PG= getYN_P(xn, yn, zn);
  double coordZ_PG= getZN_P(xn, yn, zn);

  // !cartesian! coordinates in the parent grid                                                                         
  int coordX= floor(coordX_PG/ parentLenX); //if (coordX_PG= ) QUI: check if it's the ghost
  int coordY= floor(coordY_PG/ parentLenY);
  int coordZ= floor(coordZ_PG/ parentLenZ);

  

  /* end the search when the data relative to the child grid itself start,                                                                    
     i.e. at index XLEN_mlmd[parentGrid]*YLEN_mlmd[parentGrid]*ZLEN_mlmd[parentGrid] */
  // rank in the CommToParent communicator                                                                                                    
  int rankPC=-1;
  int ChildStarts= vct->getXLEN(vct->getParentGridNum()) * vct->getYLEN(vct->getParentGridNum()) * vct->getZLEN(vct->getParentGridNum());
  //cout <<"R" << SW_rank  <<" parentGridNum: " << vct->getParentGridNum() << " vct->getXLEN(vct->getParentGridNum()): " << vct->getXLEN(vct->getParentGridNum()) << ",  vct->getYLEN(vct->getParentGridNum()): " <<  vct->getYLEN(vct->getParentGridNum()) << ", vct->getZLEN(vct->getParentGridNum()): " << vct->getZLEN(vct->getParentGridNum()) << endl;
  //cout <<"R" << SW_rank << " parentLenX " << parentLenX << " parentLenY " << parentLenY << " parentLenZ " << parentLenZ << endl; 
  for (int i=0; i< ChildStarts; i++){
    /*cout <<"R" << SW_rank << " i=" <<i << " vct->getXcoord_CommToParent(i): " << vct->getXcoord_CommToParent(i) << " coordX " << coordX << endl;
    cout <<"R" << SW_rank << " i=" <<i << " vct->getYcoord_CommToParent(i): " << vct->getYcoord_CommToParent(i) << " coordY " << coordY << endl;
    cout <<"R" << SW_rank << " i=" <<i << " vct->getZcoord_CommToParent(i): " << vct->getZcoord_CommToParent(i) << " coordZ " << coordZ << endl;*/
    if ((coordX == vct->getXcoord_CommToParent(i)) && (coordY == vct->getYcoord_CommToParent(i)) && (coordZ == vct->getZcoord_CommToParent(i)) ){
      rankPC= i;
    }
  }

  if (rankPC==-1){
    cout <<"R" << SW_rank  << " grid " << numGrid<< " Fatal error in getParentRankFromGridPoint";
    cout <<"R" << SW_rank  << " You need to change the reciprocal positions of the grids, aborting .." << endl;
    cout <<"R" << SW_rank  << " I was looking for the parent rank of point (GC coords) " << coordX_PG <<" - " << coordY_PG << " - " << coordZ_PG << endl;
    cout <<"R" << SW_rank  << " Calculated CG coords: " << coordX <<" - " << coordY << " - " << coordZ << endl; 
    cout << "Aborting now ..." << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  else {
    if (false){
      cout << "R" <<SW_rank <<": local coords: [ " << getXN(xn, yn, zn) <<", "<< getYN(xn, yn, zn)<<", "<<getZN(xn, yn,zn) <<endl;
      cout << "R" <<SW_rank <<": on parent grid: [ " << getXN_P(xn, yn, zn) <<", "<< getYN_P(xn, yn, zn)<<", "<<getZN_P(xn, yn, zn) << " (origin at [ " << Ox <<", " << Oy  <<", " << Oz <<"] )";
      cout << "R" <<SW_rank <<": hosted in parent grid core: " << rankPC <<endl;
    }
    return rankPC;
  }
}


/** fix3B **/
int Grid3DCU::getParentRankFromGridCenter(VirtualTopology3D * vct, int xn, int yn, int zn){

    /*** pay exceptional attention to this description ***/
    /* nx, ny, nz: index in the current grid, which is a child*/
    /* returns the rank IN THE PARENT-CHILD communicator of the coarse grid core where the center is hosted       
       only the active part of the parent grid is examined*/

  // rank on MPI_COMM_WORLD
  int SW_rank=vct->getSystemWide_rank();

  if (vct->getCommToParent()== MPI_COMM_NULL){
    cout << "Grid " <<numGrid  <<" R" << SW_rank << ": Fatal error in getParentRankFromGridPoint:" << endl ;
    MPI_Abort(MPI_COMM_WORLD, -1);
    return -1; // if you are not a child, i return -1 to provoke a segm fault
  }
  // coordinates in the parent grid --> here moved to centers
  double coordX_PG= getXC_P(xn, yn, zn);
  double coordY_PG= getYC_P(xn, yn, zn);
  double coordZ_PG= getZC_P(xn, yn, zn);

  // !cartesian! coordinates in the parent grid                                                                         
  int coordX= floor(coordX_PG/ parentLenX); //if (coordX_PG= ) QUI: check if it's the ghost
  int coordY= floor(coordY_PG/ parentLenY);
  int coordZ= floor(coordZ_PG/ parentLenZ);

  /* end the search when the data relative to the child grid itself start,                                                         
     i.e. at index XLEN_mlmd[parentGrid]*YLEN_mlmd[parentGrid]*ZLEN_mlmd[parentGrid] */
  // rank in the CommToParent communicator                                                                                          
  int rankPC=-1;
  int ChildStarts= vct->getXLEN(vct->getParentGridNum()) * vct->getYLEN(vct->getParentGridNum()) * vct->getZLEN(vct->getParentGridNum());

  for (int i=0; i< ChildStarts; i++){
    if ((coordX == vct->getXcoord_CommToParent(i)) && (coordY == vct->getYcoord_CommToParent(i)) && (coordZ == vct->getZcoord_CommToParent(i)) ){
      rankPC= i;
    }
  }

  if (rankPC==-1){
    cout <<"R" << SW_rank  << " grid " << numGrid<< " Fatal error in getParentRankFromGridCenter";
    cout <<"R" << SW_rank  << " You need to change the reciprocal positions of the grids, aborting .." << endl;
    cout <<"R" << SW_rank  << " I was looking for the parent rank of point (GC coords) " << coordX_PG <<" - " << coordY_PG << " - " << coordZ_PG << endl;
    cout <<"R" << SW_rank  << " Calculated CG coords: " << coordX <<" - " << coordY << " - " << coordZ << endl; 
    cout << "Aborting now ..." << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  else {
    if (false){
      cout << "R" <<SW_rank <<": local coords: [ " << getXC(xn, yn, zn) <<", "<< getYC(xn, yn, zn)<<", "<<getZC(xn, yn,zn) <<endl;
      cout << "R" <<SW_rank <<": on parent grid: [ " << getXC_P(xn, yn, zn) <<", "<< getYC_P(xn, yn, zn)<<", "<<getZC_P(xn, yn, zn) << " (origin at [ " << Ox <<", " << Oy  <<", " << Oz <<"] )";
      cout << "R" <<SW_rank <<": hosted in parent grid core: " << rankPC <<endl;
    }
    return rankPC;
  }
}

/** end fix3B **/

/** grid --> obvious     
    FACE --> (bottom, top, left, right, front, back) (needed only for debug prints)     
    DIR --> direction of the changing index (0->X, 1->Y, 2->Z)
    i0_s --> the index of the RG first point (included) in the direction to explore    
    i0_e --> the index of the RG last point (included) in the direction to explore 
    i1 --> fixed index, in the lower-order fixed direction (ordering X,Y,Z)     
    i2 --> fixed index, in the higher-order fixed direction (ordering X,Y,Z)     
    *SPXperC --> X coordinate (not index!) in the CG of the first point FOR EACH CG CORE     
    *SPYperC --> Y coordinate (not index!) in the CG of the first point FOR EACH CG CORE     
    *SPZperC --> Z coordinate (not index!) in the CG of the first point FOR EACH CG CORE     
    *NPperC --> # of point in this CG core per direction    
    *rank --> rank, IN THE PARENT-CHILD COMMUNICATOR, of the CG core     
    Ncores --> # of the CG cores involved in BC in this direction     
    *IndexFirstPointperC --> index in the RG of the first point per core in the selected direction      **/
void Grid3DCU::RGBCExploreDirection(VirtualTopology3D *vct,string FACE, int DIR, int i0_s, int i0_e, int i1, int i2, double *SPXperC, double *SPYperC, double *SPZperC, int *NPperC, int *rank, int* Ncores, int *IndexFirstPointperC){

  //cout << "Inside Explore: FACE ->" << FACE <<", DIR ->" << DIR << ", i0_s -->" << i0_s <<", i0_e -->" << i0_e <<", i1 --> " << i1 <<", i2 -->" << i2;              

  int SW_rank=vct->getSystemWide_rank();
  int rank_CommToParent;

  // NB: extremes are included; keep the <=
  if (DIR==0){ // changing index in X                                                      
    for(int i=i0_s; i<= i0_e; i++){
      int j=i1;
      int k=i2;
      //cout << "R" << SW_rank <<", " << FACE <<endl;                                    
      rank_CommToParent= getParentRankFromGridPoint(vct, i, j, k);
      if (i==i0_s){ // at the beginning                                                                       
	SPXperC[0]= getXN_P(i, j, k);
	SPYperC[0]= getYN_P(i, j, k);
	SPZperC[0]= getZN_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]= i0_s;
        *Ncores=1;
      } else { // not at the beginning                       
        if (rank_CommToParent== rank[*Ncores-1]){ // still in the same rank      
          NPperC[(*Ncores)-1]++;  // i only update the # of points               
	}else { // changed rank                                                                            
          SPXperC[*Ncores]= getXN_P(i, j, k);
          SPYperC[*Ncores]= getYN_P(i, j, k);
          SPZperC[*Ncores]= getZN_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=i;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning                                                             
    }// end for                                                                       
  } // end if dir 

  if (DIR==1){ //changing index in Y     
    for(int j=i0_s; j<= i0_e; j++){
      int i=i1;
      int k=i2;
      rank_CommToParent= getParentRankFromGridPoint(vct, i, j, k);
      if (j==i0_s){ // at the beginning                               
	SPXperC[0]= getXN_P(i, j, k);
        SPYperC[0]= getYN_P(i, j, k);
        SPZperC[0]= getZN_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]=i0_s;
        *Ncores=1;
      } else { // not at the beginning                                                               
        if (rank_CommToParent== rank[(*Ncores)-1]){ // still in the same rank     
          NPperC[(*Ncores)-1]++; // i only update the # of points                   
	}
        else { // changed rank           
	  SPXperC[*Ncores]= getXN_P(i, j, k);
          SPYperC[*Ncores]= getYN_P(i, j, k);
          SPZperC[*Ncores]= getZN_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=j;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning      
    }// end for                                                                     
  } // end if dir

  if (DIR==2){ //changing index in Z                                                  
    for(int k=i0_s; k<= i0_e; k++){
      int i=i1;
      int j=i2;
      rank_CommToParent= getParentRankFromGridPoint(vct, i, j, k);
      if (k==i0_s){ // at the beginning                                            
	SPXperC[0]= getXN_P(i, j, k);
        SPYperC[0]= getYN_P(i, j, k);
        SPZperC[0]= getZN_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]=i0_s;
        *Ncores=1;
      } else {
        if (rank_CommToParent== rank[*Ncores-1]){
          // still in the same rank                            
          NPperC[(*Ncores)-1]++; // i only update the # of points         
        }else { // changed rank                      
          SPXperC[*Ncores]= getXN_P(i, j, k);
          SPYperC[*Ncores]= getYN_P(i, j, k);
          SPZperC[*Ncores]= getZN_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=k;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning    

    }// end for                                         
  } // end if dir 

}

/** fix3B **/
void Grid3DCU::RGBCExploreDirection_Centers(VirtualTopology3D *vct,string FACE, int DIR, int i0_s, int i0_e, int i1, int i2, double *SPXperC, double *SPYperC, double *SPZperC, int *NPperC, int *rank, int* Ncores, int *IndexFirstPointperC){

  //cout << "Inside Explore: FACE ->" << FACE <<", DIR ->" << DIR << ", i0_s -->" << i0_s <<", i0_e -->" << i0_e <<", i1 --> " << i1 <<", i2 -->" << i2;              

  int SW_rank=vct->getSystemWide_rank();
  int rank_CommToParent;

  // NB: extremes are included; keep the <=
  if (DIR==0){ // changing index in X                                                      
    for(int i=i0_s; i<= i0_e; i++){
      int j=i1;
      int k=i2;
      //cout << "R" << SW_rank <<", " << FACE <<endl;                                    
      rank_CommToParent= getParentRankFromGridCenter(vct, i, j, k);
      if (i==i0_s){ // at the beginning                                                                       
	SPXperC[0]= getXC_P(i, j, k);
	SPYperC[0]= getYC_P(i, j, k);
	SPZperC[0]= getZC_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]= i0_s;
        *Ncores=1;
      } else { // not at the beginning                       
        if (rank_CommToParent== rank[*Ncores-1]){ // still in the same rank      
          NPperC[(*Ncores)-1]++;  // i only update the # of points               
	}else { // changed rank                                                                            
          SPXperC[*Ncores]= getXC_P(i, j, k);
          SPYperC[*Ncores]= getYC_P(i, j, k);
          SPZperC[*Ncores]= getZC_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=i;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning                                                             
    }// end for                                                                       
  } // end if dir 

  if (DIR==1){ //changing index in Y     
    for(int j=i0_s; j<= i0_e; j++){
      int i=i1;
      int k=i2;
      rank_CommToParent= getParentRankFromGridCenter(vct, i, j, k);
      if (j==i0_s){ // at the beginning                               
	SPXperC[0]= getXC_P(i, j, k);
        SPYperC[0]= getYC_P(i, j, k);
        SPZperC[0]= getZC_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]=i0_s;
        *Ncores=1;
      } else { // not at the beginning                                                               
        if (rank_CommToParent== rank[(*Ncores)-1]){ // still in the same rank     
          NPperC[(*Ncores)-1]++; // i only update the # of points                   
	}
        else { // changed rank           
	  SPXperC[*Ncores]= getXC_P(i, j, k);
          SPYperC[*Ncores]= getYC_P(i, j, k);
          SPZperC[*Ncores]= getZC_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=j;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning      
    }// end for                                                                     
  } // end if dir

  if (DIR==2){ //changing index in Z                                                  
    for(int k=i0_s; k<= i0_e; k++){
      int i=i1;
      int j=i2;
      rank_CommToParent= getParentRankFromGridCenter(vct, i, j, k);
      if (k==i0_s){ // at the beginning                                            
	SPXperC[0]= getXC_P(i, j, k);
        SPYperC[0]= getYC_P(i, j, k);
        SPZperC[0]= getZC_P(i, j, k);
        NPperC[0]= 1;
        rank[0]= rank_CommToParent;
        IndexFirstPointperC[0]=i0_s;
        *Ncores=1;
      } else {
        if (rank_CommToParent== rank[*Ncores-1]){
          // still in the same rank                            
          NPperC[(*Ncores)-1]++; // i only update the # of points         
        }else { // changed rank                      
          SPXperC[*Ncores]= getXC_P(i, j, k);
          SPYperC[*Ncores]= getYC_P(i, j, k);
          SPZperC[*Ncores]= getZC_P(i, j, k);
          NPperC[*Ncores]= 1;
          rank[*Ncores]= rank_CommToParent;
          IndexFirstPointperC[*Ncores]=k;
          (*Ncores)++;
        }// end changed rank 
      } // end else not at the beginning    

    }// end for                                         
  } // end if dir 

}

/** end fix3B **/


void Grid3DCU::getParentLimits(VirtualTopology3D *vct, int N, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax){
  
  int XC= vct->getXcoord_CommToParent(N);
  int YC= vct->getYcoord_CommToParent(N);
  int ZC= vct->getZcoord_CommToParent(N);

  *xmin= parentLenX* XC;
  *xmax= parentLenX* (XC+1);

  *ymin= parentLenY* YC;
  *ymax= parentLenY* (YC+1);

  *zmin= parentLenZ* ZC;
  *zmax= parentLenZ* (ZC+1);
}

/* mlmd: coordinate on your parent grid     
   NB: i want it to be able to manage also negative indexes or indexes > nxn/ nyn/ nzn    
   (for the phase 1 of particle init BC) */

double Grid3DCU::getXN_P(int X, int Y, int Z){ 
  double dx= node_xcoord[1]-node_xcoord[0];

  if (X>-1 and X<nxn) {
    return node_xcoord[X]+ Ox;} // "normal" case
  if (X<0){
    return Ox+ node_xcoord[0]+ dx*X;}
  if (X>nxn-1){
    return Ox+ node_xcoord[nxn-1] + dx*(X-(nxn-1));}

}

double Grid3DCU::getXC_P(int X, int Y, int Z){ 
  double dx= center_xcoord[1]-center_xcoord[0];

  if (X>-1 and X<nxc) {
    return center_xcoord[X]+ Ox;} // "normal" case
  if (X<0){
    return Ox+ center_xcoord[0]+ dx*X;}
  if (X>nxc-1){
    return Ox+ center_xcoord[nxc-1] + dx*(X-(nxc-1));}

}

double Grid3DCU::getYN_P(int X, int Y, int Z){ 
  double dy= node_ycoord[1]-node_ycoord[0];
  
  if (Y>-1 and Y<nyn) 
    return node_ycoord[Y]+ Oy; // "normal" case
  if (Y<0)
    return Oy+ node_ycoord[0]+ dy*Y;
  if (Y>nyn-1)
    return Oy+ node_ycoord[nyn-1] + dy*(Y-(nyn-1));
}

double Grid3DCU::getYC_P(int X, int Y, int Z){ 
  double dy= center_ycoord[1]-center_ycoord[0];
  
  if (Y>-1 and Y<nyc) 
    return center_ycoord[Y]+ Oy; // "normal" case
  if (Y<0)
    return Oy+ center_ycoord[0]+ dy*Y;
  if (Y>nyc-1)
    return Oy+ center_ycoord[nyc-1] + dy*(Y-(nyc-1));
}

double Grid3DCU::getZN_P(int X, int Y, int Z){ 
  double dz= node_zcoord[1]-node_zcoord[0];
  
  if (Z>-1 and Z<nzn) 
    return node_zcoord[Z]+ Oz; // "normal" case
  if (Z<0)
    return Oz+ node_zcoord[0]+ dz*Z;
  if (Z>nzn-1)
    return Oz+ node_zcoord[nzn-1] + dz*(Z-(nzn-1));
}

double Grid3DCU::getZC_P(int X, int Y, int Z){ 
  double dz= center_zcoord[1]-center_zcoord[0];
  
  if (Z>-1 and Z<nzc) 
    return center_zcoord[Z]+ Oz; // "normal" case
  if (Z<0)
    return Oz+ center_zcoord[0]+ dz*Z;
  if (Z>nzc-1)
    return Oz+ center_zcoord[nzn-1] + dz*(Z-(nzc-1));
}


double Grid3DCU::getXN_XT(int X, int Y, int Z){ 
  double dx= node_xcoord[1]-node_xcoord[0];

  if (X>-1 and X<nxn) {
    return node_xcoord[X];} // "normal" case
  if (X<0){
    return node_xcoord[0]+ dx*X;}
  if (X>nxn-1){
    return node_xcoord[nxn-1] + dx*(X-(nxn-1));}

}

double Grid3DCU::getYN_XT(int X, int Y, int Z){ 
  double dy= node_ycoord[1]-node_ycoord[0];
  
  if (Y>-1 and Y<nyn) 
    return node_ycoord[Y]; // "normal" case
  if (Y<0)
    return node_ycoord[0]+ dy*Y;
  if (Y>nyn-1)
    return node_ycoord[nyn-1] + dy*(Y-(nyn-1));
}

double Grid3DCU::getZN_XT(int X, int Y, int Z){ 
  double dz= node_zcoord[1]-node_zcoord[0];
  
  if (Z>-1 and Z<nzn) 
    return node_zcoord[Z]; // "normal" case
  if (Z<0)
    return node_zcoord[0]+ dz*Z;
  if (Z>nzn-1)
    return node_zcoord[nzn-1] + dz*(Z-(nzn-1));
}

void Grid3DCU::Explore3DAndCommit_Centers(int i_s, int i_e, int j_s, int j_e, int k_s, int k_e, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg, VirtualTopology3D * vct , char  dir){
  // policy:    
  // explore Z dir  
  // for on the number of cores found there: explore Y dir  
  // for on the number of cores found there: explore X dir   
  // finally, commit  NB: all faces should have the same c 
  int MS= nxc; if (nyc>MS) MS= nyc; if (nzc>MS) MS= nzc;
  int rank_CTP= vct->getRank_CommToParent();
  /*******************************************************************/
  // DIR1: starting point, in CG coordinates, per core 
  double *Dir1_SPXperC= new double[MS];
  double *Dir1_SPYperC= new double[MS];
  double *Dir1_SPZperC= new double[MS];
  // DIR1: number of Refined grid point in this direction, per core   
  int *Dir1_NPperC= new int[MS];  // this does not need to be this big, but anyway  
  // DIR1: core ranks in the CommToParent communicator 
  int *Dir1_rank= new int [MS]; // this does not need to be this big, but anyway   
  int *Dir1_IndexFirstPointperC= new int [MS];
  int Dir1_Ncores=0;

  // DIR2: starting point, in CG coordinates, per core 
  double *Dir2_SPXperC= new double[MS];
  double *Dir2_SPYperC= new double[MS];
  double *Dir2_SPZperC= new double[MS];
  // DIR2: number of Refined grid point in this direction, per core   
  int *Dir2_NPperC= new int[MS];  // this does not need to be this big, but anyway   
  // DIR2: core ranks in the CommToParent communicator 
  int *Dir2_rank= new int [MS]; // this does not need to be this big, but anyway   
  int *Dir2_IndexFirstPointperC= new int [MS];
  int Dir2_Ncores=0;

  // DIR3: starting point, in CG coordinates, per core 
  double *Dir3_SPXperC= new double[MS];
  double *Dir3_SPYperC= new double[MS];
  double *Dir3_SPZperC= new double[MS];
  // DIR3: number of Refined grid point in this direction, per core   
  int *Dir3_NPperC= new int[MS];  // this does not need to be this big, but anyway   
  // DIR3: core ranks in the CommToParent communicator
  int *Dir3_rank= new int [MS]; // this does not need to be this big, but anyway   
  int *Dir3_IndexFirstPointperC= new int [MS];
  int Dir3_Ncores=0;
  /*******************************************************************/

  string FACE="nn";
  // Z dir / Dir 1  
  RGBCExploreDirection_Centers(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

  for (int n=0; n<Dir1_Ncores; n++){ // it will find again the core in Dir 1, but it will also explore Dir 2    
    // Y dir / Dir 2     
      RGBCExploreDirection_Centers(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n],  Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC);

    for (int m=0; m< Dir2_Ncores; m++){ //it will find again the core in Dir 1, but it will also explore Dir 2
      // X dir / Dir 3
      RGBCExploreDirection_Centers(vct, FACE, 0, i_s, i_e, Dir2_IndexFirstPointperC[m],  Dir1_IndexFirstPointperC[n], Dir3_SPXperC, Dir3_SPYperC, Dir3_SPZperC, Dir3_NPperC, Dir3_rank, &Dir3_Ncores, Dir3_IndexFirstPointperC);

      for (int NN=0; NN< Dir3_Ncores; NN++){
        // using the function written for BCs, it's the same 
        Assign_RGBC_struct_Values(RGBC_Info + (*numMsg), Dir3_IndexFirstPointperC[NN], Dir2_IndexFirstPointperC[m], Dir1_IndexFirstPointperC[n], -1, Dir3_NPperC[NN], Dir2_NPperC[m], Dir1_NPperC[n], Dir3_SPXperC[NN], Dir3_SPYperC[NN], Dir3_SPZperC[NN], Dir3_rank[NN], rank_CTP, *numMsg);

        (*numMsg)++;

        int tmp= Dir3_NPperC[NN]*Dir2_NPperC[m]*Dir1_NPperC[n];
        if (tmp > (*MaxSizeMsg)) (*MaxSizeMsg)= tmp;

      } // end for (int NN=0; NN< Dir3_Ncores; NN++){  
    } // end  for (int m=0; m< Dir2_Ncores; m++){
  } // end for (int n=0; n<Dir1_Ncores; n++){                                            
  // end here, the same as Explore3DAndCommit for particles 


  // the msg signaling the end      
  RGBC_Info[*numMsg].RG_core= -1;
  RGBC_Info[*numMsg].CG_core= -1;
  // deletes                     
  delete[]Dir1_SPXperC;
  delete[]Dir1_SPYperC;
  delete[]Dir1_SPZperC;
  delete[]Dir1_NPperC;
  delete[]Dir1_rank;
  delete[]Dir1_IndexFirstPointperC;

  delete[]Dir2_SPXperC;
  delete[]Dir2_SPYperC;
  delete[]Dir2_SPZperC;
  delete[]Dir2_NPperC;
  delete[]Dir2_rank;
  delete[]Dir2_IndexFirstPointperC;

  delete[]Dir3_SPXperC;
  delete[]Dir3_SPYperC;
  delete[]Dir3_SPZperC;
  delete[]Dir3_NPperC;
  delete[]Dir3_rank;
  delete[]Dir3_IndexFirstPointperC;

}
// add one handshake msg to the list - to use in field initialisations, where you need to init more variables
void Assign_RGBC_struct_Values(RGBC_struct *s, int ix_first_tmp, int iy_first_tmp, int iz_first_tmp,int BCside_tmp, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp , int MsgID_tmp) {
  /// check the struct definition before assigning here!!!
  s->ix_first= ix_first_tmp;
  s->iy_first= iy_first_tmp;
  s->iz_first= iz_first_tmp;

  s->BCside= BCside_tmp;

  s->np_x= np_x_tmp;
  s->np_y= np_y_tmp;
  s->np_z= np_z_tmp;

  s->CG_x_first= CG_x_first_tmp;
  s->CG_y_first= CG_y_first_tmp;
  s->CG_z_first= CG_z_first_tmp;

  s->CG_core= CG_core_tmp;
  s->RG_core= RG_core_tmp;

  s->MsgID= MsgID_tmp;

  return;
}

// add one handshake msg to the list - to use in particle initialisations, where you need to init less variables
void Assign_RGBC_struct_Values(RGBC_struct *s, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp, int MsgID_tmp ) {

  s->np_x = np_x_tmp;
  s->np_y = np_y_tmp;
  s->np_z = np_z_tmp;

  s->CG_x_first = CG_x_first_tmp;
  s->CG_y_first = CG_y_first_tmp;
  s->CG_z_first = CG_z_first_tmp;

  s->CG_core = CG_core_tmp;
  
  s->RG_core = RG_core_tmp;

  s->MsgID = MsgID_tmp;

  return;
}


/*void Grid3DCU::MPI_RGBC_struct_commit(){

  RGBC_struct *a;
  MPI_Datatype type[13]={MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT};
  int blocklen[13]={1,1,1,1,1,1,1,1,1,1,1,1,1};
  // displacement in bytes                                                                                                                                                   
                                                                                                                                                                        
  MPI_Aint disp[13];

  disp[0]= (MPI_Aint) &(a->ix_first) - (MPI_Aint)a ;
  disp[1]= (MPI_Aint) &(a->iy_first) - (MPI_Aint)a ;
  disp[2]= (MPI_Aint) &(a->iz_first) - (MPI_Aint)a ;
  // BCside                                                                                                                                                                    
  disp[3]= (MPI_Aint) &(a->BCside) - (MPI_Aint)a ;
  // np_*                                                                                                                                                                      
  disp[4]= (MPI_Aint) &(a->np_x) - (MPI_Aint)a ;
  disp[5]= (MPI_Aint) &(a->np_y) - (MPI_Aint)a ;
  disp[6]= (MPI_Aint) &(a->np_z) - (MPI_Aint)a ;
  // CG_*_first                                                                                                                                                                
  disp[7]= (MPI_Aint) &(a->CG_x_first) - (MPI_Aint)a ;
  disp[8]= (MPI_Aint) &(a->CG_y_first) - (MPI_Aint)a ;
  disp[9]= (MPI_Aint) &(a->CG_z_first) - (MPI_Aint)a ;
  // the cores                                                                                                                                                                 
  disp[10]= (MPI_Aint) &(a->CG_core) - (MPI_Aint)a ;
  disp[11]= (MPI_Aint) &(a->RG_core) - (MPI_Aint)a ;
  // the msg id                                                                                                                                                                
  disp[12]= (MPI_Aint) &(a->MsgID) - (MPI_Aint)a ;

  MPI_Type_create_struct(13, blocklen, disp, type, &MPI_RGBC_struct);
  MPI_Type_commit(&MPI_RGBC_struct);

  }*/
/** end mlmd specific functions **/

