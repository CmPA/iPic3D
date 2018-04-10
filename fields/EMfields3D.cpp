#include "EMfields3D.h"
 
/*! constructor */
/*! pre-mlmd: no topology needed:
  mlmd: topology is needed for the grid communicator */
//EMfields3D::EMfields3D(Collective * col, Grid * grid) {

/*! TO DO: adapt BC to mlmd */

EMfields3D::EMfields3D(Collective * col, Grid * grid, VirtualTopology3D * vct) {

  /*! mlmd: 'grid->' values are already local to the mlmd grid */ 
  nxc = grid->getNXC();
  nxn = grid->getNXN();
  nyc = grid->getNYC();
  nyn = grid->getNYN();
  nzc = grid->getNZC();
  nzn = grid->getNZN();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  invVOL = grid->getInvVOL();
  xStart = grid->getXstart();
  xEnd = grid->getXend();
  yStart = grid->getYstart();
  yEnd = grid->getYend();
  zStart = grid->getZstart();
  zEnd = grid->getZend();
 
  /*! grid number in the mlmd grid hierarchy */
  numGrid= grid->getNumGrid(); 
  /* coordinates of origin on the parent grid */
  //Ox= grid->getOx(); Oy= grid->getOy(); Oz= grid->getOz();

  /*! number of children in the mlmd grid hierarchy */
  numChildren= vct->getNumChildren();
  dx_Ch=new double[numChildren];
  dy_Ch=new double[numChildren];
  dz_Ch=new double[numChildren];
  for (int i=0; i< numChildren; i++){
    c= vct->getChildGridNum(i);
    dx_Ch[i]= col->getDx_mlmd(c) ;
    dy_Ch[i]= col->getDy_mlmd(c);
    dz_Ch[i]= col->getDz_mlmd(c);
  }


  Lx = col->getLx_mlmd(numGrid);
  Ly = col->getLy_mlmd(numGrid);
  Lz = col->getLz_mlmd(numGrid);

  SmoothGrid = col->getSmoothGrid(numGrid);

  ns = col->getNs();
  c = col->getC();
  dt = col->getDt();
  th = col->getTh();
  ue0 = col->getU0(0);
  ve0 = col->getV0(0);
  we0 = col->getW0(0);
  /*! mlmd: adapt this to mlmd, at the moment not touched yet */
  x_center = col->getx_center(); 
  y_center = col->gety_center();
  z_center = col->getz_center();
  L_square = col->getL_square();
  /*! end mlmd: adapt this to mlmd, at the moment not touched yer */

  Fext = 0.0;

  delt = c * th * dt;
  PoissonCorrection = false;
  PoissonCorrection_RGFace = true;
  if (col->getPoissonCorrection()=="yes") PoissonCorrection = true;
 
  //if (numGrid >0) PoissonCorrection = true;

  CGtol = col->getCGtol();
  GMREStol = col->getGMREStol();
  qom = new double[ns];
  for (int i = 0; i < ns; i++)
    qom[i] = col->getQOM(i);
  // boundary conditions: PHI and EM fields
  bcPHIfaceXright = col->getBcPHIfaceXright();
  bcPHIfaceXleft = col->getBcPHIfaceXleft();
  bcPHIfaceYright = col->getBcPHIfaceYright();
  bcPHIfaceYleft = col->getBcPHIfaceYleft();
  bcPHIfaceZright = col->getBcPHIfaceZright();
  bcPHIfaceZleft = col->getBcPHIfaceZleft();

  bcEMfaceXright = col->getBcEMfaceXright();
  bcEMfaceXleft = col->getBcEMfaceXleft();
  bcEMfaceYright = col->getBcEMfaceYright();
  bcEMfaceYleft = col->getBcEMfaceYleft();
  bcEMfaceZright = col->getBcEMfaceZright();
  bcEMfaceZleft = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  // Earth Simulation
  B1x = col->getB1x();
  B1y = col->getB1y();
  B1z = col->getB1z();
  delta = col->getDelta();
  Smooth = col->getSmooth();
  // get the density background for the gem Challange
  rhoINIT   = new double[ns];
  rhoINJECT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i = 0; i < ns; i++) {
    rhoINIT  [i] = col->getRHOinit(i);
    rhoINJECT[i] = col->getRHOinject(i);
    if ((fabs(col->getW0(i)) != 0) || (fabs(col->getU0(i)) != 0)) // GEM and LHDI
      DriftSpecies[i] = true;
    else
      DriftSpecies[i] = false;
  }
  /*! parameters for GEM challenge */
  FourPI = 16 * atan(1.0);
  /*! Restart */
  restart1 = col->getRestart_status();
  RestartDirName = col->getRestartDirName();
  Case = col->getCase();

  // OpenBC
  injFieldsLeft   = new injInfoFields(nxn, nyn, nzn);
  injFieldsRight  = new injInfoFields(nxn, nyn, nzn);
  injFieldsTop    = new injInfoFields(nxn, nyn, nzn);
  injFieldsBottom = new injInfoFields(nxn, nyn, nzn);
  injFieldsFront  = new injInfoFields(nxn, nyn, nzn);
  injFieldsRear   = new injInfoFields(nxn, nyn, nzn);

  // arrays allocation: nodes
  Ex = newArr3(double, nxn, nyn, nzn);
  Ey = newArr3(double, nxn, nyn, nzn);
  Ez = newArr3(double, nxn, nyn, nzn);
  Exth = newArr3(double, nxn, nyn, nzn);
  Eyth = newArr3(double, nxn, nyn, nzn);
  Ezth = newArr3(double, nxn, nyn, nzn);
  Bxn = newArr3(double, nxn, nyn, nzn);
  Byn = newArr3(double, nxn, nyn, nzn);
  Bzn = newArr3(double, nxn, nyn, nzn);
  rhon = newArr3(double, nxn, nyn, nzn);
  Jx = newArr3(double, nxn, nyn, nzn);
  Jy = newArr3(double, nxn, nyn, nzn);
  Jz = newArr3(double, nxn, nyn, nzn);
  Jxh = newArr3(double, nxn, nyn, nzn);
  Jyh = newArr3(double, nxn, nyn, nzn);
  Jzh = newArr3(double, nxn, nyn, nzn);
  // External imposed fields
  Bx_ext = newArr3(double,nxn,nyn,nzn);
  By_ext = newArr3(double,nxn,nyn,nzn);
  Bz_ext = newArr3(double,nxn,nyn,nzn);

  // i initialize it now just to shut up valgrind
  for (int i=0; i< nxn; i++)
    for (int j=0; j< nyn; j++)
      for (int k=0; k< nzn; k++){
	Bx_ext[i][j][k]=0.0;
	By_ext[i][j][k]=0.0;
	Bz_ext[i][j][k]=0.0;
      }
  // Jx_ext = newArr3(double,nxn,nyn,nzn);
  // Jy_ext = newArr3(double,nxn,nyn,nzn);
  // Jz_ext = newArr3(double,nxn,nyn,nzn);
  // involving species
  //rhons = newArr4(double, ns, nxn, nyn, nzn);
  rhons = newArr4(double, ns*2, nxn, nyn, nzn); // to save stuff    

  rhocs = newArr4(double, ns, nxc, nyc, nzc);

  RHOINIT = newArr4(double, ns, nxc, nyc, nzc);

  // debug, to remove (you can use this to save stuff)
  Jxs = newArr4(double, ns, nxn, nyn, nzn);
  Jys = newArr4(double, ns, nxn, nyn, nzn);
  Jzs = newArr4(double, ns, nxn, nyn, nzn);

  pXXsn = newArr4(double, ns, nxn, nyn, nzn);
  pXYsn = newArr4(double, ns, nxn, nyn, nzn);
  pXZsn = newArr4(double, ns, nxn, nyn, nzn);
  pYYsn = newArr4(double, ns, nxn, nyn, nzn);
  pYZsn = newArr4(double, ns, nxn, nyn, nzn);
  pZZsn = newArr4(double, ns, nxn, nyn, nzn);
  // arrays allocation: central points 
  PHI = newArr3(double, nxc, nyc, nzc);
  Bxc = newArr3(double, nxc, nyc, nzc);
  Byc = newArr3(double, nxc, nyc, nzc);
  Bzc = newArr3(double, nxc, nyc, nzc);
  rhoc = newArr3(double, nxc, nyc, nzc);
  rhoh = newArr3(double, nxc, nyc, nzc);

  // temporary arrays
  tempXC = newArr3(double, nxc, nyc, nzc);
  tempYC = newArr3(double, nxc, nyc, nzc);
  tempZC = newArr3(double, nxc, nyc, nzc);

  tempXN = newArr3(double, nxn, nyn, nzn);
  tempYN = newArr3(double, nxn, nyn, nzn);
  tempZN = newArr3(double, nxn, nyn, nzn);
  tempC = newArr3(double, nxc, nyc, nzc);
  tempX = newArr3(double, nxn, nyn, nzn);
  tempY = newArr3(double, nxn, nyn, nzn);
  tempZ = newArr3(double, nxn, nyn, nzn);
  temp2X = newArr3(double, nxn, nyn, nzn);
  temp2Y = newArr3(double, nxn, nyn, nzn);
  temp2Z = newArr3(double, nxn, nyn, nzn);
  imageX = newArr3(double, nxn, nyn, nzn);
  imageY = newArr3(double, nxn, nyn, nzn);
  imageZ = newArr3(double, nxn, nyn, nzn);
  imageEX = newArr3(double, nxc, nyc, nzc);
  imageEY = newArr3(double, nxc, nyc, nzc);
  imageEZ = newArr3(double, nxc, nyc, nzc);
  imageBX = newArr3(double, nxn, nyn, nzn);
  imageBY = newArr3(double, nxn, nyn, nzn);
  imageBZ = newArr3(double, nxn, nyn, nzn);
  Dx = newArr3(double, nxn, nyn, nzn);
  Dy = newArr3(double, nxn, nyn, nzn);
  Dz = newArr3(double, nxn, nyn, nzn);
  vectX = newArr3(double, nxn, nyn, nzn);
  vectY = newArr3(double, nxn, nyn, nzn);
  vectZ = newArr3(double, nxn, nyn, nzn);
  divC = newArr3(double, nxc, nyc, nzc);
  arr = newArr3(double,nxn,nyn,nzn);

  Lambda = newArr3(double, nxn, nyn, nzn);

  /*! mlmd specific section */
  MLMDVerbose= col->getMLMDVerbose();

  /* wether to perform mlmd operations */
  MLMD_BC = col->getMLMD_BC();
  MLMD_PROJECTION = col->getMLMD_PROJECTION();
  ParticleREPOPULATION = col->getMLMD_ParticleREPOPULATION();
  MLMD_BCBufferArea = col->getMLMD_BCBufferArea();
  MLMD_fixBCenters = col->getMLMD_fixBCenters();
  MLMD_InterpolateOldBCell = col->getMLMD_InterpolateOldBCell();
  if (MLMD_InterpolateOldBCell){
    MLMD_fixBCenters=false;
  }

  /* # of fields that I am sending as BC: 9 (Ex, Ey, Ez, Exth, Eyth, Ezth, Bxn, Byn, Bzn) */
  NumF= 9; 
  /* # of fields to send for MLMD_fixBCenters: Bxc, Byc, Bzc */
  Numfix3B= 3; 

  /* create mlmd-related mpi datatypes */
  MPI_RGBC_struct_commit();

  /* end create mlmd-related mpi datatypes */
  /*! end mlmd specidic section */

  /* for mlmd BC communication */
  TAG_BC_GHOST= 2;
  TAG_BC_ACTIVE= 3;
  TAG_PROJ= 4;
  TAG_II= 5;
  TAG_BC_BUFFER= 7;
  TAG_BC_FIX3B= 9;
  /* stuff that it's convenient not to pass around */
  // since the number of core per grid may vary wildly and give problems
  // size message buffer on the max size
  int MaxGridCoreN= vct->getMaxGridCoreN(); // NB: this is only local

  /*if (vct->getCartesian_rank() == 0)
    cout << "In EMfields3D.cpp, grid " << numGrid << ", MaxGridCoreN is " << MaxGridCoreN << endl;*/
  
  MAX_RG_numBCMessages= (int) (MaxGridCoreN*6+1); // something smarter has to be done with this guy
  MAX_size_LevelWide= MAX_RG_numBCMessages* 4;
  
  int P= vct->getParentGridNum();
  // resolution of the parent
  if (P>-1){
    DxP= grid->getDx_mlmd(P);
    DyP= grid->getDy_mlmd(P);
    DzP= grid->getDz_mlmd(P);
    // refinement jump per grid
    RFx= DxP/dx;
    RFy= DyP/dy;
    RFz= DzP/dz;

    // NB: one cell in the inputfile means nzn=4
    if ((ceil(RFx) > nxn-1 and nxn>4)  or (ceil(RFy) > nyn-1 and nyn>4) or (ceil(RFz) > nzn-1 and nzn>4)){
      cout << "The number of cells per core must be > RF, aborting..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }else { // if these guys are used in grid 0, i try to force an error
    DxP=0; DyP=0; DzP=0;
    RFx=0; RFy=0; RFz=0;
  }

  /* length over which to do the buffering */
  BufX= col->getBuf();
  BufY= col->getBuf();
  BufZ= col->getBuf();

  BSNeeded= false;


  RepopulateBeforeMover= col->getRepopulateBeforeMover();

  if (vct->getCartesian_rank() == 0){
    cout << "Grid " << numGrid <<", nxn= " << nxn << " ,BufX= " << BufX <<", nyn= " << nyn << " ,BufY= "<< BufY <<", nzn= " << nzn << " ,BufZ= "<< BufZ << endl;
  }

  SmoothFaces= false;

  CommToParent_InDel= vct->getCommToParent();
}

/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::startEcalc(Grid * grid, VirtualTopology3D * vct, Collective *col) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid <<": E CALCULATION ***" << endl;
  double ***divE = newArr3(double, nxc, nyc, nzc);
  double ***gradPHIX = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIY = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIZ = newArr3(double, nxn, nyn, nzn);

  double *xkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  double *bkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  // set to zero all the stuff 
  eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, divE, nxc, nyc, nzc);
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, gradPHIX, nxn, nyn, nzn);
  eqValue(0.0, gradPHIY, nxn, nyn, nzn);
  eqValue(0.0, gradPHIZ, nxn, nyn, nzn);
  // Adjust E calculating laplacian(PHI) = div(E) -4*PI*rho DIVERGENCE CLEANING
  if (PoissonCorrection) {
    if (vct->getCartesian_rank() == 0)
      cout << "*** G" << numGrid <<": DIVERGENCE CLEANING ***" << endl;
    grid->divN2C(divE, Ex, Ey, Ez);
    scale(tempC, rhoc, -FourPI, nxc, nyc, nzc);
    sum(divE, tempC, nxc, nyc, nzc);
    // move to krylov space 
    phys2solver(bkrylovPoisson, divE, nxc, nyc, nzc);
    // use conjugate gradient first
    if (!CG(xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 3000, CGtol, &Field::PoissonImage, grid, vct, this)) {
      if (vct->getCartesian_rank() == 0)
        cout << "CG not Converged. Trying with GMRes. Consider to increase the number of the CG iterations" << endl;
      eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
      GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 20, 200, GMREStol, grid, vct, this);
    }
    solver2phys(PHI, xkrylovPoisson, nxc, nyc, nzc);
    communicateCenterBC(nxc, nyc, nzc, PHI, 2, 2, 2, 2, 2, 2, vct);
    // calculate the gradient
    grid->gradC2N(gradPHIX, gradPHIY, gradPHIZ, PHI);
    // sub
    sub(Ex, gradPHIX, nxn, nyn, nzn);
    sub(Ey, gradPHIY, nxn, nyn, nzn);
    sub(Ez, gradPHIZ, nxn, nyn, nzn);

  }                             // end of divergence cleaning 

  delete[]xkrylovPoisson;
  delete[]bkrylovPoisson;
  delArr3(divE, nxc, nyc);
  delArr3(gradPHIX, nxn, nyn);
  delArr3(gradPHIY, nxn, nyn);
  delArr3(gradPHIZ, nxn, nyn);

  if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){
    if (!MLMD_InterpolateOldBCell){ 
      if (vct->getCartesian_rank() == 0)
	cout << "I AM NOT APPLYING B BC" << endl;
    }
    else{ /* with this command, I am interpolating all B^n cell values */
      if (vct->getCartesian_rank() == 0)
	cout << "I am interpolating all B^n cell values" << endl;


      //cout << "RG_numBCMessages_fix3B: "<< RG_numBCMessages_fix3B << endl;
      setBC_Nodes(vct, Bxc, Byc, Bzc, Bxc_fix3B_BC, Byc_fix3B_BC, Bzc_fix3B_BC, RGBC_Info_fix3B, RG_numBCMessages_fix3B);
    }
  }
  // set B BC
  /*if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC and !MLMD_fixBCenters){ // this added here when removing the BC exchange before B
    setBC_Nodes(vct, Bxn, Byn, Bzn, Bxn_Ghost_BC, Byn_Ghost_BC, Bzn_Ghost_BC, RGBC_Info_Ghost, RG_numBCMessages_Ghost);
    setBC_Nodes(vct, Bxn, Byn, Bzn, Bxn_Active_BC, Byn_Active_BC, Bzn_Active_BC, RGBC_Info_Active, RG_numBCMessages_Active);
    
    // buffer interpolation
    if (MLMD_BCBufferArea){
      setBC_Nodes(vct, Bxn, Byn, Bzn, Bxn_Buffer_BC, Byn_Buffer_BC, Bzn_Buffer_BC, RGBC_Info_Buffer, RG_numBCMessages_Buffer);
    }

    int NcellsX=2;
    int NcellsY=2;
    int NcellsZ=2;
            
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);

    
    // fix cells
    if (vct->getXleft_neighbor() == MPI_PROC_NULL)
      fixBghostCells_Left(0, NcellsX);
    if (vct->getYleft_neighbor() == MPI_PROC_NULL)
      fixBghostCells_Left(1, NcellsY);
    if (vct->getZleft_neighbor() == MPI_PROC_NULL)
      fixBghostCells_Left(2, NcellsZ);
    
    if (vct->getXright_neighbor() == MPI_PROC_NULL)
      fixBghostCells_Right(0, NcellsX);
    if (vct->getYright_neighbor() == MPI_PROC_NULL)
      fixBghostCells_Right(1, NcellsY);
    if (vct->getZright_neighbor() == MPI_PROC_NULL)
    fixBghostCells_Right(2, NcellsZ);
  }

  
  
  if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC and MLMD_fixBCenters){ 
    setBC_Nodes(vct, Bxn, Byn, Bzn, Bxn_Ghost_BC, Byn_Ghost_BC, Bzn_Ghost_BC, RGBC_Info_Ghost, RG_numBCMessages_Ghost);

    if (vct->getCartesian_rank()==0)
      cout << "FixB in refined grid " << endl;
    // in this case it actually sets the centers
    setBC_Nodes(vct, Bxc, Byc, Bzc, Bxc_fix3B_BC, Byc_fix3B_BC, Bzc_fix3B_BC, RGBC_Info_fix3B, RG_numBCMessages_fix3B); 
  
  }
  */
  
  //correctDivB(grid, vct);


}

void EMfields3D::calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col) {

  startEcalc(grid,vct, col);
  
  double *xkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 E components
  double *bkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 components
  eqValue(0.0, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, bkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));

  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid <<": MAXWELL SOLVER ***" << endl;


  // prepare the source 
  MaxwellSource(bkrylov, grid, vct, col);

  if (!(Case== "initTestIntProj" and numGrid ==0)){
    phys2solver(xkrylov, Ex, Ey, Ez, nxn, nyn, nzn);
  }else{
    for (int i=0; i< 3 * (nxn - 2) * (nyn - 2) * (nzn - 2); i++){
      xkrylov[i]=0.0;
    }
  }
  // solver

  if (Case!="MAX_Show_RG_BC" and Case !="initTestBC" and Case!="TestBBoundary"  and Case!= "initTestIntProj"  ){
    GMRES(&Field::MaxwellImage, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2), bkrylov, 20, 200, GMREStol, grid, vct, this);
  }else{
    cout << "Grid " << numGrid <<"No GMRES" << endl;
  }

  endEcalc(xkrylov, grid, vct, col);

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;

  if (numGrid==0) return;
  return;
  
  for (int i=0; i< nxn; i++)
    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++){

	if (i==0 or j==0 or k==0 or i==nxn-1 or j==nyn-1 or k==nzn-1){
	  if (fabs(Ezth[i][j][k])> DBL_EPSILON) {
	    cout << "ABORTING" <<endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	  //else cout << "BC OK" << endl;
	}

	if (i==1 or j==1 or k==1 or i==nxn-2 or j==nyn-2 or k==nzn-2){
	  if (fabs(Ezth[i][j][k])> DBL_EPSILON) {
	    cout << "ABORTING" <<endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	  //else cout << "BC OK" << endl;
	}

	/*if (i==0 or j==0 or k==0 or i==nxn-1 or j==nyn-1 or k==nzn-1){
	  cout << "Multiplying E ghost x 2" << endl;
	  Exth[i][j][k]= Exth[i][j][k]*2;
	  Eyth[i][j][k]= Eyth[i][j][k]*2;
	  Ezth[i][j][k]= Ezth[i][j][k]*2;
	  }*/

      }

}
 
void EMfields3D::endEcalc(double* xkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col)
{
  
  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, xkrylov, nxn, nyn, nzn);

  /* here i used to set E BC, but I have removed that */
  /* I put back ghost E th BC just for the smoothing */
  if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){ // mlmd

     setBC_Nodes(vct, Exth, Eyth, Ezth, Exth_Ghost_BC, Eyth_Ghost_BC, Ezth_Ghost_BC, RGBC_Info_Ghost, RG_numBCMessages_Ghost);

  }

  // save the old E if you need to regenerate it after receiving projection
  if (numChildren>0 and MLMD_PROJECTION and ApplyProjection){
    for (int i=0; i< nxn; i++)
      for (int j=0; j< nyn; j++)
	for (int k=0; k< nzn; k++){
	  Ex_n[i][j][k]= Ex[i][j][k];
	  Ey_n[i][j][k]= Ey[i][j][k];
	  Ez_n[i][j][k]= Ez[i][j][k];
	}
  } // end if (numChildren>0 and MLMD_BC and ApplyProjection){

  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);  
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);

  // this is a copy of E BEFORE SMOOTHING to use for BC
  if (BSNeeded){
    for (int i=0; i< nxn; i++)
      for (int j=0; j< nyn; j++)
	 for (int k=0; k< nzn; k++){
	   Ex_BS[i][j][k]= Ex[i][j][k];
	   Ey_BS[i][j][k]= Ey[i][j][k];
	   Ez_BS[i][j][k]= Ez[i][j][k];
	 }
   }

   smoothE(Smooth, vct, col);

   /* this used to be only for grid 0, now for everybody */
   communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

   communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
   
   // OpenBC                                                                      
   
   if (!( vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC)){
     
     BoundaryConditionsE(Exth, Eyth, Ezth, nxn, nyn, nzn, grid, vct);
     BoundaryConditionsE(Ex, Ey, Ez, nxn, nyn, nzn, grid, vct);
     
   }


 }


 /*! Calculate sorgent for Maxwell solver */
 void EMfields3D::MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col) {

   // ME: remember bkrylov is stripped of the ghosts

   eqValue(0.0, tempC, nxc, nyc, nzc);
   eqValue(0.0, tempX, nxn, nyn, nzn);
   eqValue(0.0, tempY, nxn, nyn, nzn);
   eqValue(0.0, tempZ, nxn, nyn, nzn);
   eqValue(0.0, tempXN, nxn, nyn, nzn);
   eqValue(0.0, tempYN, nxn, nyn, nzn);
   eqValue(0.0, tempZN, nxn, nyn, nzn);
   eqValue(0.0, temp2X, nxn, nyn, nzn);
   eqValue(0.0, temp2Y, nxn, nyn, nzn);
   eqValue(0.0, temp2Z, nxn, nyn, nzn);

   // communicate
   if (! (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC)){  // non mlmd

     // "normal option"
     communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
     communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
     communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);


     if (Case=="ForceFree") fixBforcefree(grid,vct);
     if (Case=="GEM")       fixBgem(grid, vct);
     if (Case=="GEMnoPert") fixBgem(grid, vct);

     // OpenBC:
     BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc,grid,vct);
   }

   if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC and MLMD_fixBCenters){
     if (vct->getCartesian_rank()==0)
       cout << "FixB in refined grid " << endl;
     // in this case it actually sets the centers
     setBC_Nodes(vct, Bxc, Byc, Bzc, Bxc_fix3B_BC, Byc_fix3B_BC, Bzc_fix3B_BC, RGBC_Info_fix3B, RG_numBCMessages_fix3B);

   }

   // prepare curl of B for known term of Maxwell solver: for the source term
   grid->curlC2N(tempXN, tempYN, tempZN, Bxc, Byc, Bzc);   // the ghost node here is not ok for the mlmd, but it's never used and scraped before entering GMRES

   scale(temp2X, Jxh, -FourPI / c, nxn, nyn, nzn);
   scale(temp2Y, Jyh, -FourPI / c, nxn, nyn, nzn);
   scale(temp2Z, Jzh, -FourPI / c, nxn, nyn, nzn);

   // -- dipole SOURCE version using J_ext
   // addscale(-FourPI/c,temp2X,Jx_ext,nxn,nyn,nzn);
   // addscale(-FourPI/c,temp2Y,Jy_ext,nxn,nyn,nzn);
   // addscale(-FourPI/c,temp2Z,Jz_ext,nxn,nyn,nzn);
   // -- end of dipole SOURCE version using J_ext

   // curl(B) - 4pi J_hat
   sum(temp2X, tempXN, nxn, nyn, nzn);
   sum(temp2Y, tempYN, nxn, nyn, nzn);
   sum(temp2Z, tempZN, nxn, nyn, nzn);  

   // cdt [curl(B) - 4pi J_hat]
   scale(temp2X, delt, nxn, nyn, nzn);
   scale(temp2Y, delt, nxn, nyn, nzn);
   scale(temp2Z, delt, nxn, nyn, nzn);


   communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct);
   grid->gradC2N(tempX, tempY, tempZ, rhoh);  

   // cdt^2 * 4pi * grad(rho_hat)
   scale(tempX, -delt * delt * FourPI, nxn, nyn, nzn);
   scale(tempY, -delt * delt * FourPI, nxn, nyn, nzn);
   scale(tempZ, -delt * delt * FourPI, nxn, nyn, nzn); 
   // end tmp

   // sum E, past values:  E + cdt^2 * 4pi * grad(rho_hat)
   sum(tempX, Ex, nxn, nyn, nzn);
   sum(tempY, Ey, nxn, nyn, nzn);
   sum(tempZ, Ez, nxn, nyn, nzn);

   // cdt [curl(B) - 4pi J_hat]  +  E + cdt^2 * 4pi * grad(rho_hat)
   sum(tempX, temp2X, nxn, nyn, nzn);
   sum(tempY, temp2Y, nxn, nyn, nzn);
   sum(tempZ, temp2Z, nxn, nyn, nzn);

   if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){ // mlmd

     // put the BC on the soure and in the image = intermediate solution
     setBC_Nodes(vct, tempX, tempY, tempZ, Exth_Active_BC, Eyth_Active_BC, Ezth_Active_BC, RGBC_Info_Active, RG_numBCMessages_Active);

     //     divECorrection_AllFaces(tempX, tempY, tempZ, grid, vct);

     if (MLMD_BCBufferArea){
      setBC_Nodes(vct, tempX, tempY, tempZ, Exth_Buffer_BC, Eyth_Buffer_BC, Ezth_Buffer_BC, RGBC_Info_Buffer, RG_numBCMessages_Buffer);
     }

     // smoothFaces
     if (SmoothFaces and Smooth> 0.9999){
       smoothFaces(Smooth, tempX, grid, vct, 1);
       smoothFaces(Smooth, tempY, grid, vct, 1);
       smoothFaces(Smooth, tempZ, grid, vct, 1);
     }
   } else { // normal stuff


     // Boundary condition in the known term
     // boundary condition: Xleft
     if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
       perfectConductorLeftS(tempX, tempY, tempZ, 0);
     // boundary condition: Xright
     if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
       perfectConductorRightS(tempX, tempY, tempZ, 0);
     // boundary condition: Yleft
     if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
       perfectConductorLeftS(tempX, tempY, tempZ, 1);
     // boundary condition: Yright
     if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
       perfectConductorRightS(tempX, tempY, tempZ, 1);
     // boundary condition: Zleft
     if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
       perfectConductorLeftS(tempX, tempY, tempZ, 2);
     // boundary condition: Zright
     if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
       perfectConductorRightS(tempX, tempY, tempZ, 2);

   }


   // physical space -> Krylov space (cuts the guard cells)
   phys2solver(bkrylov, tempX, tempY, tempZ, nxn, nyn, nzn);


 }
 /*! Mapping of Maxwell image to give to solver */
 void EMfields3D::MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct) {
   eqValue(0.0, im, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
   eqValue(0.0, imageX, nxn, nyn, nzn);
   eqValue(0.0, imageY, nxn, nyn, nzn);
   eqValue(0.0, imageZ, nxn, nyn, nzn);
   eqValue(0.0, tempX, nxn, nyn, nzn);
   eqValue(0.0, tempY, nxn, nyn, nzn);
   eqValue(0.0, tempZ, nxn, nyn, nzn);
   eqValue(0.0, Dx, nxn, nyn, nzn);
   eqValue(0.0, Dy, nxn, nyn, nzn);
   eqValue(0.0, Dz, nxn, nyn, nzn);


   // move from krylov space to physical space
   solver2phys(vectX, vectY, vectZ, vector, nxn, nyn, nzn); 

   grid->lapN2N(imageX, vectX ,vct); 
   grid->lapN2N(imageY, vectY, vct);
   grid->lapN2N(imageZ, vectZ, vct);

   neg(imageX, nxn, nyn, nzn);
   neg(imageY, nxn, nyn, nzn);
   neg(imageZ, nxn, nyn, nzn);
   // grad(div(mu dot E(n + theta)) mu dot E(n + theta) = D
   // MLMD: here ghost nodes is compromised by mu
   MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, grid, vct);
   // MLMD: here ghost center is compromised by ghost node
   grid->divN2C(divC, Dx, Dy, Dz);
   // communicate you should put BC 
   // think about the Physics 
   // communicateCenterBC(nxc,nyc,nzc,divC,1,1,1,1,1,1,vct);



   if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){
     communicateCenter(nxc, nyc, nzc, divC, vct);
   }else{
     communicateCenterBC(nxc, nyc, nzc, divC, 2, 2, 2, 2, 2, 2, vct);  // GO with Neumann, now then go with rho
   }

   // MLMD: here first active node is compromised
   grid->gradC2N(tempX, tempY, tempZ, divC);

   // -lap(E(n +theta)) - grad(div(mu dot E(n + theta))
   sub(imageX, tempX, nxn, nyn, nzn);
   sub(imageY, tempY, nxn, nyn, nzn);
   sub(imageZ, tempZ, nxn, nyn, nzn);

   // scale delt*delt
   scale(imageX, delt * delt, nxn, nyn, nzn);
   scale(imageY, delt * delt, nxn, nyn, nzn);
   scale(imageZ, delt * delt, nxn, nyn, nzn);

   // -lap(E(n +theta)) - grad(div(mu dot E(n + theta))) + eps dot E(n + theta)
   sum(imageX, Dx, nxn, nyn, nzn);
   sum(imageY, Dy, nxn, nyn, nzn);
   sum(imageZ, Dz, nxn, nyn, nzn);
   sum(imageX, vectX, nxn, nyn, nzn);
   sum(imageY, vectY, nxn, nyn, nzn);
   sum(imageZ, vectZ, nxn, nyn, nzn);

   // Temporal damping
   sumscalprod(imageX, delt, vectX, Lambda, nxn, nyn, nzn);
   sumscalprod(imageY, delt, vectY, Lambda, nxn, nyn, nzn);
   sumscalprod(imageZ, delt, vectZ, Lambda, nxn, nyn, nzn);

   // MLMD: first active node is compromised; but then replaced with BC

   if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){
     //cout << "setBC_NodesImage in Active " << endl;
     //set what you have to set, in the active nodes
     setBC_NodesImage(vct, imageX, imageY, imageZ, vectX, vectY, vectZ, Exth_Active_BC, Eyth_Active_BC, Ezth_Active_BC, RGBC_Info_Active, RG_numBCMessages_Active);

     if (MLMD_BCBufferArea){ 
       setBC_NodesImage(vct, imageX, imageY, imageZ, vectX, vectY, vectZ, Exth_Buffer_BC, Eyth_Buffer_BC, Ezth_Buffer_BC, RGBC_Info_Buffer, RG_numBCMessages_Buffer);
     }

   }else {

     // boundary condition: Xleft
     if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
       perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 0, grid);
     // boundary condition: Xright
     if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
       perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 0, grid);
     // boundary condition: Yleft
     if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
       perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 1, grid);
     // boundary condition: Yright
     if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
       perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 1, grid);
     // boundary condition: Zleft
     if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
       perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 2, grid);
     // boundary condition: Zright
     if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
       perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 2, grid);

     // OpenBC
     BoundaryConditionsEImage(imageX, imageY, imageZ, vectX, vectY, vectZ, nxn, nyn, nzn, vct, grid);
   }// end else


    // move from physical space to krylov space
   phys2solver(im, imageX, imageY, imageZ, nxn, nyn, nzn);


 }

 /*! Calculate PI dot (vectX, vectY, vectZ) */
 void EMfields3D::PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid * grid) {
   double beta, edotb, omcx, omcy, omcz, denom;
   beta = .5 * qom[ns] * dt / c;
   for (int i = 1; i < nxn - 1; i++)
     for (int j = 1; j < nyn - 1; j++)
       for (int k = 1; k < nzn - 1; k++) {
	 omcx = beta * (Bxn[i][j][k] + Fext*Bx_ext[i][j][k]);
	 omcy = beta * (Byn[i][j][k] + Fext*By_ext[i][j][k]);
	 omcz = beta * (Bzn[i][j][k] + Fext*Bz_ext[i][j][k]);
	 edotb = vectX[i][j][k] * omcx + vectY[i][j][k] * omcy + vectZ[i][j][k] * omcz;
	 denom = 1 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
	 PIdotX[i][j][k] += (vectX[i][j][k] + (vectY[i][j][k] * omcz - vectZ[i][j][k] * omcy + edotb * omcx)) * denom;
	 PIdotY[i][j][k] += (vectY[i][j][k] + (vectZ[i][j][k] * omcx - vectX[i][j][k] * omcz + edotb * omcy)) * denom;
	 PIdotZ[i][j][k] += (vectZ[i][j][k] + (vectX[i][j][k] * omcy - vectY[i][j][k] * omcx + edotb * omcz)) * denom;
       }


 }
 /*! Calculate MU dot (vectX, vectY, vectZ) */
 void EMfields3D::MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid * grid, VirtualTopology3D * vct) {

   double beta, edotb, omcx, omcy, omcz, denom;
   for (int i = 0; i < nxn ; i++)
     for (int j = 0; j < nyn ; j++)
       for (int k = 0; k < nzn ; k++) {
	 MUdotX[i][j][k] = 0.0;
	 MUdotY[i][j][k] = 0.0;
	 MUdotZ[i][j][k] = 0.0;
       }
   for (int is = 0; is < ns; is++) {
     beta = .5 * qom[is] * dt / c;
     for (int i = 0; i < nxn ; i++)
       for (int j = 0; j < nyn ; j++)
	 for (int k = 0; k < nzn ; k++) {
	   omcx = beta * (Bxn[i][j][k] + Fext*Bx_ext[i][j][k]);
	   omcy = beta * (Byn[i][j][k] + Fext*By_ext[i][j][k]);
	   omcz = beta * (Bzn[i][j][k] + Fext*Bz_ext[i][j][k]);
	   edotb = vectX[i][j][k] * omcx + vectY[i][j][k] * omcy + vectZ[i][j][k] * omcz; 


	   denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);


	   MUdotX[i][j][k] += (vectX[i][j][k] + (vectY[i][j][k] * omcz - vectZ[i][j][k] * omcy + edotb * omcx)) * denom;
	   MUdotY[i][j][k] += (vectY[i][j][k] + (vectZ[i][j][k] * omcx - vectX[i][j][k] * omcz + edotb * omcy)) * denom;
	   MUdotZ[i][j][k] += (vectZ[i][j][k] + (vectX[i][j][k] * omcy - vectY[i][j][k] * omcx + edotb * omcz)) * denom;
	 }

   }

 }
 /* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
 void EMfields3D::smooth(double value, double ***vector, int type, Grid * grid, VirtualTopology3D * vct) {

   if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

   int nvolte = 6;
   for (int icount = 1; icount < nvolte + 1; icount++) {

     if (value != 1.0) {
       double alpha;
       int nx, ny, nz;
       switch (type) {
       case (0):
	 nx = grid->getNXC();
	 ny = grid->getNYC();
	 nz = grid->getNZC();
	 communicateCenterBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);

	 break;
       case (1):
	 nx = grid->getNXN();
	 ny = grid->getNYN();
	 nz = grid->getNZN();
	 communicateNodeBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);
	 break;
       }
       double ***temp = newArr3(double, nx, ny, nz);
       if (icount % 2 == 1) {
	 value = 0.;
       }
       else {
	 value = 0.5;
       }
       alpha = (1.0 - value) / 6;
       for (int i = 1; i < nx - 1; i++)
	 for (int j = 1; j < ny - 1; j++)
	   for (int k = 1; k < nz - 1; k++)
	     temp[i][j][k] = value * vector[i][j][k] + alpha * (vector[i - 1][j][k] + vector[i + 1][j][k] + vector[i][j - 1][k] + vector[i][j + 1][k] + vector[i][j][k - 1] + vector[i][j][k + 1]);
       for (int i = 1; i < nx - 1; i++)
	 for (int j = 1; j < ny - 1; j++)
	   for (int k = 1; k < nz - 1; k++)
	     vector[i][j][k] = temp[i][j][k];
       delArr3(temp, nx, ny);
     }
   }
 }
 /* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
 void EMfields3D::smoothE(double value, VirtualTopology3D * vct, Collective *col) {

   if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

   //NB: ghost is not smoothed
   bool ExtraSmoothing= false;

   ExtraSmoothing= false; // extra smoothing is actually worse, so DO NOT enable it 

   if (vct->getCartesian_rank()==0 and numGrid >0 and ExtraSmoothing){
     cout << "Grid " << numGrid <<" is doing extra boundary smoothing" << endl;
   }
   int nvolte = 6;

   int i_s, i_e;
   int j_s, j_e;
   int k_s, k_e;
   for (int icount = 1; icount < nvolte + 1; icount++) {
     if (value != 1.0) {
       double alpha;
       
       /* without setting ghosts in RG 
       communicateNodeBoxStencilBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
	 communicateNodeBoxStencilBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
	 communicateNodeBoxStencilBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);*/
	 
       /* if ghosts in the RG are set*/
       if (! MLMD_BC){ // CG 
	 communicateNodeBoxStencilBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
	 communicateNodeBoxStencilBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
	 communicateNodeBoxStencilBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
       } 
       else { //RG
	 communicateNodeBoxStencil(nxn, nyn, nzn, Ex, vct);
	 communicateNodeBoxStencil(nxn, nyn, nzn, Ey, vct);
	 communicateNodeBoxStencil(nxn, nyn, nzn, Ez, vct);
       }

	
       double ***temp = newArr3(double, nxn, nyn, nzn);
       if (icount % 2 == 1) {
	 value = 0.;
       }
       else {
	 value = 0.5;
       }
       alpha = (1.0 - value) / 6;

       // RG BC need to be smoothed

       i_s=1; i_e=nxn-1;
       j_s=1; j_e=nyn-1;
       k_s=1; k_e=nzn-1;
       

       int BBX= RFx*4;
       int BBY= RFy*4;
       int BBZ= RFz*4;

       // Exth
       for (int i = i_s; i < i_e; i++)
	 for (int j = j_s; j < j_e; j++)
	   for (int k = k_s; k < k_e; k++){

	     if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){

	       value= 0;

	     }
	     else{
	       if (icount % 2 == 1) {
		 value = 0.;
	       }
	       else {
		 value = 0.5;
	       }

	     }
	     alpha = (1.0 - value) / 6;

	     temp[i][j][k] = value * Ex[i][j][k] + alpha * (Ex[i - 1][j][k] + Ex[i + 1][j][k] + Ex[i][j - 1][k] + Ex[i][j + 1][k] + Ex[i][j][k - 1] + Ex[i][j][k + 1]);}

       for (int i = i_s; i < i_e; i++)
         for (int j = j_s; j < j_e; j++)
           for (int k = k_s; k < k_e; k++)
	     Ex[i][j][k] = temp[i][j][k];
       // Eyth
       for (int i = i_s; i < i_e; i++)
         for (int j = j_s; j < j_e; j++)
           for (int k = k_s; k < k_e; k++){
	     
	     if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){
	       
	       
	       value= 0;
	       
	     }
	     else{
	       if (icount % 2 == 1) {
		 value = 0.;
	       }
	       else {
		 value = 0.5;
	       }
	       
	     }
	     alpha = (1.0 - value) / 6;
	     
	     temp[i][j][k] = value * Ey[i][j][k] + alpha * (Ey[i - 1][j][k] + Ey[i + 1][j][k] + Ey[i][j - 1][k] + Ey[i][j + 1][k] + Ey[i][j][k - 1] + Ey[i][j][k + 1]);
	   }
       for (int i = i_s; i < i_e; i++)
	 for (int j = j_s; j < j_e; j++)
	   for (int k = k_s; k < k_e; k++){
	     
	     Ey[i][j][k] = temp[i][j][k];
     }
     // Ezth
     for (int i = i_s; i < i_e; i++)
       for (int j = j_s; j < j_e; j++)
	 for (int k = k_s; k < k_e; k++){
	   
	   
	   if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){
	     
	     
	     value= 0;
	     
	   }
	   else{
	     if (icount % 2 == 1) {
	       value = 0.;
	     }
	     else {
	       value = 0.5;
	     }
	     
	   }
	   alpha = (1.0 - value) / 6;
	   
	   temp[i][j][k] = value * Ez[i][j][k] + alpha * (Ez[i - 1][j][k] + Ez[i + 1][j][k] + Ez[i][j - 1][k] + Ez[i][j + 1][k] + Ez[i][j][k - 1] + Ez[i][j][k + 1]);
	 }
     for (int i = i_s; i < i_e; i++)
       for (int j = j_s; j < j_e; j++)
	 for (int k = k_s; k < k_e; k++)
	   
	   
	   Ez[i][j][k] = temp[i][j][k];
     
     
     delArr3(temp, nxn, nyn);
   }
 }
}

/**/
void EMfields3D::smoothE_NoComm(double value, VirtualTopology3D * vct, Collective *col) {

  if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }


   bool ExtraSmoothing= false;


   if (numGrid>0) ExtraSmoothing= true;

   ExtraSmoothing= false; // extra smoothing is actually worse, so DO NOT enable it 

   if (vct->getCartesian_rank()==0 and numGrid >0 and ExtraSmoothing){
     cout << "Grid " << numGrid <<" is doing extra boundary smoothing" << endl;
   }
   int nvolte = 6;
   for (int icount = 1; icount < nvolte + 1; icount++) {
     if (value != 1.0) {
       double alpha;

       double ***temp = newArr3(double, nxn, nyn, nzn);
       if (icount % 2 == 1) {
	 value = 0.;
       }
       else {
	 value = 0.5;
       }
       alpha = (1.0 - value) / 6;

       // not to blur active node solution in the RG
       int i_s=1, i_e= nxn-1;
       int j_s=1, j_e= nyn-1;
       int k_s=1, k_e= nzn-1;

       int BBX= RFx*4;
       int BBY= RFy*4;
       int BBZ= RFz*4;

       // Exth
       for (int i = i_s; i < i_e; i++)
	 for (int j = j_s; j < j_e; j++)
	   for (int k = k_s; k < k_e; k++){

	     if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){

	       value= 0;

	     }
	     else{
	       if (icount % 2 == 1) {
		 value = 0.;
	       }
	       else {
		 value = 0.5;
	       }

	     }
	     alpha = (1.0 - value) / 6;

	     temp[i][j][k] = value * Ex[i][j][k] + alpha * (Ex[i - 1][j][k] + Ex[i + 1][j][k] + Ex[i][j - 1][k] + Ex[i][j + 1][k] + Ex[i][j][k - 1] + Ex[i][j][k + 1]);}

       for (int i = 1; i < nxn - 1; i++)
	 for (int j = 1; j < nyn - 1; j++)
	   for (int k = 1; k < nzn - 1; k++)
	     Ex[i][j][k] = temp[i][j][k];
       // Eyth
       for (int i = 1; i < nxn - 1; i++)
	 for (int j = 1; j < nyn - 1; j++)
	   for (int k = 1; k < nzn - 1; k++){

	     if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){


	       value= 0;

	     }
	     else{
	       if (icount % 2 == 1) {
		 value = 0.;
	       }
	       else {
		 value = 0.5;
	       }

	     }
	     alpha = (1.0 - value) / 6;

	     temp[i][j][k] = value * Ey[i][j][k] + alpha * (Ey[i - 1][j][k] + Ey[i + 1][j][k] + Ey[i][j - 1][k] + Ey[i][j + 1][k] + Ey[i][j][k - 1] + Ey[i][j][k + 1]);
	   }

       for (int i = 1; i < nxn - 1; i++)
	 for (int j = 1; j < nyn - 1; j++)
	   for (int k = 1; k < nzn - 1; k++){

	     Ey[i][j][k] = temp[i][j][k];
	   }
       // Ezth
       for (int i = 1; i < nxn - 1; i++)
	 for (int j = 1; j < nyn - 1; j++)
	   for (int k = 1; k < nzn - 1; k++){

	     if ( ExtraSmoothing and  ((i < BBX and vct->getXleft_neighbor() == MPI_PROC_NULL) or (i>nxn-BBX and vct->getXright_neighbor() == MPI_PROC_NULL) or (j < BBY and vct->getYleft_neighbor() == MPI_PROC_NULL) or (j>nyn-BBY and vct->getYright_neighbor() == MPI_PROC_NULL) or (k < BBZ and vct->getZleft_neighbor() == MPI_PROC_NULL) or (k>nzn-BBZ and vct->getZright_neighbor() == MPI_PROC_NULL))){


	       value= 0;

	     }
	     else{
	       if (icount % 2 == 1) {
		 value = 0.;
	       }
	       else {
		 value = 0.5;
	       }

	     }
	     alpha = (1.0 - value) / 6;

	     temp[i][j][k] = value * Ez[i][j][k] + alpha * (Ez[i - 1][j][k] + Ez[i + 1][j][k] + Ez[i][j - 1][k] + Ez[i][j + 1][k] + Ez[i][j][k - 1] + Ez[i][j][k + 1]);
	   }
       for (int i = 1; i < nxn - 1; i++)
	 for (int j = 1; j < nyn - 1; j++)
	   for (int k = 1; k < nzn - 1; k++)
	     Ez[i][j][k] = temp[i][j][k];


       delArr3(temp, nxn, nyn);
     }
   }
 }

 /**/


 /* SPECIES: Interpolation smoothing TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector type = 1 --> node based vector */
 void EMfields3D::smooth(double value, double ****vector, int is, int type, Grid * grid, VirtualTopology3D * vct) {
   cout << "Smoothing for Species not implemented in 3D" << endl;
 }

 /*! fix the B boundary when running gem */
 void EMfields3D::fixBgem(Grid * grid, VirtualTopology3D * vct) {
   if (vct->getYright_neighbor() == MPI_PROC_NULL) {
     for (int i = 0; i < nxc; i++)
       for (int k = 0; k < nzc; k++) {
	 Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
	 Bxc[i][nyc - 2][k] = Bxc[i][nyc - 1][k];
	 Bxc[i][nyc - 3][k] = Bxc[i][nyc - 1][k];
	 Byc[i][nyc - 1][k] = B0y;
	 Bzc[i][nyc - 1][k] = B0z;
	 Bzc[i][nyc - 2][k] = B0z;
	 Bzc[i][nyc - 3][k] = B0z;
       }
   }
   if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
     for (int i = 0; i < nxc; i++)
       for (int k = 0; k < nzc; k++) {
	 Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
	 Bxc[i][1][k] = Bxc[i][0][k];
	 Bxc[i][2][k] = Bxc[i][0][k];
	 Byc[i][0][k] = B0y;
	 Bzc[i][0][k] = B0z;
	 Bzc[i][1][k] = B0z;
	 Bzc[i][2][k] = B0z;
       }
   }

 }

 /*! fix the B boundary when running forcefree */
 void EMfields3D::fixBforcefree(Grid * grid, VirtualTopology3D * vct) {
   if (vct->getYright_neighbor() == MPI_PROC_NULL) {
     for (int i = 0; i < nxc; i++)
       for (int k = 0; k < nzc; k++) {
	 Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
	 Byc[i][nyc - 1][k] = B0y;
	 Bzc[i][nyc - 1][k] = B0z / cosh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);;
	 Bzc[i][nyc - 2][k] = B0z / cosh((grid->getYC(i, nyc - 2, k) - Ly / 2) / delta);;
	 Bzc[i][nyc - 3][k] = B0z / cosh((grid->getYC(i, nyc - 3, k) - Ly / 2) / delta);
       }
   }
   if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
     for (int i = 0; i < nxc; i++)
       for (int k = 0; k < nzc; k++) {
	 Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
	 Byc[i][0][k] = B0y;
	 Bzc[i][0][k] = B0z / cosh((grid->getYC(i, 0, k) - Ly / 2) / delta);
	 Bzc[i][1][k] = B0z / cosh((grid->getYC(i, 1, k) - Ly / 2) / delta);
	 Bzc[i][2][k] = B0z / cosh((grid->getYC(i, 2, k) - Ly / 2) / delta);
       }
   }

 }


 /*! adjust densities on boundaries that are not periodic */
 void EMfields3D::adjustNonPeriodicDensities(int is, int bcPfaceXright, int bcPfaceXleft, int bcPfaceYright, int bcPfaceYleft, int bcPfaceZright, int bcPfaceZleft, VirtualTopology3D * vct) {
   /* this function is needed only if i do not have particles in GC; 
      with mlmd, bcPface... <0 they are there, hence skip it */
   if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceXleft <0 ))) {
     for (int i = 1; i < nyn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][1][i][k] += rhons[is][1][i][k];
	 Jxs  [is][1][i][k] += Jxs  [is][1][i][k];
	 Jys  [is][1][i][k] += Jys  [is][1][i][k];
	 Jzs  [is][1][i][k] += Jzs  [is][1][i][k];
	 pXXsn[is][1][i][k] += pXXsn[is][1][i][k];
	 pXYsn[is][1][i][k] += pXYsn[is][1][i][k];
	 pXZsn[is][1][i][k] += pXZsn[is][1][i][k];
	 pYYsn[is][1][i][k] += pYYsn[is][1][i][k];
	 pYZsn[is][1][i][k] += pYZsn[is][1][i][k];
	 pZZsn[is][1][i][k] += pZZsn[is][1][i][k];
       }
   }
   if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceYleft <0 ))) {
     for (int i = 1; i < nxn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][i][1][k] += rhons[is][i][1][k];
	 Jxs  [is][i][1][k] += Jxs  [is][i][1][k];
	 Jys  [is][i][1][k] += Jys  [is][i][1][k];
	 Jzs  [is][i][1][k] += Jzs  [is][i][1][k];
	 pXXsn[is][i][1][k] += pXXsn[is][i][1][k];
	 pXYsn[is][i][1][k] += pXYsn[is][i][1][k];
	 pXZsn[is][i][1][k] += pXZsn[is][i][1][k];
	 pYYsn[is][i][1][k] += pYYsn[is][i][1][k];
	 pYZsn[is][i][1][k] += pYZsn[is][i][1][k];
	 pZZsn[is][i][1][k] += pZZsn[is][i][1][k];
       }
   }
   if (vct->getZleft_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceZleft <0 ))) {
     for (int i = 1; i < nxn - 1; i++)
       for (int j = 1; j < nyn - 1; j++) {
	 rhons[is][i][j][1] += rhons[is][i][j][1];
	 Jxs  [is][i][j][1] += Jxs  [is][i][j][1];
	 Jys  [is][i][j][1] += Jys  [is][i][j][1];
	 Jzs  [is][i][j][1] += Jzs  [is][i][j][1];
	 pXXsn[is][i][j][1] += pXXsn[is][i][j][1];
	 pXYsn[is][i][j][1] += pXYsn[is][i][j][1];
	 pXZsn[is][i][j][1] += pXZsn[is][i][j][1];
	 pYYsn[is][i][j][1] += pYYsn[is][i][j][1];
	 pYZsn[is][i][j][1] += pYZsn[is][i][j][1];
	 pZZsn[is][i][j][1] += pZZsn[is][i][j][1];
       }
   }
   if (vct->getXright_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceXright <0 ))) {

     for (int i = 1; i < nyn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][nxn - 2][i][k] += rhons[is][nxn - 2][i][k];
	 Jxs  [is][nxn - 2][i][k] += Jxs  [is][nxn - 2][i][k];
	 Jys  [is][nxn - 2][i][k] += Jys  [is][nxn - 2][i][k];
	 Jzs  [is][nxn - 2][i][k] += Jzs  [is][nxn - 2][i][k];
	 pXXsn[is][nxn - 2][i][k] += pXXsn[is][nxn - 2][i][k];
	 pXYsn[is][nxn - 2][i][k] += pXYsn[is][nxn - 2][i][k];
	 pXZsn[is][nxn - 2][i][k] += pXZsn[is][nxn - 2][i][k];
	 pYYsn[is][nxn - 2][i][k] += pYYsn[is][nxn - 2][i][k];
	 pYZsn[is][nxn - 2][i][k] += pYZsn[is][nxn - 2][i][k];
	 pZZsn[is][nxn - 2][i][k] += pZZsn[is][nxn - 2][i][k];
       }
   }
   if (vct->getYright_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceYright <0 ))) {

     for (int i = 1; i < nxn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][i][nyn - 2][k] += rhons[is][i][nyn - 2][k];
	 Jxs  [is][i][nyn - 2][k] += Jxs  [is][i][nyn - 2][k];
	 Jys  [is][i][nyn - 2][k] += Jys  [is][i][nyn - 2][k];
	 Jzs  [is][i][nyn - 2][k] += Jzs  [is][i][nyn - 2][k];
	 pXXsn[is][i][nyn - 2][k] += pXXsn[is][i][nyn - 2][k];
	 pXYsn[is][i][nyn - 2][k] += pXYsn[is][i][nyn - 2][k];
	 pXZsn[is][i][nyn - 2][k] += pXZsn[is][i][nyn - 2][k];
	 pYYsn[is][i][nyn - 2][k] += pYYsn[is][i][nyn - 2][k];
	 pYZsn[is][i][nyn - 2][k] += pYZsn[is][i][nyn - 2][k];
	 pZZsn[is][i][nyn - 2][k] += pZZsn[is][i][nyn - 2][k];
       }
   }
   if (vct->getZright_neighbor_P() == MPI_PROC_NULL and (! (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceZright <0 ))) {

     for (int i = 1; i < nxn - 1; i++)
       for (int j = 1; j < nyn - 1; j++) {
	 rhons[is][i][j][nzn - 2] += rhons[is][i][j][nzn - 2];
	 Jxs  [is][i][j][nzn - 2] += Jxs  [is][i][j][nzn - 2];
	 Jys  [is][i][j][nzn - 2] += Jys  [is][i][j][nzn - 2];
	 Jzs  [is][i][j][nzn - 2] += Jzs  [is][i][j][nzn - 2];
	 pXXsn[is][i][j][nzn - 2] += pXXsn[is][i][j][nzn - 2];
	 pXYsn[is][i][j][nzn - 2] += pXYsn[is][i][j][nzn - 2];
	 pXZsn[is][i][j][nzn - 2] += pXZsn[is][i][j][nzn - 2];
	 pYYsn[is][i][j][nzn - 2] += pYYsn[is][i][j][nzn - 2];
	 pYZsn[is][i][j][nzn - 2] += pYZsn[is][i][j][nzn - 2];
	 pZZsn[is][i][j][nzn - 2] += pZZsn[is][i][j][nzn - 2];
       }
   }


   /*cout << "APPLYING MOMENT PATCH, PROBABLY TO REMOVE: " << endl;
     // as expected, this does not influence results 

   // careful: I have changed the ! in the condition  and the index 

   if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceXleft <0 ))) {
     for (int i = 1; i < nyn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][0][i][k] += rhons[is][0][i][k];
	 Jxs  [is][0][i][k] += Jxs  [is][0][i][k];
	 Jys  [is][0][i][k] += Jys  [is][0][i][k];
	 Jzs  [is][0][i][k] += Jzs  [is][0][i][k];
	 pXXsn[is][0][i][k] += pXXsn[is][0][i][k];
	 pXYsn[is][0][i][k] += pXYsn[is][0][i][k];
	 pXZsn[is][0][i][k] += pXZsn[is][0][i][k];
	 pYYsn[is][0][i][k] += pYYsn[is][0][i][k];
	 pYZsn[is][0][i][k] += pYZsn[is][0][i][k];
	 pZZsn[is][0][i][k] += pZZsn[is][0][i][k];
       }
   }
   if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceYleft <0 ))) {
     for (int i = 1; i < nxn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][i][0][k] += rhons[is][i][0][k];
	 Jxs  [is][i][0][k] += Jxs  [is][i][0][k];
	 Jys  [is][i][0][k] += Jys  [is][i][0][k];
	 Jzs  [is][i][0][k] += Jzs  [is][i][0][k];
	 pXXsn[is][i][0][k] += pXXsn[is][i][0][k];
	 pXYsn[is][i][0][k] += pXYsn[is][i][0][k];
	 pXZsn[is][i][0][k] += pXZsn[is][i][0][k];
	 pYYsn[is][i][0][k] += pYYsn[is][i][0][k];
	 pYZsn[is][i][0][k] += pYZsn[is][i][0][k];
	 pZZsn[is][i][0][k] += pZZsn[is][i][0][k];
       }
   }
   if (vct->getZleft_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceZleft <0 ))) {
     for (int i = 1; i < nxn - 1; i++)
       for (int j = 1; j < nyn - 1; j++) {
	 rhons[is][i][j][0] += rhons[is][i][j][0];
	 Jxs  [is][i][j][0] += Jxs  [is][i][j][0];
	 Jys  [is][i][j][0] += Jys  [is][i][j][0];
	 Jzs  [is][i][j][0] += Jzs  [is][i][j][0];
	 pXXsn[is][i][j][0] += pXXsn[is][i][j][0];
	 pXYsn[is][i][j][0] += pXYsn[is][i][j][0];
	 pXZsn[is][i][j][0] += pXZsn[is][i][j][0];
	 pYYsn[is][i][j][0] += pYYsn[is][i][j][0];
	 pYZsn[is][i][j][0] += pYZsn[is][i][j][0];
	 pZZsn[is][i][j][0] += pZZsn[is][i][j][0];
       }
   }
   if (vct->getXright_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceXright <0 ))) {

     for (int i = 1; i < nyn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][nxn - 1][i][k] += rhons[is][nxn - 1][i][k];
	 Jxs  [is][nxn - 1][i][k] += Jxs  [is][nxn - 1][i][k];
	 Jys  [is][nxn - 1][i][k] += Jys  [is][nxn - 1][i][k];
	 Jzs  [is][nxn - 1][i][k] += Jzs  [is][nxn - 1][i][k];
	 pXXsn[is][nxn - 1][i][k] += pXXsn[is][nxn - 1][i][k];
	 pXYsn[is][nxn - 1][i][k] += pXYsn[is][nxn - 1][i][k];
	 pXZsn[is][nxn - 1][i][k] += pXZsn[is][nxn - 1][i][k];
	 pYYsn[is][nxn - 1][i][k] += pYYsn[is][nxn - 1][i][k];
	 pYZsn[is][nxn - 1][i][k] += pYZsn[is][nxn - 1][i][k];
	 pZZsn[is][nxn - 1][i][k] += pZZsn[is][nxn - 1][i][k];
       }
   }
   if (vct->getYright_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceYright <0 ))) {

     for (int i = 1; i < nxn - 1; i++)
       for (int k = 1; k < nzn - 1; k++) {
	 rhons[is][i][nyn - 1][k] += rhons[is][i][nyn - 1][k];
	 Jxs  [is][i][nyn - 1][k] += Jxs  [is][i][nyn - 1][k];
	 Jys  [is][i][nyn - 1][k] += Jys  [is][i][nyn - 1][k];
	 Jzs  [is][i][nyn - 1][k] += Jzs  [is][i][nyn - 1][k];
	 pXXsn[is][i][nyn - 1][k] += pXXsn[is][i][nyn - 1][k];
	 pXYsn[is][i][nyn - 1][k] += pXYsn[is][i][nyn - 1][k];
	 pXZsn[is][i][nyn - 1][k] += pXZsn[is][i][nyn - 1][k];
	 pYYsn[is][i][nyn - 1][k] += pYYsn[is][i][nyn - 1][k];
	 pYZsn[is][i][nyn - 1][k] += pYZsn[is][i][nyn - 1][k];
	 pZZsn[is][i][nyn - 1][k] += pZZsn[is][i][nyn - 1][k];
       }
   }
   if (vct->getZright_neighbor_P() == MPI_PROC_NULL and ( (vct->getCommToParent_P(is)!= MPI_COMM_NULL and ParticleREPOPULATION and bcPfaceZright <0 ))) {

     for (int i = 1; i < nxn - 1; i++)
       for (int j = 1; j < nyn - 1; j++) {
	 rhons[is][i][j][nzn - 1] += rhons[is][i][j][nzn - 1];
	 Jxs  [is][i][j][nzn - 1] += Jxs  [is][i][j][nzn - 1];
	 Jys  [is][i][j][nzn - 1] += Jys  [is][i][j][nzn - 1];
	 Jzs  [is][i][j][nzn - 1] += Jzs  [is][i][j][nzn - 1];
	 pXXsn[is][i][j][nzn - 1] += pXXsn[is][i][j][nzn - 1];
	 pXYsn[is][i][j][nzn - 1] += pXYsn[is][i][j][nzn - 1];
	 pXZsn[is][i][j][nzn - 1] += pXZsn[is][i][j][nzn - 1];
	 pYYsn[is][i][j][nzn - 1] += pYYsn[is][i][j][nzn - 1];
	 pYZsn[is][i][j][nzn - 1] += pYZsn[is][i][j][nzn - 1];
	 pZZsn[is][i][j][nzn - 1] += pZZsn[is][i][j][nzn - 1];
       }
       }*/

 }

 void EMfields3D::ConstantChargeOpenBCv2(Grid * grid, VirtualTopology3D * vct) {

   int nx = grid->getNXN();
   int ny = grid->getNYN();
   int nz = grid->getNZN();

   for (int is = 0; is < ns; is++) {

     if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
       for (int j=0; j < ny;j++)
	 for (int k=0; k < nz;k++){
	   rhons[is][0][j][k] = rhons[is][4][j][k];
	   rhons[is][1][j][k] = rhons[is][4][j][k];
	   rhons[is][2][j][k] = rhons[is][4][j][k];
	   rhons[is][3][j][k] = rhons[is][4][j][k];
	 }
     }

     if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
       for (int j=0; j < ny;j++)
	 for (int k=0; k < nz;k++){
	   rhons[is][nx-4][j][k] = rhons[is][nx-5][j][k];
	   rhons[is][nx-3][j][k] = rhons[is][nx-5][j][k];
	   rhons[is][nx-2][j][k] = rhons[is][nx-5][j][k];
	   rhons[is][nx-1][j][k] = rhons[is][nx-5][j][k];
	 }
     }

     if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
       for (int i=0; i < nx;i++)
	 for (int k=0; k < nz;k++){
	   rhons[is][i][0][k] = rhons[is][i][4][k];
	   rhons[is][i][1][k] = rhons[is][i][4][k];
	   rhons[is][i][2][k] = rhons[is][i][4][k];
	   rhons[is][i][3][k] = rhons[is][i][4][k];
	 }
     }

     if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
       for (int i=0; i < nx;i++)
	 for (int k=0; k < nz;k++){
	   rhons[is][i][ny-4][k] = rhons[is][i][ny-5][k];
	   rhons[is][i][ny-3][k] = rhons[is][i][ny-5][k];
	   rhons[is][i][ny-2][k] = rhons[is][i][ny-5][k];
	   rhons[is][i][ny-1][k] = rhons[is][i][ny-5][k];
	 }
     }

     if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
       for (int i=0; i < nx;i++)
	 for (int j=0; j < ny;j++){
	   rhons[is][i][j][0] = rhons[is][i][j][4];
	   rhons[is][i][j][1] = rhons[is][i][j][4];
	   rhons[is][i][j][2] = rhons[is][i][j][4];
	   rhons[is][i][j][3] = rhons[is][i][j][4];
	 }
     }


     if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
       for (int i=0; i < nx;i++)
	 for (int j=0; j < ny;j++){
	   rhons[is][i][j][nz-4] = rhons[is][i][j][nz-5];
	   rhons[is][i][j][nz-3] = rhons[is][i][j][nz-5];
	   rhons[is][i][j][nz-2] = rhons[is][i][j][nz-5];
	   rhons[is][i][j][nz-1] = rhons[is][i][j][nz-5];
	 }
     }
   }

 }

 void EMfields3D::ConstantChargeOpenBC(Grid * grid, VirtualTopology3D * vct) {

   int nx = grid->getNXN();
   int ny = grid->getNYN();
   int nz = grid->getNZN();

   for (int is = 0; is < ns; is++) {

     double ff = (qom[is] / fabs(qom[is]));

     if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcEMfaceXleft ==2)) {
       for (int j=0; j < ny;j++)
	 for (int k=0; k < nz;k++){
	   rhons[is][0][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][1][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][2][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][3][j][k] = ff * rhoINIT[is] / FourPI;
	 }
     }

     if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcEMfaceXright ==2)) {
       for (int j=0; j < ny;j++)
	 for (int k=0; k < nz;k++){
	   rhons[is][nx-4][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][nx-3][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][nx-2][j][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][nx-1][j][k] = ff * rhoINIT[is] / FourPI;
	 }
     }

     if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcEMfaceYleft ==2))  {
       for (int i=0; i < nx;i++)
	 for (int k=0; k < nz;k++){
	   rhons[is][i][0][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][1][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][2][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][3][k] = ff * rhoINIT[is] / FourPI;
	 }
     }

     if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcEMfaceYright ==2))  {
       for (int i=0; i < nx;i++)
	 for (int k=0; k < nz;k++){
	   rhons[is][i][ny-4][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][ny-3][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][ny-2][k] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][ny-1][k] = ff * rhoINIT[is] / FourPI;
	 }
     }

     if(vct->getZleft_neighbor()==MPI_PROC_NULL && (bcEMfaceZleft ==2))  {
       for (int i=0; i < nx;i++)
	 for (int j=0; j < ny;j++){
	   rhons[is][i][j][0] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][1] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][2] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][3] = ff * rhoINIT[is] / FourPI;
	 }
     }


     if(vct->getZright_neighbor()==MPI_PROC_NULL && (bcEMfaceZright ==2))  {
       for (int i=0; i < nx;i++)
	 for (int j=0; j < ny;j++){
	   rhons[is][i][j][nz-4] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][nz-3] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][nz-2] = ff * rhoINIT[is] / FourPI;
	   rhons[is][i][j][nz-1] = ff * rhoINIT[is] / FourPI;
	 }
     }
   }

 }

 void EMfields3D::ConstantChargePlanet(Grid * grid, VirtualTopology3D * vct, double R, double x_center, double y_center, double z_center) {

   double xd;
   double yd;
   double zd;

   for (int is = 0; is < ns; is++) {
     for (int i = 1; i < nxn; i++) {
       for (int j = 1; j < nyn; j++) {
	 for (int k = 1; k < nzn; k++) {

	   xd = fabs(grid->getXN(i,j,k) - x_center) - dx;
	   yd = fabs(grid->getYN(i,j,k) - y_center) - dy;
	   zd = fabs(grid->getZN(i,j,k) - z_center) - dz;

	   if ((xd*xd+yd*yd+zd*zd) < R*R) {
	     rhons[is][i][j][k] = (qom[is] / fabs(qom[is])) * rhoINJECT[is] / FourPI;
	   }

	 }
       }
     }
   }

 }

 /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
void EMfields3D::calculateB(Grid * grid, VirtualTopology3D * vct, Collective *col, int cycle) {

  /* left here from a test
  if (numGrid==1){
    grid->interpC2N(Bxn, Bxc);
    grid->interpC2N(Byn, Byc);
    grid->interpC2N(Bzn, Bzc);

    cout << "Grid 1: i am not calculating B; i am interpolating what i received";

    return;
    }else {return;}
  */
  if (vct->getCartesian_rank() == 0)
     //cout << "*** B CALCULATION ***" << endl;
     cout << "*** G" << numGrid <<": B CALCULATION ***" << endl;

   // calculate the curl of Eth
   if (Case != "TestFix3B"){
     grid->curlN2C(tempXC, tempYC, tempZC, Exth, Eyth, Ezth);
     
     if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC and (MLMD_fixBCenters or MLMD_InterpolateOldBCell)){       
       if (vct->getCartesian_rank() == 0)
	 cout << "I am NOT updating the first active cell in B" << endl;
       
       for (int i=0; i< nxc; i++)
	 for (int j=0; j< nyc; j++)
	   for (int k=0; k< nzc; k++){
	     bool CASEX1= (vct->getXleft_neighbor() == MPI_PROC_NULL and (i==0 or i==1));
	     bool CASEX2= (vct->getXright_neighbor() == MPI_PROC_NULL and (i==nxc-1 or i==nxc-2));
	     bool CASEY1= (vct->getYleft_neighbor() == MPI_PROC_NULL and (j==0 or j==1));
	     bool CASEY2= (vct->getYright_neighbor() == MPI_PROC_NULL and (j==nyc-1 or j==nyc-2));
	     bool CASEZ1= (vct->getZleft_neighbor() == MPI_PROC_NULL and (k==0 or k==1));
	     bool CASEZ2= (vct->getZright_neighbor() == MPI_PROC_NULL and (k==nzc-1 or k==nzc-2));

	     if (CASEX1 or CASEX2 or CASEY1 or CASEY2 or CASEZ1 or CASEZ2){
	       tempXC[i][j][k]=0.0;
	       tempYC[i][j][k]=0.0;
	       tempZC[i][j][k]=0.0;
	     }
	   }
     } // end of non-updating curl E at the boundaries


     // update the magnetic field
     addscale(-c * dt, 1, Bxc, tempXC, nxc, nyc, nzc);
     addscale(-c * dt, 1, Byc, tempYC, nxc, nyc, nzc);
     addscale(-c * dt, 1, Bzc, tempZC, nxc, nyc, nzc);
     
     // communicate ghost 
     /*communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
     communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
     communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);*/

     communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct);
     communicateCenterBC(nxc, nyc, nzc, Byc, 2, 2, 2, 2, 2, 2, vct);
     communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct);
     
     if (! (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC)){
       
       if (Case=="ForceFree") fixBforcefree(grid,vct);
       if (Case=="GEM")       fixBgem(grid, vct);
       if (Case=="GEMnoPert") fixBgem(grid, vct);

       // OpenBC:
       BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc,grid,vct);
     }

   } 

   grid->interpC2N(Bxn, Bxc);
   grid->interpC2N(Byn, Byc);
   grid->interpC2N(Bzn, Bzc);

   //cout << "Before communicateNodeBC"<<endl;
   communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
   communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
   //cout << "End communicateNodeBC"<<endl;
   
   /* do not put B BC here, they would be at old state */ 

   /* HERE I PRINT E N+1, B N+1--
      IF I PRINT IN ENDECALC, IT WOULD BE B N */
   /*double maxintE=0.0; double maxintB=0.0;
   double intE, intB;
   for (int i=1; i< nxn-1; i++)
     for (int j=1; j< nyn-1; j++)
       for (int k=1; k< nzn-1; k++){
	 intE= Ex[i][j][k]*Ex[i][j][k]+ Ey[i][j][k]*Ey[i][j][k]+ Ez[i][j][k]*Ez[i][j][k];

	 if (intE> maxintE) maxintE=intE;
       }

   for (int i=1; i< nxc-1; i++)
     for (int j=1; j< nyc-1; j++)
       for (int k=1; k< nzc-1; k++){
	 intB= Bxc[i][j][k]*Bxc[i][j][k]+ Byc[i][j][k]*Byc[i][j][k]+ Bzc[i][j][k]*Bzc[i][j][k];

	 if (intB> maxintB) maxintB=intB;
       }

       cout << "Grid " << numGrid<< ", Cycle " << cycle <<": maxintE: " << maxintE << ", maxintB: " << maxintB <<endl */
   return;

 }

 void EMfields3D::fixBghostCells_Left(int Dir, int NCells){

   int i_s= 0, i_e= nxc-1;
   int j_s= 0, j_e= nyc-1;
   int k_s= 0, k_e= nzc-1;

   if (Dir==0) i_e= NCells;
   if (Dir==1) j_e= NCells;
   if (Dir==2) k_e= NCells;


   for (int i=i_s; i< i_e; i++)
     for (int j=j_s; j<j_e; j++)
       for(int k=k_s; k<k_e; k++){
	 Bxc[i][j][k] = .125 * (Bxn[i][j][k] + Bxn[i + 1][j][k] + Bxn[i][j + 1][k] + Bxn[i][j][k + 1] + Bxn[i + 1][j + 1][k] + Bxn[i + 1][j][k + 1] + Bxn[i][j + 1][k + 1] + Bxn[i + 1][j + 1][k + 1]);
	 Byc[i][j][k] = .125 * (Byn[i][j][k] + Byn[i + 1][j][k] + Byn[i][j + 1][k] + Byn[i][j][k + 1] + Byn[i + 1][j + 1][k] + Byn[i + 1][j][k + 1] + Byn[i][j + 1][k + 1] + Byn[i + 1][j + 1][k + 1]);
	 Bzc[i][j][k] = .125 * (Bzn[i][j][k] + Bzn[i + 1][j][k] + Bzn[i][j + 1][k] + Bzn[i][j][k + 1] + Bzn[i + 1][j + 1][k] + Bzn[i + 1][j][k + 1] + Bzn[i][j + 1][k + 1] + Bzn[i + 1][j + 1][k + 1]);
       }
 }


 void EMfields3D::fixBghostCells_Right(int Dir, int NCells){

   int i_s= 0, i_e= nxc-1;
   int j_s= 0, j_e= nyc-1;
   int k_s= 0, k_e= nzc-1;

   if (Dir==0) i_s= nxc-1-NCells;
   if (Dir==1) j_s= nyc-1-NCells;
   if (Dir==2) k_s= nzc-1-NCells;


   for (int i=i_s; i< i_e; i++)
     for (int j=j_s; j<j_e; j++)
       for(int k=k_s; k<k_e; k++){
	 Bxc[i][j][k] = .125 * (Bxn[i][j][k] + Bxn[i + 1][j][k] + Bxn[i][j + 1][k] + Bxn[i][j][k + 1] + Bxn[i + 1][j + 1][k] + Bxn[i + 1][j][k + 1] + Bxn[i][j + 1][k + 1] + Bxn[i + 1][j + 1][k + 1]);
	 Byc[i][j][k] = .125 * (Byn[i][j][k] + Byn[i + 1][j][k] + Byn[i][j + 1][k] + Byn[i][j][k + 1] + Byn[i + 1][j + 1][k] + Byn[i + 1][j][k + 1] + Byn[i][j + 1][k + 1] + Byn[i + 1][j + 1][k + 1]);
	 Bzc[i][j][k] = .125 * (Bzn[i][j][k] + Bzn[i + 1][j][k] + Bzn[i][j + 1][k] + Bzn[i][j][k + 1] + Bzn[i + 1][j + 1][k] + Bzn[i + 1][j][k + 1] + Bzn[i][j + 1][k + 1] + Bzn[i + 1][j + 1][k + 1]);
       }
 }

 /*! initialize EM field with transverse electric waves 1D and rotate anticlockwise (theta degrees) */
 void EMfields3D::initEM_rotate(VirtualTopology3D * vct, Grid * grid, Collective *col, double B, double theta) {
   // initialize E and rhos on nodes
   for (int i = 0; i < nxn; i++)
     for (int j = 0; j < nyn; j++) {
       Ex[i][j][0] = 0.0;
       Ey[i][j][0] = 0.0;
       Ez[i][j][0] = 0.0;
       Bxn[i][j][0] = B * cos(theta * M_PI / 180);
       Byn[i][j][0] = B * sin(theta * M_PI / 180);
       Bzn[i][j][0] = 0.0;
       rhons[0][i][j][0] = 0.07957747154595; // electrons: species is now first index
       rhons[1][i][j][0] = 0.07957747154595; // protons: species is now first index
     }
   // initialize B on centers
   grid->interpN2C(Bxc, Bxn);
   grid->interpN2C(Byc, Byn);
   grid->interpN2C(Bzc, Bzn);


   for (int is = 0; is < ns; is++)
     grid->interpN2C_GC(rhocs, is, rhons);

 }
 /* initialize a light wave -MLMD ready */
 void EMfields3D::initLightWave(VirtualTopology3D * vct, Grid * grid, Collective *col){

   if (numGrid >0 and col->getMLMD_InitialInterpolation()) {
     int rr= vct->getCartesian_rank();

     if (rr==0){
       cout << "I will be initialising the RG by interpolation" << endl;
     }

     return; // i want to initialise RG differently 
   }

   // to use with periodic BC

   if (restart1 == 0) {
     // initialize
     if (vct->getCartesian_rank() == 0) {
       cout << "------------------------------------------" << endl;
       cout << "Initialize Light Wave in MLMD system" << endl;
       cout << "------------------------------------------" << endl;
     }

     double PI=FourPI/4.0;
     // density does not need to be initialised   
     double intensity= 0.1;
     double KX= 2*PI/(col->getLx_mlmd(0)); // may need a + dx
     double KY= 2*PI/(col->getLy_mlmd(0)); // may need a + dy
     double K= sqrt(KX*KX+ KY*KY);
     double startingT= 0.25/K;

     double globalX, globalY, globalZ;

     /*cout << "Grid " << numGrid << " Ox "<< col->getOx_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oy "<< col->getOy_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oz "<< col->getOz_SW(numGrid) <<endl;*/

     cout << "iPic: intensity: " << intensity << " KX "<< KX << " KY " << KY << " K " << K << " startingT " << startingT << endl;    

     for (int i = 0; i < nxn; i++)
       for (int j = 0; j < nyn; j++)
	 for (int k = 0; k < nzn; k++) {
	   // initialize the density for species

	   globalX= grid->getXN(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYN(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZN(i,j,k)+ col->getOz_SW(numGrid) ;
	   // electric field
	   Ex[i][j][k] = intensity*KY*cos(KX*globalX)*cos(KY*globalY)*sin(K*startingT)/KX;
	   Ey[i][j][k] = intensity*sin(KX*globalX)*sin(KY*globalY) *sin(K*startingT); 
	   Ez[i][j][k] = 0.0;
	   // Magnetic field
	   Bxn[i][j][k] = 0.0;
	   Byn[i][j][k] = 0.0;
	   
	   Bzn[i][j][k] = intensity*K*cos(KX*globalX)*sin(KY*globalY) * cos(K*startingT)/KX;
	   /* left here from a test
	   Bxn[i][j][k]=globalX;
	   Byn[i][j][k]=globalY;
	     
	   if (numGrid==1){
	     cout << "SONO QUI " << endl;
	     Bzn[i][j][k]=0.0;
	     Byn[i][j][k]=0.0;
	     Bxn[i][j][k]=0.0;
	     } */
	   
	 }


     // initialize B on centers
     for (int i = 0; i < nxc; i++)
       for (int j = 0; j < nyc; j++)
	 for (int k = 0; k < nzc; k++) {

	   globalX= grid->getXC(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYC(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZC(i,j,k)+ col->getOz_SW(numGrid) ;

	   Bxc[i][j][k] = 0.0;
	   Byc[i][j][k] = 0.0;
	   Bzc[i][j][k] = intensity*K*cos(KX*globalX)*sin(KY*globalY)*cos(K*startingT)/KX; 

	   /* QUA left here from a test
	   Bxc[i][j][k]=globalX;
	   Byc[i][j][k]=globalY;
	     
	   if (numGrid==1){
	     cout << "SONO QUI 2" << endl;
	     Bzc[i][j][k]=0.0;
	     Byc[i][j][k]=0.0;
	     Bxc[i][j][k]=0.0;
	     }*/
	   


	 }

     // initialize B on centers
     /*grid->interpN2C_GC(Bxc, Bxn);
     grid->interpN2C_GC(Byc, Byn);
     grid->interpN2C_GC(Bzc, Bzn);*/


     for (int is = 0; is < ns; is++)
       grid->interpN2C_GC(rhocs, is, rhons);

   }
   else {
     init(vct, grid, col);            // use the fields from restart file
   }

 }

void EMfields3D::initTestIntProj(VirtualTopology3D * vct, Grid * grid, Collective *col){

   if (numGrid >0 and col->getMLMD_InitialInterpolation()) {
     int rr= vct->getCartesian_rank();

     if (rr==0){
       cout << "I will be initialising the RG by interpolation" << endl;
     }

     return; // i want to initialise RG differently 
   }

   // to use with periodic BC

   if (restart1 == 0) {
     // initialize
     if (vct->getCartesian_rank() == 0) {
       cout << "------------------------------------------" << endl;
       cout << "initTestIntProj in MLMD system" << endl;
       cout << "------------------------------------------" << endl;
     }

     double PI=FourPI/4.0;
     // density does not need to be initialised   
     double intensity= 0.1;
     double KX= 2*PI/(col->getLx_mlmd(0)); // may need a + dx
     double KY= 2*PI/(col->getLy_mlmd(0)); // may need a + dy
     double K= sqrt(KX*KX+ KY*KY);
     double startingT= 0.25/K;

     double globalX, globalY, globalZ;

     /*cout << "Grid " << numGrid << " Ox "<< col->getOx_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oy "<< col->getOy_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oz "<< col->getOz_SW(numGrid) <<endl;*/

     cout << "iPic: intensity: " << intensity << " KX "<< KX << " KY " << KY << " K " << K << " startingT " << startingT << endl;    

     for (int i = 0; i < nxn; i++)
       for (int j = 0; j < nyn; j++)
	 for (int k = 0; k < nzn; k++) {
	   // initialize the density for species

	   globalX= grid->getXN(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYN(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZN(i,j,k)+ col->getOz_SW(numGrid) ;
	   // electric field
	   Ex[i][j][k] = 0.0; //intensity*KY*cos(KX*globalX)*cos(KY*globalY)*sin(K*startingT)/KX;
	   Ey[i][j][k] = 0.0; //intensity*sin(KX*globalX)*sin(KY*globalY) *sin(K*startingT); 
	   Ez[i][j][k] = 0.0;
	   // Magnetic field
	   Bxn[i][j][k] = 0.0;
	   Byn[i][j][k] = 0.0;
	   Bzn[i][j][k] = 0.0; //intensity*K*cos(KX*globalX)*sin(KY*globalY) * cos(K*startingT)/KX;
	 }


     // initialize B on centers
     for (int i = 0; i < nxc; i++)
       for (int j = 0; j < nyc; j++)
	 for (int k = 0; k < nzc; k++) {

	   globalX= grid->getXC(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYC(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZC(i,j,k)+ col->getOz_SW(numGrid) ;

	   Bxc[i][j][k] = 0.0;
	   Byc[i][j][k] = 0.0;
	   Bzc[i][j][k] = 0.0; //intensity*K*cos(KX*globalX)*sin(KY*globalY)*cos(K*startingT)/KX; 

	   }

     // initialize B on centers
     /*grid->interpN2C_GC(Bxc, Bxn);
     grid->interpN2C_GC(Byc, Byn);
     grid->interpN2C_GC(Bzc, Bzn);*/


     for (int is = 0; is < ns; is++)
       grid->interpN2C_GC(rhocs, is, rhons);

     /* i put some b.s. in the RG */
     if (numGrid >0){
       for (int i=0; i< nxn; i++)
	 for (int j=0; j< nyn; j++)
	   for (int k=0; k< nzn; k++){
	     Ex[i][j][k]= 0.0;
	     Ey[i][j][k]= 0.0;
	     Ez[i][j][k]= 0.0;
	     Exth[i][j][k]= 5; //i;
	     Eyth[i][j][k]= 6; //10*i;
	     Ezth[i][j][k]= 7; //100*i;

	     Bxn[i][j][k]= 1000;
	     Byn[i][j][k]= 1000;
	     Bzn[i][j][k]= 1000;
	   }
       for (int i=0; i< nxc; i++)
	 for (int j=0; j< nyc; j++)
	   for (int k=0; k< nzc; k++){
	     Bxc[i][j][k]= 1000;
	     Byc[i][j][k]= 1000;
	     Bzc[i][j][k]= 1000;
	   }

     }


   }
   else {
     init(vct, grid, col);            // use the fields from restart file
   }

}


void EMfields3D::initTestBBoundary(VirtualTopology3D * vct, Grid * grid, Collective *col){

   if (numGrid >0 and col->getMLMD_InitialInterpolation()) {
     int rr= vct->getCartesian_rank();

     if (rr==0){
       cout << "I will be initialising the RG by interpolation" << endl;
     }

     return; // i want to initialise RG differently 
   }

   // to use with periodic BC

   if (restart1 == 0) {
     // initialize
     if (vct->getCartesian_rank() == 0) {
       cout << "------------------------------------------" << endl;
       cout << "Initialize Light Wave in MLMD system" << endl;
       cout << "------------------------------------------" << endl;
     }

     double PI=FourPI/4.0;
     // density does not need to be initialised   
     double intensity= 0.1;
     double KX= 2*PI/(col->getLx_mlmd(0)); // may need a + dx
     double KY= 2*PI/(col->getLy_mlmd(0)); // may need a + dy
     double K= sqrt(KX*KX+ KY*KY);
     double startingT= 0.25/K;

     double globalX, globalY, globalZ;

     /*cout << "Grid " << numGrid << " Ox "<< col->getOx_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oy "<< col->getOy_SW(numGrid) <<endl;
     cout << "Grid " <<numGrid<< " Oz "<< col->getOz_SW(numGrid) <<endl;*/

     for (int i = 0; i < nxn; i++)
       for (int j = 0; j < nyn; j++)
	 for (int k = 0; k < nzn; k++) {
	   // initialize the density for species

	   globalX= grid->getXN(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYN(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZN(i,j,k)+ col->getOz_SW(numGrid) ;
	   // electric field
	   Exth[i][j][k] = intensity*KY*cos(KX*globalX)*cos(KY*globalY)*sin(K*startingT)/KX;
	   Eyth[i][j][k] = intensity*sin(KX*globalX)*sin(KY*globalY) *sin(K*startingT); 
	   Ezth[i][j][k] = 0.0;

	   Ex[i][j][k] = intensity*KY*cos(KX*globalX)*cos(KY*globalY)*sin(K*startingT)/KX;
	   Ey[i][j][k] = intensity*sin(KX*globalX)*sin(KY*globalY) *sin(K*startingT); 
	   Ez[i][j][k] = 0.0;

	   // Magnetic field
	   Bxn[i][j][k] = 0.0;
	   Byn[i][j][k] = 0.0;
	   Bzn[i][j][k] = intensity*K*cos(KX*globalX)*sin(KY*globalY) * cos(K*startingT)/KX;
	 }


     // initialize B on centers
     for (int i = 0; i < nxc; i++)
       for (int j = 0; j < nyc; j++)
	 for (int k = 0; k < nzc; k++) {

	   globalX= grid->getXC(i,j,k)+ col->getOx_SW(numGrid) ;
	   globalY= grid->getYC(i,j,k)+ col->getOy_SW(numGrid) ;
	   globalZ= grid->getZC(i,j,k)+ col->getOz_SW(numGrid) ;

	   Bxc[i][j][k] = 0.0;
	   Byc[i][j][k] = 0.0;
	   Bzc[i][j][k] = intensity*K*cos(KX*globalX)*sin(KY*globalY)*cos(K*startingT)/KX; 

	   }

     // initialize B on centers
     /*grid->interpN2C_GC(Bxc, Bxn);
     grid->interpN2C_GC(Byc, Byn);
     grid->interpN2C_GC(Bzc, Bzn);*/


     for (int is = 0; is < ns; is++)
       grid->interpN2C_GC(rhocs, is, rhons);

   }
   else {
     init(vct, grid, col);            // use the fields from restart file
   }

 }

 /*!Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
 void EMfields3D::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid) {

   double alpha;
   alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

   ne_mod *= alpha;
   ni_mod *= alpha;
   // cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
   for (int i = 0; i < nxn; i++)
     for (int j = 0; j < nyn; j++) {
       rhons[0][i][j][0] += ne_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ne_phase);
       rhons[1][i][j][0] += ni_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ni_phase);
     }

   for (int is = 0; is < ns; is++)
     grid->interpN2C_GC(rhocs, is, rhons);
 }


 /*!Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
 void EMfields3D::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid) {

   double alpha;

   alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

   Ex_mod *= alpha;
   Ey_mod *= alpha;
   Ez_mod *= alpha;
   Bx_mod *= alpha;
   By_mod *= alpha;
   Bz_mod *= alpha;

   for (int i = 0; i < nxn; i++)
     for (int j = 0; j < nyn; j++) {
       Ex[i][j][0] += Ex_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ex_phase);
       Ey[i][j][0] += Ey_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ey_phase);
       Ez[i][j][0] += Ez_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ez_phase);
       Bxn[i][j][0] += Bx_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bx_phase);
       Byn[i][j][0] += By_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + By_phase);
       Bzn[i][j][0] += Bz_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bz_phase);

     }

   // initialize B on centers
   grid->interpN2C(Bxc, Bxn);
   grid->interpN2C(Byc, Byn);
   grid->interpN2C(Bzc, Bzn);


 }


 /*! Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
 void EMfields3D::calculateHatFunctions(Grid * grid, VirtualTopology3D * vct) {
   // smoothing

   smooth(Smooth, rhoc, 0, grid, vct);

   // calculate j hat

   for (int is = 0; is < ns; is++) {
     grid->divSymmTensorN2C(tempXC, tempYC, tempZC, pXXsn, pXYsn, pXZsn, pYYsn, pYZsn, pZZsn, is);

     scale(tempXC, -dt / 2.0, nxc, nyc, nzc);
     scale(tempYC, -dt / 2.0, nxc, nyc, nzc);
     scale(tempZC, -dt / 2.0, nxc, nyc, nzc);
     // communicate before interpolating -- with 2, it copies the same value in the ghost cells
     communicateCenterBC_P(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
     communicateCenterBC_P(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
     communicateCenterBC_P(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);

     grid->interpC2N(tempXN, tempXC);
     grid->interpC2N(tempYN, tempYC);
     grid->interpC2N(tempZN, tempZC);
     sum(tempXN, Jxs, nxn, nyn, nzn, is);
     sum(tempYN, Jys, nxn, nyn, nzn, is);
     sum(tempZN, Jzs, nxn, nyn, nzn, is);
     // PIDOT
     PIdot(Jxh, Jyh, Jzh, tempXN, tempYN, tempZN, is, grid);

   }
   smooth(Smooth, Jxh, 1, grid, vct);
   smooth(Smooth, Jyh, 1, grid, vct);
   smooth(Smooth, Jzh, 1, grid, vct);
   // calculate rho hat = rho - (dt*theta)div(jhat)
  grid->divN2C(tempXC, Jxh, Jyh, Jzh);
  scale(tempXC, -dt * th, nxc, nyc, nzc);
  sum(tempXC, rhoc, nxc, nyc, nzc);
  eq(rhoh, tempXC, nxc, nyc, nzc);
  // communicate rhoh
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct);
}
/*! Image of Poisson Solver */
void EMfields3D::PoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct) {
  // allocate 2 three dimensional service vectors
  double ***temp = newArr3(double, nxc, nyc, nzc);
  double ***im = newArr3(double, nxc, nyc, nzc);
  eqValue(0.0, image, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, nxc, nyc, nzc);
  // calculate the laplacian
  grid->lapC2Cpoisson(im, temp, vct);
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);
  // deallocate temporary array and objects
  delArr3(temp, nxc, nyc);
  delArr3(im, nxc, nyc);
}
void EMfields3D::PoissonImage_2D(double *image, double *vector, Grid * grid, VirtualTopology3D * vct, int nxc, int nyc, int nzc) {

  cout << "PoissonImage_2D: nxc= " << nxc << ", nyc= " << nyc << ", nzc= " << nzc << endl;
  // allocate 2 three dimensional service vectors
  double ***temp = newArr3(double, nxc, nyc, nzc);
  double ***im = newArr3(double, nxc, nyc, nzc);
  eqValue(0.0, image, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, nxc, nyc, nzc);
  // calculate the laplacian
  grid->lapC2Cpoisson(im, temp, vct, nxc, nyc, nzc);

  // I am putting BC here: df=0 at first active , 
  // this correction ok only for X left
  for (int i=0; i< nxc; i++)
    for (int j=0; j< nyc; j++)
      for (int k=0; k< nzc; k++){
	if (j==1 or k==1 or j== nyc-2 or k== nzc-2)
	  im[i][j][k]=0.0;
    }
      
  
  
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);
  // deallocate temporary array and objects
  delArr3(temp, nxc, nyc);
  delArr3(im, nxc, nyc);
}
/*! interpolate charge density and pressure density from node to center */
void EMfields3D::interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid) {
  // do we need communication or not really?
  grid->interpN2C_GC(rhoc, rhon);
}
/*! communicate ghost for grid -> Particles interpolation */
void EMfields3D::communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct, Grid * grid) {
  // NB: these bc passed are the particle ones
  // interpolate adding common nodes among processors
  communicateInterp(nxn, nyn, nzn, ns, rhons, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jxs, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jys, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jzs, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXXsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXYsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXZsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pYYsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pYZsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pZZsn, 0, 0, 0, 0, 0, 0, vct);

   
  // calculate the correct densities on the boundaries
  // mlmd: not to do if bcPface... <0 (the mlmd BC conditions)
  adjustNonPeriodicDensities(ns, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);
  // put the correct values on ghost cells

  communicateNode_P(nxn, nyn, nzn, rhons, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jxs, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jys, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jzs, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXXsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXYsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXZsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pYYsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pYZsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pZZsn, ns, vct);

}


/** add an amount of charge density to charge density field at node X,Y */
void Moments::addRho(double weight[][2][2], int X, int Y, int Z) {

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	rho[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
void Moments::addJx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	Jx[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Y to current density field on the node */
void Moments::addJy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	Jy[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Z to current density field on the node */
void Moments::addJz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	Jz[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XX to current density field on the node */
void Moments::addPxx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;
	
	pXX[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XY to current density field on the node*/
void Moments::addPxy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	pXY[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XZ to current density field on the node */
void Moments::addPxz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	pXZ[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction YY to current density field on the node*/
void Moments::addPyy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	pYY[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction YZ to current density field on the node */
void Moments::addPyz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;
	
	pYZ[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
void Moments::addPzz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
	const double temp = weight[i][j][k] * invVOL;

	pZZ[X - i][Y - j][Z - k] += temp;
      }
}

void EMfields3D::addToSpeciesMoments(const Moments & in, int is) {
  assert_eq(in.get_nx(), nxn);
  assert_eq(in.get_ny(), nyn);
  assert_eq(in.get_nz(), nzn);
  for (register int i = 0; i < nxn; i++) {
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
	rhons[is][i][j][k] += in.get_rho(i, j, k);
	Jxs[is][i][j][k] += in.get_Jx(i, j, k);
	Jys[is][i][j][k] += in.get_Jy(i, j, k);
	Jzs[is][i][j][k] += in.get_Jz(i, j, k);
	pXXsn[is][i][j][k] += in.get_pXX(i, j, k);
	pXYsn[is][i][j][k] += in.get_pXY(i, j, k);
	pXZsn[is][i][j][k] += in.get_pXZ(i, j, k);
	pYYsn[is][i][j][k] += in.get_pYY(i, j, k);
	pYZsn[is][i][j][k] += in.get_pYZ(i, j, k);
	pZZsn[is][i][j][k] += in.get_pZZ(i, j, k);
      }
  }
}

/*! add an amount of charge density to charge density field at node X,Y */
void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is) {
   
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	rhons[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}

void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
   
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles                
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	rhons[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}

/*! add an amount of charge density to charge density field at node X,Y */
void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is, VirtualTopology3D * vct, double xp, double yp, double zp) {
   
  int R= vct->getCartesian_rank();
  int C0= vct->getCoordinates(0);
  int C1= vct->getCoordinates(1);
  int C2= vct->getCoordinates(2);
  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  if (X<0 or X>nxn-1){
    cout << "Grid " << numGrid << " R " << R << " particle tries to accumulate outside of grid " << endl;
    cout << "Grid " << numGrid << " R " << R << "Inside add rho: X " << X << " of " << nxn << endl;
    cout << "Grid " << numGrid << " R " << R << " x: " << xp << " -dx " << -dx <<" Lx+dx " << Lx+dx << endl;
    cout << "Grid " << numGrid << " R " << R << " " << C0 << "/ " << XLEN << endl;
    //return;
    int j; for(int i=0; i<10000; i++) j++;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (Y<0 or Y>nyn-1){
    cout << "Grid " << numGrid << " R " << R << " particle tries to accumulate outside of grid " << endl;
    cout << "Grid " << numGrid << " R " << R << "Inside add rho: Y " << Y << " of " << nyn << endl;
    cout << "Grid " << numGrid << " R " << R <<" y: " << yp << " -dy " << -dy <<" Ly+dy " << Ly+dy << endl;    
    cout << "Grid " << numGrid << " R " << R << " " << C1 << "/ " << YLEN << endl;
    //return;
    int j; for(int i=0; i<10000; i++) j++;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if (Z<0 or Z>nzn-1){ 
    cout << "Grid " << numGrid << " R " << R << " particle tries to accumulate outside of grid " << endl;
    cout << "Grid " << numGrid << " R " << R << "Inside add rho: Z " << Z << " of " << nzn << endl;
    cout << "Grid " << numGrid << " R " << R <<" z: " << zp << " -dz " << -dz <<" Lz+dz " << Lz+dz << endl;
    cout << "Grid " << numGrid << " R " << R << " " << C2 << "/ " << ZLEN << endl;
    //return;
    int j; for (int i=0; i<10000; i++) j++;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	rhons[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of charge density to current density - direction X to current density field on the node */
void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	Jxs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Y to current density field on the node */
void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	Jys[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of current density - direction Z to current density field on the node */
void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	Jzs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XX to current density field on the node */
void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pXXsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XY to current density field on the node */
void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pXYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction XZ to current density field on the node */
void EMfields3D::addPxz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pXZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YY to current density field on the node */
void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pYYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction YZ to current density field on the node */
void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pYZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}
/*! add an amount of pressure density - direction ZZ to current density field on the node */
void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	pZZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
}

/*! add an amount of charge density to current density - direction X to current density field on the node */
void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles                  
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;
	
	Jxs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of current density - direction Y to current density field on the node */
void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles            
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	Jys[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of current density - direction Z to current density field on the node */
void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles                
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	Jzs[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction XX to current density field on the node */
void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles         
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pXXsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction XY to current density field on the node */
void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles            
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pXYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction XZ to current density field on the node */
void EMfields3D::addPxz(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pXZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction YY to current density field on the node */
void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pYYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction YZ to current density field on the node */
void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pYZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}
/*! add an amount of pressure density - direction ZZ to current density field on the node */
void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is, int nxn, int nyn, int nzn) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

	// added, to avoid segm fault of repopulated particles
	if (X-i<0 or Y - j<0 or Z - k<0 ) continue;
	if (X-i>nxn-1 or Y - j> nyn-1 or Z - k> nzn-1) continue;

	pZZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
      }
}

/*! set to 0 all the densities fields */
void EMfields3D::setZeroDensities() {
  for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
	Jx  [i][j][k] = 0.0;
	Jxh [i][j][k] = 0.0;
	Jy  [i][j][k] = 0.0;
	Jyh [i][j][k] = 0.0;
	Jz  [i][j][k] = 0.0;
	Jzh [i][j][k] = 0.0;
	rhon[i][j][k] = 0.0;
      }
  for (register int i = 0; i < nxc; i++)
    for (register int j = 0; j < nyc; j++)
      for (register int k = 0; k < nzc; k++) {
	rhoc[i][j][k] = 0.0;
	rhoh[i][j][k] = 0.0;
      }
  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
	for (register int k = 0; k < nzn; k++) {
	  rhons[kk][i][j][k] = 0.0;
	  Jxs  [kk][i][j][k] = 0.0;
	  Jys  [kk][i][j][k] = 0.0;
	  Jzs  [kk][i][j][k] = 0.0;
	  pXXsn[kk][i][j][k] = 0.0;
	  pXYsn[kk][i][j][k] = 0.0;
	  pXZsn[kk][i][j][k] = 0.0;
	  pYYsn[kk][i][j][k] = 0.0;
	  pYZsn[kk][i][j][k] = 0.0;
	  pZZsn[kk][i][j][k] = 0.0;
	}

}
/*!SPECIES: Sum the charge density of different species on NODES */
void EMfields3D::sumOverSpecies(VirtualTopology3D * vct) {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
	for (register int k = 0; k < nzn; k++)
	  rhon[i][j][k] += rhons[is][i][j][k];
}

/*!SPECIES: Sum current density for different species */
void EMfields3D::sumOverSpeciesJ() {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
	for (register int k = 0; k < nzn; k++) {
	  Jx[i][j][k] += Jxs[is][i][j][k];
	  Jy[i][j][k] += Jys[is][i][j][k];
	  Jz[i][j][k] += Jzs[is][i][j][k];
	}
}



/*! initialize Magnetic and Electric Field with initial configuration */
void EMfields3D::init(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++) {
      for (int j = 0; j < nyn; j++) {
	for (int k = 0; k < nzn; k++) {
	  for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  Bxn[i][j][k] = B0x;
	  Byn[i][j][k] = B0y;
	  Bzn[i][j][k] = B0z;
	}
      }
    }

    // initialize B on centers
    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);

    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {                        // READING FROM RESTART
    if (vct->getCartesian_rank() == 0)
      cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;


    stringstream ss;
    stringstream ss1;
    ss << vct->getCartesian_rank();
    ss1 << numGrid ;
    string name_file = RestartDirName + "/restart" + ss.str() + "_G" + ss1.str() + ".hdf";
    
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << name_file << endl;
      cout << "RESTART NOT POSSIBLE" << endl;
    }

    dataset_id = H5Dopen2(file_id, "/fields/Bx/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);



    // Bxn
    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    dataset_id = H5Dopen2(file_id, "/fields/By/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    dataset_id = H5Dopen2(file_id, "/fields/Bz/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    dataset_id = H5Dopen2(file_id, "/fields/Ex/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey 
    dataset_id = H5Dopen2(file_id, "/fields/Ey/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez 
    dataset_id = H5Dopen2(file_id, "/fields/Ez/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int jj = 1; jj < nzn - 1; jj++)
	  Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species

    stringstream *species_name = new stringstream[ns];
    for (int is = 0; is < ns; is++) {
      species_name[is] << is;
      string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
      dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
	for (int j = 1; j < nyn - 1; j++)
	  for (int jj = 1; jj < nzn - 1; jj++)
	    rhons[is][i][j][jj] = temp_storage[k++];

      communicateNode_P(nxn, nyn, nzn, rhons, is, vct);
      status = H5Dclose(dataset_id);

    }

    if (col->getCase()=="Dipole") {
      ConstantChargePlanet(grid, vct, col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
    }

    // ConstantChargeOpenBC(grid, vct);

    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // communicate E
    communicateNodeBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
    delete[]species_name;
  }
  


}

/*! initiliaze EM for GEM challange */
void EMfields3D::initBATSRUS(VirtualTopology3D * vct, Grid * grid, Collective *col) {
#ifdef BATSRUS
  cout << "------------------------------------------" << endl;
  cout << "         Initialize from BATSRUS          " << endl;
  cout << "------------------------------------------" << endl;

  // loop over species and cell centers: fill in charge density
  for (int is=0; is < ns; is++)
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for (int k=0; k < nzc; k++)
	  {
	    // WARNING getFluidRhoCenter contains "case" statment
	    rhocs[is][i][j][k] = col->getFluidRhoCenter(i,j,k,is);
	  }

  // loop over cell centers and fill in magnetic and electric fields
  for (int i=0; i < nxc; i++)
    for (int j=0; j < nyc; j++)
      for (int k=0; k < nzc; k++)
	{
	  // WARNING getFluidRhoCenter contains "case" statment
	  col->setFluidFieldsCenter(&Ex[i][j][k],&Ey[i][j][k],&Ez[i][j][k],
				    &Bxc[i][j][k],&Byc[i][j][k],&Bzc[i][j][k],i,j,k);
	}

  // interpolate from cell centers to nodes (corners of cells)
  for (int is=0 ; is<ns; is++)
    grid->interpC2N(rhons[is],rhocs[is]);
  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);
#endif
}

/*! initiliaze EM for GEM challange */
void EMfields3D::initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;


  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;


  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge (MLMD-ready) with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {

	  globalx= grid->getXN(i, j, k) + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + grid->getOz_SW();
         
	  double xpert;
	  double ypert;

	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((globaly - Ly / 2) / delta) * cosh((globaly - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * tanh((globaly - Ly / 2) / delta);
	  // add the initial GEM perturbation
	  // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
	  Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
	  // add the initial X perturbation
	  xpert = globalx - Lx / 2;
	  ypert = globaly - Ly / 2;
	  exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
	  Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
   
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);

    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


/* non-mlmd version
void EMfields3D::initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
	  // add the initial GEM perturbation
	  // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
	  Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
	  // add the initial X perturbation
	  xpert = grid->getXN(i, j, k) - Lx / 2;
	  ypert = grid->getYN(i, j, k) - Ly / 2;
	  exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
	  Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  // Magnetic field
	  Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
	  // add the initial GEM perturbation
	  // Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly );
	  Byc[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
	  // add the initial X perturbation
	  xpert = grid->getXC(i, j, k) - Lx / 2;
	  ypert = grid->getYC(i, j, k) - Ly / 2;
	  exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
	  Bxc[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  Byc[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  // guide field
	  Bzc[i][j][k] = B0z;

	}
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
    }
    }*/


void EMfields3D::initOriginalGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  const double yM = grid->getYN(i, j, k) - .5 * Ly;
	  Bxn[i][j][k] = B0x * tanh(yM / delta);
	  // add the initial GEM perturbation
	  const double xM = grid->getXN(i, j, k) - .5 * Lx;
	  Bxn[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
	  Byn[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  // Magnetic field
	  const double yM = grid->getYC(i, j, k) - .5 * Ly;
	  Bxc[i][j][k] = B0x * tanh(yM / delta);
	  // add the initial GEM perturbation
	  const double xM = grid->getXC(i, j, k) - .5 * Lx;
	  Bxc[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
	  Byc[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
	  Bzc[i][j][k] = B0z;
	}
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double deltax = 8. * delta;
  const double deltay = 4. * delta;
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  const double xM = grid->getXN(i, j, k) - .5 * Lx;
	  const double yB = grid->getYN(i, j, k) - .25 * Ly;
	  const double yT = grid->getYN(i, j, k) - .75 * Ly;
	  const double yBd = yB / delta;
	  const double yTd = yT / delta;
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is]) {
	      const double sech_yBd = 1. / cosh(yBd);
	      const double sech_yTd = 1. / cosh(yTd);
	      rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
	      rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
	    }
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
	  // add the initial GEM perturbation
	  Bxn[i][j][k] += 0.;
	  Byn[i][j][k] = B0y;
	  // add the initial X perturbation
	  const double xMdx = xM / deltax;
	  const double yBdy = yB / deltay;
	  const double yTdy = yT / deltay;
	  const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
	  Bxn[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
	  Byn[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
	  // add the second initial X perturbation
	  const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
	  Bxn[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
	  Byn[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);

	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  const double xM = grid->getXN(i, j, k) - .5 * Lx;
	  const double yB = grid->getYN(i, j, k) - .25 * Ly;
	  const double yT = grid->getYN(i, j, k) - .75 * Ly;
	  const double yBd = yB / delta;
	  const double yTd = yT / delta;
	  Bxc[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
	  // add the initial GEM perturbation
	  Bxc[i][j][k] += 0.;
	  Byc[i][j][k] = B0y;
	  // add the initial X perturbation
	  const double xMdx = xM / deltax;
	  const double yBdy = yB / deltay;
	  const double yTdy = yT / deltay;
	  const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
	  Bxc[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
	  Byc[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
	  // add the second initial X perturbation
	  const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
	  Bxc[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
	  Byc[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);
	  // guide field
	  Bzc[i][j][k] = B0z;
	}
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


/*! initialize GEM challenge with no Perturbation with dipole-like tail topology */
void EMfields3D::initGEMDipoleLikeTailNoPert(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  // parameters controling the field topology
  // e.g., x1=Lx/5,x2=Lx/4 give 'separated' fields, x1=Lx/4,x2=Lx/3 give 'reconnected' topology

  double x1 = Lx / 6.0;         // minimal position of the gaussian peak 
  double x2 = Lx / 4.0;         // maximal position of the gaussian peak (the one closer to the center)
  double sigma = Lx / 15;       // base sigma of the gaussian - later it changes with the grid
  double stretch_curve = 2.0;   // stretch the sin^2 function over the x dimension - also can regulate the number of 'knots/reconnecitons points' if less than 1
  double skew_parameter = 0.50; // skew of the shape of the gaussian
  double pi = 3.1415927;
  double r1, r2, delta_x1x2;

  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field

	  delta_x1x2 = x1 - x2 * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

	  r1 = (grid->getYN(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
	  r2 = (grid->getYN(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

	  // tail-like field topology
	  Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

	  Byn[i][j][k] = B0y;
	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  // Magnetic field

	  delta_x1x2 = x1 - x2 * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

	  r1 = (grid->getYC(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
	  r2 = (grid->getYC(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

	  // tail-like field topology
	  Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

	  Byc[i][j][k] = B0y;
	  // guide field
	  Bzc[i][j][k] = B0z;

	}
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }

}

/*! initialize GEM challenge with no Perturbation */
void EMfields3D::initGEMnoPert(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
	  Byn[i][j][k] = B0y;
	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  // Magnetic field
	  Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
	  Byc[i][j][k] = B0y;
	  // guide field
	  Bzc[i][j][k] = B0z;

	}
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initRandomField(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  double **modes_seed = newArr2(double, 7, 7);
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    double phixy;
    double phix;
    double phiy;
    double phiz;
    double kx;
    double ky;
    phixy = rand() / (double) RAND_MAX;
    phiz = rand() / (double) RAND_MAX;
    phix = rand() / (double) RAND_MAX;
    phiy = rand() / (double) RAND_MAX;
    for (int m = -3; m < 4; m++)
      for (int n = -3; n < 4; n++) {
	modes_seed[m + 3][n + 3] = rand() / (double) RAND_MAX;
      }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = 0.0;
	  Byn[i][j][k] = 0.0;
	  Bzn[i][j][k] = 0.0;
	  for (int m = -3; m < 4; m++)
	    for (int n = -3; n < 4; n++) {

	      kx = 2.0 * M_PI * m / Lx;
	      ky = 2.0 * M_PI * n / Ly;
	      Bxn[i][j][k] += -B0x * ky * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
	      Byn[i][j][k] += B0x * kx * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
	      Bzn[i][j][k] += B0x * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
	    }

	  /* for (int m=1; m < 4; m++) for (int n=1; n < 4; n++){ kx=2.0*M_PI*m/Lx; ky=2.0*M_PI*n/Ly; Bxn[i][j][k] += B0x/kx*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phixy); Byn[i][j][k] += B0x/ky*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phixy); Bzn[i][j][k] += B0x/(kx+ky)*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phiz); } for(int n=1; n < 4; n++){ ky=2.0*M_PI*n/Ly; Bxn[i][j][k] += B0x/(2.0*M_PI/Lx)*cos(grid->getYN(i,j,k)*ky+2.0*M_PI*phix); } for(int m=1; m < 4; m++){ kx=2.0*M_PI*m/Lx; Byn[i][j][k] += B0x/(2.0*M_PI/Ly)*cos(grid->getXN(i,j,k)*kx+2.0*M_PI*phiy); } */
	}
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
  delArr2(modes_seed, 7);
}

/*! Init Force Free (JxB=0) */
void EMfields3D::initForceFree(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------" << endl;
      cout << "Initialize Force Free with Perturbation" << endl;
      cout << "----------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
      }
      cout << "Smoothing Factor = " << Smooth << endl;
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
	  // add the initial GEM perturbation
	  Bxn[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXN(i, j, k) / Lx) * sin(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
	  Byn[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXN(i, j, k) / Lx) * cos(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
	  // guide field
	  Bzn[i][j][k] = B0z / cosh((grid->getYN(i, j, k) - Ly / 2) / delta);
	}
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
	  // add the perturbation
	  Bxc[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXC(i, j, k) / Lx) * sin(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
	  Byc[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXC(i, j, k) / Lx) * cos(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
	  // guide field
	  Bzc[i][j][k] = B0z / cosh((grid->getYC(i, j, k) - Ly / 2) / delta);
	}

    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}
/*! Initialize the EM field with constants values or from restart */
void EMfields3D::initBEAM(VirtualTopology3D * vct, Grid * grid, Collective *col, double x_center, double y_center, double z_center, double radius) {
  double distance;
  // initialize E and rhos on nodes
  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  Bxn[i][j][k] = 0.0;
	  Byn[i][j][k] = 0.0;
	  Bzn[i][j][k] = 0.0;
	  distance = (grid->getXN(i, j, k) - x_center) * (grid->getXN(i, j, k) - x_center) / (radius * radius) + (grid->getYN(i, j, k) - y_center) * (grid->getYN(i, j, k) - y_center) / (radius * radius) + (grid->getZN(i, j, k) - z_center) * (grid->getZN(i, j, k) - z_center) / (4 * radius * radius);
	  // plasma
	  rhons[0][i][j][k] = rhoINIT[0] / FourPI;  // initialize with constant density
	  // electrons
	  rhons[1][i][j][k] = rhoINIT[1] / FourPI;
	  // beam
	  if (distance < 1.0)
	    rhons[2][i][j][k] = rhoINIT[2] / FourPI;
	  else
	    rhons[2][i][j][k] = 0.0;
	}
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
	for (int k = 0; k < nzc; k++) {
	  // Magnetic field
	  Bxc[i][j][k] = 0.0;
	  Byc[i][j][k] = 0.0;
	  Bzc[i][j][k] = 0.0;


	}
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {                        // EM initialization from RESTART
    init(vct, grid, col);            // use the fields from restart file
  }

}

void EMfields3D::UpdateCycle(int cycle){
  currentCycle= cycle;
}

void EMfields3D::UpdateFext(int cycle){


  if (numGrid >0 ) return; // only CG MAY do this

  /* -- NOTE: Hardcoded option -- */
  enum   {LINEAR,STAIRCASE};
  int    utype = STAIRCASE;
  /* -- END NOTE -- */

  double t_beg = 500.0;
  double t_end = 4500.0;
  double Fmin  = 0.1;
  double Fmax  = 1.0;

  double m     = (Fmax - Fmin) / (t_end - t_beg);
  double b     = Fmax - m*t_end;

  if (utype==LINEAR) {
    Fext = m * cycle + b;
  }
  else {
    // Staircase function in 10 steps:
    if (cycle%int((t_end-t_beg)/10) == 0) Fext += (Fmax-Fmin)/10.0;
  }

  if (cycle < t_beg) Fext = Fmin;
  if (cycle > t_end) Fext = Fmax;

}

double EMfields3D::getFext(){
  return (Fext);
}

void EMfields3D::SetDipole_3Bext(VirtualTopology3D *vct, Grid *grid, Collective *col){

  double BE0 = sqrt(B1x*B1x + B1y*B1y + B1z*B1z);
  double a=delta;

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){

	double xc=x_center;
	double yc=y_center;
	double zc=z_center;

	double x = grid->getXN(i,j,k);
	double y = grid->getYN(i,j,k);
	double z = grid->getZN(i,j,k);

	double rx = x-xc;
	//3d: double ry = y-yc;
	double ry = 0.0;
	double rz = z-zc;

	double r         = sqrt(rx*rx + ry*ry + rz*rz);
	double costheta  = fabs(rz) / r;
	double sintheta  = sqrt( 1 - costheta*costheta);
	//3d: double sinphi    = fabs(rx) / (r*sintheta);
	//3d: double cosphi    = fabs(ry) / (r*sintheta);
	double sinphi    = 1.0;
	double cosphi    = 0.0;

	double aor = a/r;
	double Br = -2 * BE0 * aor*aor*aor * costheta;
	double Bt =    - BE0 * aor*aor*aor * sintheta;

	if (r>a) {
	  Bx_ext[i][j][k] =   Bt * costheta * sinphi + Br * sintheta * sinphi;
	  By_ext[i][j][k] = - Bt * costheta * cosphi + Br * sintheta * cosphi;
	  Bz_ext[i][j][k] =   Bt * sintheta          + Br * costheta;
	}
	else {
	  Bx_ext[i][j][k] = B0x;
	  By_ext[i][j][k] = B0y;
	  Bz_ext[i][j][k] = B0z;
	}

	//if (r<1.0*a) cout << " >> " << r << " " << ry << " ct=" << costheta << " st=" << sintheta << " cp=" << cosphi << " sp=" << sinphi << " ar=" << aor << " Br=" << Br << " Bt=" << Bt << " Bx=" << Bx_ext[i][j][k] << " By=" << By_ext[i][j][k] << " Bz=" << Bz_ext[i][j][k] << endl;

      }
    }
  }

  //UpdateRHOcs(grid);

}

void EMfields3D::SetDipole_2Bext(VirtualTopology3D *vct, Grid *grid, Collective *col){

  /* -- NOTE: Hardcoded option */
  bool twodim = false;
  /* -- END NOTE -- */

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){

	double a=0.75*delta; // 0.75 x To avoid problems with de-centered dipoles

	double xc=x_center;
	double yc=y_center;
	double zc=z_center;

	double x = grid->getXN(i,j,k);
	double y = grid->getYN(i,j,k);
	double z = grid->getZN(i,j,k);

	double rx = x-xc;
	double ry = y-yc;
	double rz = z-zc;
	if (twodim) ry = 0.0;

	double r      = sqrt(rx*rx + ry*ry + rz*rz);

	double Mx = B1x;
	double My = B1y;
	double Mz = B1z;

	if (r < 1.0*a) {
	  rz = sqrt(a*a - (rx*rx + ry*ry));
	  double one_r3 = 1.0/(a*a*a);
	  double rhx = rx/a;
	  double rhy = ry/a;
	  double rhz = rz/a;

	  double Bxe = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhx - Mx );
	  double Bye = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhy - My );
	  double Bze = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhz - Mz );
	  double Bme = sqrt(Bxe*Bxe + Bye*Bye + Bze*Bze);

	  Bx_ext[i][j][k] = 0.0;
	  By_ext[i][j][k] = 0.0;
	  Bz_ext[i][j][k] = (B1z/fabs(B1z)) * Bme;
	}
	else {

	  double one_r3 = 1.0/(r*r*r);

	  double rhx = rx/r;
	  double rhy = ry/r;
	  double rhz = rz/r;

	  Bx_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhx - Mx );
	  By_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhy - My );
	  Bz_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhz - Mz );
	}

      }
    }
  }

  //UpdateRHOcs(grid);

}

void EMfields3D::initDipole_2(VirtualTopology3D *vct, Grid *grid, Collective *col){

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){

	double a=delta;

	double xc=x_center;
	double yc=y_center;
	double zc=z_center;

	double x = grid->getXN(i,j,k);
	double y = grid->getYN(i,j,k);
	double z = grid->getZN(i,j,k);

	double rx = x-xc;
	double ry = y-yc;
	double rz = z-zc;

	double r      = sqrt(rx*rx + ry*ry + rz*rz);

	if (r > a/5.0) {
	  double one_r3 = 1.0/(r*r*r);

	  double rhx = rx/r;
	  double rhy = ry/r;
	  double rhz = rz/r;

	  double Mx = B1x;
	  double My = B1y;
	  double Mz = B1z;

	  Bx_ext[i][j][k] = B0x + one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhx - Mx );
	  By_ext[i][j][k] = B0y + one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhy - My );
	  Bz_ext[i][j][k] = B0z + one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhz - Mz );
	}
	else{
	  Bx_ext[i][j][k] = B0x;
	  By_ext[i][j][k] = B0y;
	  Bz_ext[i][j][k] = B0z;
	}

	Bxn[i][j][k] = B0x;
	Byn[i][j][k] = B0y;
	Bzn[i][j][k] = B0z;
	//Bxn[i][j][k] += Bx_ext[i][j][k];
	//Byn[i][j][k] += By_ext[i][j][k];
	//Bzn[i][j][k] += Bz_ext[i][j][k];

      }
    }
  }

  grid->interpN2C(Bxc,Bxn);
  grid->interpN2C(Byc,Byn);
  grid->interpN2C(Bzc,Bzn);

  communicateCenterBC_P(nxc,nyc,nzc,Bxc,col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5],vct);
  communicateCenterBC_P(nxc,nyc,nzc,Byc,col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5],vct);
  communicateCenterBC_P(nxc,nyc,nzc,Bzc,col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5],vct);

  UpdateRHOcs(grid);

}

void EMfields3D::UpdateRHOcs(Grid * grid){

  double r = 1.0;

  double xmin = 0.0;
  double xmax = Lx;
  double rmin = 1.0;
  double rmax = 1.0;

  for (int is = 0; is < ns; is++)
    for (int i=0; i<nxc-1; i++)
      for (int j=0; j<nyc-1; j++)
	for (int k=0; k<nzc-1; k++){
	  double x = grid->getXN(i,j,k);
	  r = rmin + (rmax-rmin) * (xmax - x) / (xmax - xmin);
	  if (r<xmin) r = rmax;
	  if (r>xmax) r = rmin;
	  rhocs[is][i][j][k] = (qom[is]/fabs(qom[is]))*(r/FourPI);
	}
  //grid->interpN2C(rhocs, is, rhons);
}

void EMfields3D::SetLambda(Grid *grid, VirtualTopology3D * vct){

  bool SetDamping= true;
  
  if (numGrid >0 ) SetDamping= false;

  if (numGrid==0) {BufX=4; BufY=4; BufZ=4;}

  if (vct->getCartesian_rank() ==0 and SetDamping)
    cout << "Grid " << numGrid << " is initialising a Lambda layer " << endl;

  double Fac=1.0;

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){

	Lambda[i][j][k]= 0.0;

	// the buffering area left starts at 2, BufLen included
	// the buffering area right starts at nxn-3, BufLen included
	if (SetDamping){

	  // X
	  if (vct->getXleft_neighbor()== MPI_PROC_NULL and i>1 and i< 2+ BufX){
	    //Fac= BufferFactor('L', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dx*Fac;
	  }

	  if (vct->getXright_neighbor()== MPI_PROC_NULL and i< nxn -2 and i> nxn-3-BufX){
	    //Fac= BufferFactor('R', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dx * Fac;
	  }

	  // Y
	  if (vct->getYleft_neighbor()== MPI_PROC_NULL and j>1 and j< 2+ BufY){
	    //Fac= BufferFactor('F', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dy * Fac;
	  }

	  if (vct->getYright_neighbor()== MPI_PROC_NULL and j< nyn -2 and j> nyn-3-BufY){
	    //Fac= BufferFactor('b', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dy * Fac;
	  }

	  // Z
	  if (vct->getZleft_neighbor()== MPI_PROC_NULL and k>1 and k< 2+ BufZ){
	    //Fac= BufferFactor('B', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dz * Fac;
	  }

	  if (vct->getZright_neighbor()== MPI_PROC_NULL and k< nzn -2 and k> nzn-3-BufZ){
	    //Fac= BufferFactor('T', i, j, k, BufX, BufY, BufZ);
	    Lambda[i][j][k]= 2.0*M_PI/dz * Fac;
	    }

	} // end SetDamping

	/*double x = grid->getXN(i,j,k);
	double y = grid->getYN(i,j,k);
	double z = grid->getZN(i,j,k);

	double xmin_r = Lx - 75.0 * dx;
	double xmax_r = Lx - 25.0  * dx;

	Lambda[i][j][k] = 0.0;

	if (x > xmin_r) {
	  if (x < xmax_r) Lambda[i][j][k] = ((x - xmin_r) /  (xmax_r - xmin_r)) * 4.0 * M_PI / dx;
	  else            Lambda[i][j][k] = 4.0 * M_PI / dx;
	  }*/

      }
    }
  }

}

double*** EMfields3D::GetLambda(){
  return(Lambda);
}

/*! Initialise a combination of magnetic dipoles */
void EMfields3D::initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col){

  double distance;

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){
	for (int is=0; is < ns; is++){
	  rhons[is][i][j][k] = rhoINIT[is]/FourPI;
	}
	Ex[i][j][k] = ebc[0];
	Ey[i][j][k] = ebc[1];
	Ez[i][j][k] = ebc[2];

	double blp[3];
	// Set coil diameter
	double a=delta;

	double xc=x_center;
	double yc=y_center;
	double zc=z_center;

	double x = grid->getXN(i,j,k);
	double y = grid->getYN(i,j,k);
	double z = grid->getZN(i,j,k);

	double r2 = ((x-xc)*(x-xc)) + ((y-yc)*(y-yc)) + ((z-zc)*(z-zc));

	// Compute dipolar field B_ext

	if (r2 > delta*delta) {
	  loopZ(blp, x, y, z, a, xc, yc, zc, B1z);
	  Bx_ext[i][j][k]  = blp[0];
	  By_ext[i][j][k]  = blp[1];
	  Bz_ext[i][j][k]  = blp[2];
	  loopX(blp, x, y, z, a, xc, yc, zc, B1x);
	  Bx_ext[i][j][k] += blp[0];
	  By_ext[i][j][k] += blp[1];
	  Bz_ext[i][j][k] += blp[2];
	  loopY(blp, x, y, z, a, xc, yc, zc, B1y);
	  Bx_ext[i][j][k] += blp[0];
	  By_ext[i][j][k] += blp[1];
	  Bz_ext[i][j][k] += blp[2];
	}
	else {
	  Bx_ext[i][j][k]  = 0.0;
	  By_ext[i][j][k]  = 0.0;
	  Bz_ext[i][j][k]  = 0.0;
	}

	Bxn[i][j][k] = B0x + Bx_ext[i][j][k];
	Byn[i][j][k] = B0y + By_ext[i][j][k];
	Bzn[i][j][k] = B0z + Bz_ext[i][j][k];

	// -- Uncomment if using the J_ext method:
	// Bx_ext[i][j][k]  = 0.0;
	// By_ext[i][j][k]  = 0.0;
	// Bz_ext[i][j][k]  = 0.0;
	// -- end Uncomment
      }
    }
  }

  grid->interpN2C(Bxc,Bxn);
  grid->interpN2C(Byc,Byn);
  grid->interpN2C(Bzc,Bzn);

  communicateCenterBC_P(nxc,nyc,nzc,Bxc,col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5],vct);
  communicateCenterBC_P(nxc,nyc,nzc,Byc,col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5],vct);
  communicateCenterBC_P(nxc,nyc,nzc,Bzc,col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5],vct);

  // -- initialize J_ext =c/4*pi curl(B) on nodes (current due to the dipole)
  // -- WHY IS THIS WORKING ANYWAYS?? grid->curlC2N(tempXN,tempYN,tempZN,Bxc_ext,Byc_ext,Bzc_ext);
  // grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc);
  // scale(Jx_ext,tempXN,c/FourPI,nxn,nyn,nzn);
  // scale(Jy_ext,tempYN,c/FourPI,nxn,nyn,nzn);
  // scale(Jz_ext,tempZN,c/FourPI,nxn,nyn,nzn);
  // for (int i=0; i < nxn; i++){
  //   for (int j=0; j < nyn; j++){
  //     for (int k=0; k < nzn; k++){
  //       if (i<3 || i>nxn-4) {
  //         Jx_ext[i][j][k] = 0.0;
  //         Jy_ext[i][j][k] = 0.0;
  //         Jz_ext[i][j][k] = 0.0;
  //       }
  //       if (j<3 || j>nyn-4) {
  //         Jx_ext[i][j][k] = 0.0;
  //         Jy_ext[i][j][k] = 0.0;
  //         Jz_ext[i][j][k] = 0.0;
  //       }
  //       if (k<3 || k>nzn-4) {
  //         Jx_ext[i][j][k] = 0.0;
  //         Jy_ext[i][j][k] = 0.0;
  //         Jz_ext[i][j][k] = 0.0;
  //       }
  //     }
  //   }
  // }
  // -- end J_ext

  for (int is=0 ; is<ns; is++)
    grid->interpN2C_GC(rhocs,is,rhons);

  if (restart1 != 0) { // EM initialization from RESTART
    init(vct,grid,col);  // use the fields from restart file
  }

}


/*! Calculate the susceptibility on the X boundary */
void EMfields3D::sustensorX(double **susxx, double **susxy, double **susxz, int N) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susxy[j][k] = 0.0;
      susxz[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
	omcx = beta * (Bxn[N][j][k] + Fext*Bx_ext[N][j][k]);
	omcy = beta * (Byn[N][j][k] + Fext*By_ext[N][j][k]);
	omcz = beta * (Bzn[N][j][k] + Fext*Bz_ext[N][j][k]);
	denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][N][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
	susxx[j][k] += (  1.0 + omcx * omcx) * denom;
	susxy[j][k] += ( omcz + omcx * omcy) * denom;
	susxz[j][k] += (-omcy + omcx * omcz) * denom;
      }
  }

}

/*! Calculate the susceptibility on the Y boundary */
void EMfields3D::sustensorY(double **susyx, double **susyy, double **susyz, int N) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susyx[i][k] = 0.0;
      susyy[i][k] = 1.0;
      susyz[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
	omcx = beta * (Bxn[i][N][k] + Fext*Bx_ext[i][N][k]);
	omcy = beta * (Byn[i][N][k] + Fext*By_ext[i][N][k]);
	omcz = beta * (Bzn[i][N][k] + Fext*Bz_ext[i][N][k]);
	denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][N][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
	susyx[i][k] += (-omcz + omcx * omcy) * denom;
	susyy[i][k] += (  1.0 + omcy * omcy) * denom;
	susyz[i][k] += (+omcx + omcy * omcz) * denom;
      }
  }

}

/*! Calculate the susceptibility on the Z boundary */
void EMfields3D::sustensorZ(double **suszx, double **suszy, double **suszz, int N) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      suszx[i][j] = 0.0;
      suszy[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
	omcx = beta * (Bxn[i][j][N] + Fext*Bx_ext[i][j][N]);
	omcy = beta * (Byn[i][j][N] + Fext*By_ext[i][j][N]);
	omcz = beta * (Bzn[i][j][N] + Fext*Bz_ext[i][j][N]);
	denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][N] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
	suszx[i][j] += ( omcy + omcx * omcz) * denom;
	suszy[i][j] += (-omcx + omcy * omcz) * denom;
	suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }

}

/*! Perfect conductor boundary conditions: LEFT wall */
void EMfields3D::perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  // June 10: put known term to the right
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
  case 0:  // boundary condition on X-DIRECTION 
    susxx = newArr2(double,nyn,nzn);
    susxy = newArr2(double,nyn,nzn);
    susxz = newArr2(double,nyn,nzn);
    sustensorX(susxx, susxy, susxz, 1);
    for (int i=1; i <  nyn-1;i++)
      for (int j=1; j <  nzn-1;j++){
	//imageX[1][i][j] = vectorX[1][i][j] - (Ex[1][i][j] - susxy[i][j]*vectorY[1][i][j] - susxz[i][j]*vectorZ[1][i][j] - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j];
	imageX[1][i][j] = vectorX[1][i][j] - ( - susxy[i][j]*vectorY[1][i][j] - susxz[i][j]*vectorZ[1][i][j] )/susxx[i][j];
	imageY[1][i][j] = vectorY[1][i][j] - 0.0*vectorY[2][i][j];
	imageZ[1][i][j] = vectorZ[1][i][j] - 0.0*vectorZ[2][i][j];
      }
    delArr2(susxx,nxn);
    delArr2(susxy,nxn);
    delArr2(susxz,nxn);
    break;
  case 1: // boundary condition on Y-DIRECTION
    susyx = newArr2(double,nxn,nzn);
    susyy = newArr2(double,nxn,nzn);
    susyz = newArr2(double,nxn,nzn);
    sustensorY(susyx, susyy, susyz, 1);
    for (int i=1; i < nxn-1;i++)
      for (int j=1; j <  nzn-1;j++){
	imageX[i][1][j] = vectorX[i][1][j] - 0.0*vectorX[i][2][j];
	//imageY[i][1][j] = vectorY[i][1][j] - (Ey[i][1][j] - susyx[i][j]*vectorX[i][1][j] - susyz[i][j]*vectorZ[i][1][j] - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j];
	imageY[i][1][j] = vectorY[i][1][j] - (- susyx[i][j]*vectorX[i][1][j] - susyz[i][j]*vectorZ[i][1][j] )/susyy[i][j];
	imageZ[i][1][j] = vectorZ[i][1][j] - 0.0*vectorZ[i][2][j];
      }
    delArr2(susyx,nxn);
    delArr2(susyy,nxn);
    delArr2(susyz,nxn);
    break;
  case 2: // boundary condition on Z-DIRECTION
    suszx = newArr2(double,nxn,nyn);
    suszy = newArr2(double,nxn,nyn);
    suszz = newArr2(double,nxn,nyn);
    sustensorZ(suszx, suszy, suszz, 1);
    for (int i=1; i <  nxn-1;i++)
      for (int j=1; j <  nyn-1;j++){
	imageX[i][j][1] = vectorX[i][j][1];
	imageY[i][j][1] = vectorX[i][j][1];
	//imageZ[i][j][1] = vectorZ[i][j][1] - (Ez[i][j][1] - suszx[i][j]*vectorX[i][j][1] - suszy[i][j]*vectorY[i][j][1] - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j];
	imageZ[i][j][1] = vectorZ[i][j][1] - (- suszx[i][j]*vectorX[i][j][1] - suszy[i][j]*vectorY[i][j][1] )/suszz[i][j];
      }
    delArr2(suszx,nxn);
    delArr2(suszy,nxn);
    delArr2(suszz,nxn);
    break;
  }
}


/*! Perfect conductor boundary conditions: RIGHT wall */
void EMfields3D::perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  double beta, omcx, omcy, omcz, denom;
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
  case 0: // boundary condition on X-DIRECTION RIGHT
    susxx = newArr2(double,nyn,nzn);
    susxy = newArr2(double,nyn,nzn);
    susxz = newArr2(double,nyn,nzn);
    sustensorX(susxx, susxy, susxz, nxn-2);
    for (int i=1; i < nyn-1;i++)
      for (int j=1; j <  nzn-1;j++){
	//imageX[nxn-2][i][j] = vectorX[nxn-2][i][j] - (Ex[nxn-2][i][j] - susxy[i][j]*vectorY[nxn-2][i][j] - susxz[i][j]*vectorZ[nxn-2][i][j] - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j];
	imageX[nxn-2][i][j] = vectorX[nxn-2][i][j] - ( - susxy[i][j]*vectorY[nxn-2][i][j] - susxz[i][j]*vectorZ[nxn-2][i][j] )/susxx[i][j];
	imageY[nxn-2][i][j] = vectorY[nxn-2][i][j] - 0.0 * vectorY[nxn-3][i][j];
	imageZ[nxn-2][i][j] = vectorZ[nxn-2][i][j] - 0.0 * vectorZ[nxn-3][i][j];
      }
    delArr2(susxx,nxn);
    delArr2(susxy,nxn);       
    delArr2(susxz,nxn);
    break;
  case 1: // boundary condition on Y-DIRECTION RIGHT
    susyx = newArr2(double,nxn,nzn);
    susyy = newArr2(double,nxn,nzn);
    susyz = newArr2(double,nxn,nzn);
    sustensorY(susyx, susyy, susyz, nyn-2);
    for (int i=1; i < nxn-1;i++)
      for (int j=1; j < nzn-1;j++){
	imageX[i][nyn-2][j] = vectorX[i][nyn-2][j] - 0.0*vectorX[i][nyn-3][j];
	//imageY[i][nyn-2][j] = vectorY[i][nyn-2][j] - (Ey[i][nyn-2][j] - susyx[i][j]*vectorX[i][nyn-2][j] - susyz[i][j]*vectorZ[i][nyn-2][j] - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
	imageY[i][nyn-2][j] = vectorY[i][nyn-2][j] - ( - susyx[i][j]*vectorX[i][nyn-2][j] - susyz[i][j]*vectorZ[i][nyn-2][j] )/susyy[i][j];
	imageZ[i][nyn-2][j] = vectorZ[i][nyn-2][j] - 0.0*vectorZ[i][nyn-3][j];
      }
    delArr2(susyx,nxn);
    delArr2(susyy,nxn);
    delArr2(susyz,nxn);
    break;
  case 2: // boundary condition on Z-DIRECTION RIGHT
    suszx = newArr2(double,nxn,nyn);
    suszy = newArr2(double,nxn,nyn);
    suszz = newArr2(double,nxn,nyn);
    sustensorZ(suszx, suszy, suszz, nzn-2);
    for (int i=1; i < nxn-1;i++)
      for (int j=1; j < nyn-1;j++){
	imageX[i][j][nzn-2] = vectorX[i][j][nzn-2];
	imageY[i][j][nzn-2] = vectorY[i][j][nzn-2];
	//imageZ[i][j][nzn-2] = vectorZ[i][j][nzn-2] - (Ez[i][j][nzn-2] - suszx[i][j]*vectorX[i][j][nzn-2] - suszy[i][j]*vectorY[i][j][nzn-2] - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j];
	imageZ[i][j][nzn-2] = vectorZ[i][j][nzn-2] - (- suszx[i][j]*vectorX[i][j][nzn-2] - suszy[i][j]*vectorY[i][j][nzn-2] )/suszz[i][j];
      }
    delArr2(suszx,nxn);
    delArr2(suszy,nxn);       
    delArr2(suszz,nxn);
    break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;


  switch(dir){
  case 0: // boundary condition on X-DIRECTION LEFT
    susxx = newArr2(double,nyn,nzn);
    susxy = newArr2(double,nyn,nzn);
    susxz = newArr2(double,nyn,nzn);
    sustensorX(susxx, susxy, susxz, 1);

    for (int i=1; i < nyn-1;i++)
      for (int j=1; j < nzn-1;j++){
	vectorX[1][i][j] = -(  - (Ex[1][i][j]  - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j] );
	vectorY[1][i][j] = ebc[1];
	vectorZ[1][i][j] = ebc[2];

	/*vectorX[1][i][j] = 0.0;
	vectorY[1][i][j] = ebc[1];
	vectorZ[1][i][j] = ebc[2];*/
	//+//          vectorX[1][i][j] = 0.0;
	//+//          vectorY[1][i][j] = 0.0;
	//+//          vectorZ[1][i][j] = 0.0;
      }
    delArr2(susxx,nxn);
    delArr2(susxy,nxn);
    delArr2(susxz,nxn);

    break;
  case 1: // boundary condition on Y-DIRECTION LEFT
    susyx = newArr2(double,nxn,nzn);
    susyy = newArr2(double,nxn,nzn);
    susyz = newArr2(double,nxn,nzn);
    sustensorY(susyx, susyy, susyz, 1);

    for (int i=1; i < nxn-1;i++)
      for (int j=1; j < nzn-1;j++){
	vectorX[i][1][j] = ebc[0];
	vectorY[i][1][j] = -( - (Ey[i][1][j]  - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j] );
	vectorZ[i][1][j] = ebc[2];
	
	/*vectorX[i][1][j] = ebc[0];
	vectorY[i][1][j] = 0.0;
	vectorZ[i][1][j] = ebc[2];*/
	//+//          vectorX[i][1][j] = 0.0;
	//+//          vectorY[i][1][j] = 0.0;
	//+//          vectorZ[i][1][j] = 0.0;
      }

    delArr2(susyx,nxn);
    delArr2(susyy,nxn);
    delArr2(susyz,nxn);
    break;
  case 2: // boundary condition on Z-DIRECTION LEFT
    suszx = newArr2(double,nxn,nyn);
    suszy = newArr2(double,nxn,nyn);
    suszz = newArr2(double,nxn,nyn);
    sustensorZ(suszx, suszy, suszz, 1);
    
    for (int i=1; i < nxn-1;i++)
      for (int j=1; j <  nyn-1;j++){
	vectorX[i][j][1] = ebc[0];
	vectorY[i][j][1] = ebc[1];
	vectorZ[i][j][1] = -( - (Ez[i][j][1] - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j]);
	/* vectorX[i][j][1] = ebc[0];
	vectorY[i][j][1] = ebc[1];
	vectorZ[i][j][1] = 0.0; */
	
	//+//          vectorX[i][j][1] = 0.0;
	//+//          vectorY[i][j][1] = 0.0;
	//+//          vectorZ[i][j][1] = 0.0;
      }
    delArr2(suszx,nxn);
    delArr2(suszy,nxn);
    delArr2(suszz,nxn);
    break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::MLMDSourceLeft(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {


  switch(dir){
  case 0: // boundary condition on X-DIRECTION LEFT
    for (int i=0; i < nyn;i++)
      for (int j=0; j < nzn;j++){

	vectorX[1][i][j] = 0.0;
	vectorY[1][i][j] = 0.0;
	vectorZ[1][i][j] = 0.0;
      }
    break;
  case 1: // boundary condition on Y-DIRECTION LEFT
    for (int i=0; i < nxn;i++)
      for (int j=0; j < nzn;j++){
	vectorX[i][1][j] = 0.0;
	vectorY[i][1][j] = 0.0;
	vectorZ[i][1][j] = 0.0;
      }
    break;
  case 2: // boundary condition on Z-DIRECTION LEFT
    for (int i=0; i < nxn;i++)
      for (int j=0; j <  nyn;j++){
	vectorX[i][j][1] = 0.0;
	vectorY[i][j][1] = 0.0;
	vectorZ[i][j][1] = 0.0;
      }
    break;
  }
}


/*! Perfect conductor boundary conditions for source: RIGHT WALL */
void EMfields3D::perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
  case 0: // boundary condition on X-DIRECTION RIGHT

    susxx = newArr2(double,nyn,nzn);
    susxy = newArr2(double,nyn,nzn);
    susxz = newArr2(double,nyn,nzn);
    sustensorX(susxx, susxy, susxz, nxn-2);
    for (int i=1; i < nyn-1;i++)
      for (int j=1; j < nzn-1;j++){
	vectorX[nxn-2][i][j] = - ( - (Ex[nxn-2][i][j] - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j]);
	vectorY[nxn-2][i][j] = ebc[1];
	vectorZ[nxn-2][i][j] = ebc[2];
	/*vectorX[nxn-2][i][j] = 0.0;
	vectorY[nxn-2][i][j] = ebc[1];
	vectorZ[nxn-2][i][j] = ebc[2];*/
	//+//          vectorX[nxn-2][i][j] = 0.0;
	//+//          vectorY[nxn-2][i][j] = 0.0;
	//+//          vectorZ[nxn-2][i][j] = 0.0;
      }
    delArr2(susxx,nxn);
    delArr2(susxy,nxn);
    delArr2(susxz,nxn);
    break;
  case 1: // boundary condition on Y-DIRECTION RIGHT
    susyx = newArr2(double,nxn,nzn);
    susyy = newArr2(double,nxn,nzn);
    susyz = newArr2(double,nxn,nzn);
    sustensorY(susyx, susyy, susyz, nyn-2);

    for (int i=1; i < nxn-1;i++)
      for (int j=1; j < nzn-1;j++){
	vectorX[i][nyn-2][j] = ebc[0];
	vectorY[i][nyn-2][j] = - ( - (Ey[i][nyn-2][j] - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j]);
	vectorZ[i][nyn-2][j] = ebc[2];
	/*vectorX[i][nyn-2][j] = ebc[0];
	vectorY[i][nyn-2][j] = 0.0;
	vectorZ[i][nyn-2][j] = ebc[2];*/
	//+//          vectorX[i][nyn-2][j] = 0.0;
	//+//          vectorY[i][nyn-2][j] = 0.0;
	//+//          vectorZ[i][nyn-2][j] = 0.0;
      }
    delArr2(susyx,nxn);
    delArr2(susyy,nxn);
    delArr2(susyz,nxn);
    break;
  case 2:
    suszx = newArr2(double,nxn,nyn);
    suszy = newArr2(double,nxn,nyn);
    suszz = newArr2(double,nxn,nyn);
    sustensorZ(suszx, suszy, suszz, nzn-2);

    for (int i=1; i <  nxn-1;i++)
      for (int j=1; j <  nyn-1;j++){
	vectorX[i][j][nzn-2] = ebc[0];
	vectorY[i][j][nzn-2] = ebc[1];
	vectorZ[i][j][nzn-2] = - ( - (Ez[i][j][nzn-2] - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j]);
	/*vectorX[i][j][nzn-2] = ebc[0];
	vectorY[i][j][nzn-2] = ebc[1];
	vectorZ[i][j][nzn-2] = 0.0;*/
	//+//          vectorX[i][j][nzn-2] = 0.0;
	//+//          vectorY[i][j][nzn-2] = 0.0;
	//+//          vectorZ[i][j][nzn-2] = 0.0;
      }
    delArr2(suszx,nxn);
    delArr2(suszy,nxn);
    delArr2(suszz,nxn);
    break;
  }
}

void EMfields3D::MLMDSourceRight(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  switch(dir){
  case 0: // boundary condition on X-DIRECTION RIGHT
    for (int i=0; i < nyn;i++)
      for (int j=0; j < nzn;j++){
	vectorX[nxn-2][i][j] = 0.0;
	vectorY[nxn-2][i][j] = 0.0;
	vectorZ[nxn-2][i][j] = 0.0;
      }
    break;
  case 1: // boundary condition on Y-DIRECTION RIGHT
    for (int i=0; i < nxn;i++)
      for (int j=0; j < nzn;j++){
	vectorX[i][nyn-2][j] = 0.0;
	vectorY[i][nyn-2][j] = 0.0;
	vectorZ[i][nyn-2][j] = 0.0;
      }
    break;
  case 2:
    for (int i=0; i <  nxn;i++)
      for (int j=0; j <  nyn;j++){
	vectorX[i][j][nzn-2] = 0.0;
	vectorY[i][j][nzn-2] = 0.0;
	vectorZ[i][j][nzn-2] = 0.0;
      }
    break;
  }
}


// OpenBCs

injInfoFields* EMfields3D::get_InfoFieldsTop() {return injFieldsTop;}
injInfoFields* EMfields3D::get_InfoFieldsBottom() {return injFieldsBottom;}
injInfoFields* EMfields3D::get_InfoFieldsLeft() {return injFieldsLeft;}
injInfoFields* EMfields3D::get_InfoFieldsRight() {return injFieldsRight;}
injInfoFields* EMfields3D::get_InfoFieldsFront() {return injFieldsFront;}
injInfoFields* EMfields3D::get_InfoFieldsRear() {return injFieldsRear;}

// Open Boundary conditions implementation

void EMfields3D::updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col){


  if (numGrid >0) return;
  /* -- NOTE: Hardcoded option -- */
  bool XRightOutflow = false;
  /* -- END NOTE --*/

  double u_0, v_0, w_0;
  u_0=col->getU0(0);
  v_0=col->getV0(0);
  w_0=col->getW0(0);

  if (vct->getXleft_neighbor() == MPI_PROC_NULL)
    {
      for (int i=0; i< 3;i++)
	for (int j=0; j<nyn;j++)
	  for (int k=0; k<nzn;k++){

	    double Bxb;
	    double Byb;
	    double Bzb;
	    double Exb;
	    double Eyb;
	    double Ezb;

	    if (!XRightOutflow) {
	      Bxb = 0.0;
	      Byb = 0.0;
	      Bzb = 0.0;
	      Exb = 0.0;
	      Eyb = 0.0;
	      Ezb = 0.0;
	    }
	    else {
	      Bxb = Bxn[i][j][k];
	      Byb = Byn[i][j][k];
	      Bzb = Bzn[i][j][k];
	      Exb = w_0*Byb-v_0*Bzb;
	      Eyb = u_0*Bzb-w_0*Bxb;
	      Ezb = v_0*Bxb-u_0*Byb;
	    }

	    injFieldsLeft->ExITemp[i][j][k]=Exb;
	    injFieldsLeft->EyITemp[i][j][k]=Eyb;
	    injFieldsLeft->EzITemp[i][j][k]=Ezb;

	    injFieldsLeft->BxITemp[i][j][k]=Bxb;
	    injFieldsLeft->ByITemp[i][j][k]=Byb;
	    injFieldsLeft->BzITemp[i][j][k]=Bzb;
	  }
    }

  if (vct->getXright_neighbor() == MPI_PROC_NULL)
    {
      for (int i=nxn-3; i< nxn; i++)
	for (int j=0; j<nyn; j++)
	  for (int k=0; k<nzn; k++){

	    double Bxb;
	    double Byb;
	    double Bzb;
	    double Exb;
	    double Eyb;
	    double Ezb;

	    if (!XRightOutflow) {
	      Bxb = B0x;
	      Byb = B0y;
	      Bzb = B0z;
	      Exb = w_0*Byb-v_0*Bzb;
	      Eyb = u_0*Bzb-w_0*Bxb;
	      Ezb = v_0*Bxb-u_0*Byb;
	    }
	    else {
	      Bxb = 0.0;
	      Byb = 0.0;
	      Bzb = 0.0;
	      Exb = 0.0;
	      Eyb = 0.0;
	      Ezb = 0.0;
	    }

	    injFieldsRight->ExITemp[i][j][k]=Exb;
	    injFieldsRight->EyITemp[i][j][k]=Eyb;
	    injFieldsRight->EzITemp[i][j][k]=Ezb;

	    injFieldsRight->BxITemp[i][j][k]=Bxb;
	    injFieldsRight->ByITemp[i][j][k]=Byb;
	    injFieldsRight->BzITemp[i][j][k]=Bzb;

	  }

    }

  if (vct->getYleft_neighbor() == MPI_PROC_NULL)
    {
      for (int i=0; i< nxn;i++)
	for (int j=0; j<3;j++)
	  for (int k=0; k<nzn;k++){

	    injFieldsBottom->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
	    injFieldsBottom->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
	    injFieldsBottom->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

	    injFieldsBottom->BxITemp[i][j][k]=B0x;
	    injFieldsBottom->ByITemp[i][j][k]=B0y;
	    injFieldsBottom->BzITemp[i][j][k]=B0z;
	  }

    }
  if (vct->getYright_neighbor() == MPI_PROC_NULL)
    {
      for (int i=0; i< nxn;i++)
	for (int j=nyn-3; j<nyn;j++)
	  for (int k=0; k<nzn;k++){

	    injFieldsTop->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
	    injFieldsTop->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
	    injFieldsTop->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

	    injFieldsTop->BxITemp[i][j][k]=B0x;
	    injFieldsTop->ByITemp[i][j][k]=B0y;
	    injFieldsTop->BzITemp[i][j][k]=B0z;
	  }

    }
  if (vct->getZleft_neighbor() == MPI_PROC_NULL)
    {
      for (int i=0; i< nxn;i++)
	for (int j=0; j<nyn;j++)
	  for (int k=0; k<3;k++){

	    injFieldsRear->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
	    injFieldsRear->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
	    injFieldsRear->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

	    injFieldsRear->BxITemp[i][j][k]=B0x;
	    injFieldsRear->ByITemp[i][j][k]=B0y;
	    injFieldsRear->BzITemp[i][j][k]=B0z;
	  }

    }

  if (vct->getZright_neighbor() == MPI_PROC_NULL)
    {
      for (int i=0; i< nxn;i++)
	for (int j=0; j<nyn;j++)
	  for (int k=nzn-3; k<nzn;k++){

	    injFieldsFront->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
	    injFieldsFront->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
	    injFieldsFront->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

	    injFieldsFront->BxITemp[i][j][k]=B0x;
	    injFieldsFront->ByITemp[i][j][k]=B0y;
	    injFieldsFront->BzITemp[i][j][k]=B0z;
	  }
    }

}

void EMfields3D::BoundaryConditionsEImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz, VirtualTopology3D *vct,Grid *grid){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
	imageX[0][j][k] = vectorX[0][j][k] - injFieldsLeft->ExITemp[0][j][k];
	imageY[0][j][k] = vectorY[0][j][k] - injFieldsLeft->EyITemp[0][j][k];
	imageZ[0][j][k] = vectorZ[0][j][k] - injFieldsLeft->EzITemp[0][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
	imageX[nx-1][j][k] = vectorX[nx-1][j][k]- injFieldsRight->ExITemp[nx-1][j][k];
	imageY[nx-1][j][k] = vectorY[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];
	imageZ[nx-1][j][k] = vectorZ[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];

      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
	imageX[i][0][k] = vectorX[i][0][k]-injFieldsBottom->ExITemp[i][0][k];
	imageY[i][0][k] = vectorY[i][0][k]-injFieldsBottom->EyITemp[i][0][k];
	imageZ[i][0][k] = vectorZ[i][0][k]-injFieldsBottom->EzITemp[i][0][k];
      }

  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
	imageX[i][ny-1][k] = vectorX[i][ny-1][k]-injFieldsTop->ExITemp[i][ny-1][k];
	imageY[i][ny-1][k] = vectorY[i][ny-1][k]-injFieldsTop->EyITemp[i][ny-1][k];
	imageZ[i][ny-1][k] = vectorZ[i][ny-1][k]-injFieldsTop->EzITemp[i][ny-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
	imageX[i][j][0] = vectorX[i][j][0]-injFieldsFront->ExITemp[i][j][0];
	imageY[i][j][0] = vectorY[i][j][0]-injFieldsFront->EyITemp[i][j][0];
	imageZ[i][j][0] = vectorZ[i][j][0]-injFieldsFront->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
	imageX[i][j][nz-1] = vectorX[i][j][nz-1]-injFieldsRear->ExITemp[i][j][nz-1];
	imageY[i][j][nz-1] = vectorY[i][j][nz-1]-injFieldsRear->EyITemp[i][j][nz-1];
	imageZ[i][j][nz-1] = vectorZ[i][j][nz-1]-injFieldsRear->EzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsB(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
	vectorX[0][j][k] = injFieldsLeft->BxITemp[0][j][k];
	vectorY[0][j][k] = injFieldsLeft->ByITemp[0][j][k];
	vectorZ[0][j][k] = injFieldsLeft->BzITemp[0][j][k];

	//      vectorX[1][j][k] = injFieldsLeft->BxITemp[1][j][k];
	//      vectorY[1][j][k] = injFieldsLeft->ByITemp[1][j][k];
	//      vectorZ[1][j][k] = injFieldsLeft->BzITemp[1][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
	//      vectorX[nx-2][j][k] = injFieldsRight->BxITemp[nx-2][j][k];
	//      vectorY[nx-2][j][k] = injFieldsRight->ByITemp[nx-2][j][k];
	//      vectorZ[nx-2][j][k] = injFieldsRight->BzITemp[nx-2][j][k];

	vectorX[nx-1][j][k] = injFieldsRight->BxITemp[nx-1][j][k];
	vectorY[nx-1][j][k] = injFieldsRight->ByITemp[nx-1][j][k];
	vectorZ[nx-1][j][k] = injFieldsRight->BzITemp[nx-1][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//      vectorX[i][1][k] = injFieldsBottom->BxITemp[i][1][k];
	//      vectorY[i][1][k] = injFieldsBottom->ByITemp[i][1][k];
	//      vectorZ[i][1][k] = injFieldsBottom->BzITemp[i][1][k];

	vectorX[i][0][k] = injFieldsBottom->BxITemp[i][0][k];
	vectorY[i][0][k] = injFieldsBottom->ByITemp[i][0][k];
	vectorZ[i][0][k] = injFieldsBottom->BzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//      vectorX[i][ny-2][k] = injFieldsTop->BxITemp[i][ny-2][k];
	//      vectorY[i][ny-2][k] = injFieldsTop->ByITemp[i][ny-2][k];
	//      vectorZ[i][ny-2][k] = injFieldsTop->BzITemp[i][ny-2][k];

	vectorX[i][ny-1][k] = injFieldsTop->BxITemp[i][ny-1][k];
	vectorY[i][ny-1][k] = injFieldsTop->ByITemp[i][ny-1][k];
	vectorZ[i][ny-1][k] = injFieldsTop->BzITemp[i][ny-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//      vectorX[i][j][1] = injFieldsRear->BxITemp[i][j][1];
	//      vectorY[i][j][1] = injFieldsRear->ByITemp[i][j][1];
	//      vectorZ[i][j][1] = injFieldsRear->BzITemp[i][j][1];

	vectorX[i][j][0] = injFieldsRear->BxITemp[i][j][0];
	vectorY[i][j][0] = injFieldsRear->ByITemp[i][j][0];
	vectorZ[i][j][0] = injFieldsRear->BzITemp[i][j][0];
      }
  }


  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//      vectorX[i][j][nz-2] = injFieldsFront->BxITemp[i][j][nz-2];
	//      vectorY[i][j][nz-2] = injFieldsFront->ByITemp[i][j][nz-2];
	//      vectorZ[i][j][nz-2] = injFieldsFront->BzITemp[i][j][nz-2];

	vectorX[i][j][nz-1] = injFieldsFront->BxITemp[i][j][nz-1];
	vectorY[i][j][nz-1] = injFieldsFront->ByITemp[i][j][nz-1];
	vectorZ[i][j][nz-1] = injFieldsFront->BzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsE(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
	//      vectorX[1][j][k] = injFieldsLeft->ExITemp[1][j][k];
	//      vectorY[1][j][k] = injFieldsLeft->EyITemp[1][j][k];
	//      vectorZ[1][j][k] = injFieldsLeft->EzITemp[1][j][k];

	vectorX[0][j][k] = injFieldsLeft->ExITemp[0][j][k];
	vectorY[0][j][k] = injFieldsLeft->EyITemp[0][j][k];
	vectorZ[0][j][k] = injFieldsLeft->EzITemp[0][j][k];
      } 
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){

	//      vectorX[nx-2][j][k] = injFieldsRight->ExITemp[nx-2][j][k];
	//      vectorY[nx-2][j][k] = injFieldsRight->EyITemp[nx-2][j][k];
	//      vectorZ[nx-2][j][k] = injFieldsRight->EzITemp[nx-2][j][k];

	vectorX[nx-1][j][k] = injFieldsRight->ExITemp[nx-1][j][k];
	vectorY[nx-1][j][k] = injFieldsRight->EyITemp[nx-1][j][k];
	vectorZ[nx-1][j][k] = injFieldsRight->EzITemp[nx-1][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//      vectorX[i][1][k] = injFieldsBottom->ExITemp[i][1][k];
	//      vectorY[i][1][k] = injFieldsBottom->EyITemp[i][1][k];
	//      vectorZ[i][1][k] = injFieldsBottom->EzITemp[i][1][k];

	vectorX[i][0][k] = injFieldsBottom->ExITemp[i][0][k];
	vectorY[i][0][k] = injFieldsBottom->EyITemp[i][0][k];
	vectorZ[i][0][k] = injFieldsBottom->EzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//      vectorX[i][ny-2][k] = injFieldsTop->ExITemp[i][ny-2][k];
	//      vectorY[i][ny-2][k] = injFieldsTop->EyITemp[i][ny-2][k];
	//      vectorZ[i][ny-2][k] = injFieldsTop->EzITemp[i][ny-2][k];

	vectorX[i][ny-1][k] = injFieldsTop->ExITemp[i][ny-1][k];
	vectorY[i][ny-1][k] = injFieldsTop->EyITemp[i][ny-1][k];
	vectorZ[i][ny-1][k] = injFieldsTop->EzITemp[i][ny-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//      vectorX[i][j][1] = injFieldsRear->ExITemp[i][j][1];
	//      vectorY[i][j][1] = injFieldsRear->EyITemp[i][j][1];
	//      vectorZ[i][j][1] = injFieldsRear->EzITemp[i][j][1];

	vectorX[i][j][0] = injFieldsRear->ExITemp[i][j][0];
	vectorY[i][j][0] = injFieldsRear->EyITemp[i][j][0];
	vectorZ[i][j][0] = injFieldsRear->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//      vectorX[i][j][nz-2] = injFieldsFront->ExITemp[i][j][nz-2];
	//      vectorY[i][j][nz-2] = injFieldsFront->EyITemp[i][j][nz-2];
	//      vectorZ[i][j][nz-2] = injFieldsFront->EzITemp[i][j][nz-2];

	vectorX[i][j][nz-1] = injFieldsFront->ExITemp[i][j][nz-1];
	vectorY[i][j][nz-1] = injFieldsFront->EyITemp[i][j][nz-1];
	vectorZ[i][j][nz-1] = injFieldsFront->EzITemp[i][j][nz-1];
      }
  }
}

/*! get Potential array ** */
double ***EMfields3D::getPHI() {
  return (PHI);
}
/*! get Ex(X,Y,Z) */
double &EMfields3D::getEx(int indexX, int indexY, int indexZ) const {
  return (Ex[indexX][indexY][indexZ]);
}
/*! get Electric field component X array */
double ***EMfields3D::getEx() {
  return (Ex);
}
double ***EMfields3D::getExth() {
  return (Exth);
}
/*! get Electric Field component X array cell without the ghost cells */
double ***EMfields3D::getExc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Ex[i][j][k]     + 
			       Ex[i + 1][j][k] + 
			       Ex[i][j + 1][k] + 
			       Ex[i][j][k + 1] + 
			       Ex[i + 1][j + 1][k] + 
			       Ex[i + 1][j][k + 1] + 
			       Ex[i][j + 1][k + 1] + 
			       Ex[i + 1][j + 1][k + 1]);

  return arr;
}
/*! get Ey(X,Y,Z) */
double &EMfields3D::getEy(int indexX, int indexY, int indexZ) const {
  return (Ey[indexX][indexY][indexZ]);
}
/*! get Electric field component Y array */
double ***EMfields3D::getEy() {
  return (Ey);
}
double ***EMfields3D::getEyth() {
  return (Eyth);
}
/*! get Electric Field component Y array cell without the ghost cells */
double ***EMfields3D::getEyc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Ey[i][j][k]     + 
			       Ey[i + 1][j][k] + 
			       Ey[i][j + 1][k] + 
			       Ey[i][j][k + 1] + 
			       Ey[i + 1][j + 1][k] + 
			       Ey[i + 1][j][k + 1] + 
			       Ey[i][j + 1][k + 1] + 
			       Ey[i + 1][j + 1][k + 1]);
  return arr;
}
/*! get Ez(X,Y,Z) */
double &EMfields3D::getEz(int indexX, int indexY, int indexZ) const {
  return (Ez[indexX][indexY][indexZ]);
}
/*! get Electric field component Z array */
double ***EMfields3D::getEz() {
  return (Ez);
}
double ***EMfields3D::getEzth() {
  return (Ezth);
}
/*! get Electric Field component Z array cell without the ghost cells */
double ***EMfields3D::getEzc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Ez[i][j][k]     +
			       Ez[i + 1][j][k] +
			       Ez[i][j + 1][k] +
			       Ez[i][j][k + 1] +
			       Ez[i + 1][j + 1][k] +
			       Ez[i + 1][j][k + 1] +
			       Ez[i][j + 1][k + 1] +
			       Ez[i + 1][j + 1][k + 1]);

  return arr;
}
/*! get Bx(X,Y,Z) */
double &EMfields3D::getBx(int indexX, int indexY, int indexZ) const {
  return (Bxn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component X array */
double ***EMfields3D::getBx() {
  return (Bxn);
}
/*! get Magnetic Field component X array cell without the ghost cells */
double EMfields3D::getBxc(int i, int j, int k) {
  return (Bxc[i][j][k]);
}
/*! get Magnetic Field component X array cell without the ghost cells */
double ***EMfields3D::getBxc() {
  return (Bxc);
}
/*! get By(X,Y,Z) */
double &EMfields3D::getBy(int indexX, int indexY, int indexZ) const {
  return (Byn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Y array */
double ***EMfields3D::getBy() {
  return (Byn);
}
/*! get Magnetic Field component Y array cell without the ghost cells */
double EMfields3D::getByc(int i, int j, int k) {
  return (Byc[i][j][k]);
}
/*! get Magnetic Field component Y array cell without the ghost cells */
double ***EMfields3D::getByc() {
  return (Byc);
}
/*! get Bz(X,Y,Z) */
double &EMfields3D::getBz(int indexX, int indexY, int indexZ) const {
  return (Bzn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Z array */
double ***EMfields3D::getBz() {
  return (Bzn);
}
/*! get Magnetic Field component Z array cell without the ghost cells */
double EMfields3D::getBzc(int i, int j, int k) {
  return (Bzc[i][j][k]);
}
/*! get Magnetic Field component Z array cell without the ghost cells */
double ***EMfields3D::getBzc() {
  return (Bzc);
}
/*! get rhoc(X,Y,Z) */
double &EMfields3D::getRHOc(int indexX, int indexY, int indexZ) const {
  return (rhoc[indexX][indexY][indexZ]);
} double ***EMfields3D::getRHOc() {
  return (rhoc);
}
/*! get density on node(indexX,indexY,indexZ) */
double &EMfields3D::getRHOn(int indexX, int indexY, int indexZ) const {
  return (rhon[indexX][indexY][indexZ]);
}
/*! get density array defined on nodes */
double ***EMfields3D::getRHOn() {
  return (rhon);
}
/*! get rhos(X,Y,Z) : density for species */
double &EMfields3D::getRHOns(int indexX, int indexY, int indexZ, int is) const {
  return (rhons[is][indexX][indexY][indexZ]);
}
/*! SPECIES: get density array defined on center cells */
double &EMfields3D::getRHOcs(int indexX, int indexY, int indexZ, int is) const {
  return (rhocs[is][indexX][indexY][indexZ]);
}
/*! SPECIES: get density array defined on center cells */
double ****& EMfields3D::getRHOcs() {
  return (rhocs);
}
/*! get density array defined on nodes */
double ****EMfields3D::getRHOns() {
  return (rhons);
}
/*! get density array defined on nodes for one species */
double ***& EMfields3D::getRHOns(int is) {
  return (rhons[is]);
}
/*! get species density component X array cell without the ghost cells */
double ***EMfields3D::getRHOcs(int is) {
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (rhons[is][i][j][k]     +
			       rhons[is][i + 1][j][k] +
			       rhons[is][i][j + 1][k] +
			       rhons[is][i][j][k + 1] +
			       rhons[is][i + 1][j + 1][k] +
			       rhons[is][i + 1][j][k + 1] +
			       rhons[is][i][j + 1][k + 1] +
			       rhons[is][i + 1][j + 1][k + 1]);

  return arr;
}

double ***& EMfields3D::getRHOcs(int is, int dummy) {
  return (rhocs[is]);
}

/*! get Bx_ext(X,Y,Z)  */
double &EMfields3D::getBx_ext(int indexX, int indexY, int indexZ) const{
  return(Bx_ext[indexX][indexY][indexZ]);
}
/*!  get By_ext(X,Y,Z) */
double &EMfields3D::getBy_ext(int indexX, int indexY, int indexZ) const{
  return(By_ext[indexX][indexY][indexZ]);
}
/*!  get Bz_ext(X,Y,Z) */
double &EMfields3D::getBz_ext(int indexX, int indexY, int indexZ) const{
  return(Bz_ext[indexX][indexY][indexZ]);
}

/*! get Bx_ext  */
double ***EMfields3D::getBx_ext() {
  return(Bx_ext);
}
/*!  get By_ext */
double ***EMfields3D::getBy_ext() {
  return(By_ext);
}
/*!  get Bz_ext */
double ***EMfields3D::getBz_ext() {
  return(Bz_ext);
}

double ***&EMfields3D::getBxTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
	arr[i][j][k] = Bxn[i][j][k] + Fext * Bx_ext[i][j][k];

  return arr;
}
double ***&EMfields3D::getByTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
	arr[i][j][k] = Byn[i][j][k] + Fext * By_ext[i][j][k];

  return arr;
}
double ***&EMfields3D::getBzTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
	arr[i][j][k] = Bzn[i][j][k] + Fext * Bz_ext[i][j][k];

  return arr;
}

/*! SPECIES: get pressure tensor component XX defined on nodes */
double ****EMfields3D::getpXXsn() {
  return (pXXsn);
}
/*! SPECIES: get pressure tensor component XY defined on nodes */
double ****EMfields3D::getpXYsn() {
  return (pXYsn);
}
/*! SPECIES: get pressure tensor component XZ defined on nodes */
double ****EMfields3D::getpXZsn() {
  return (pXZsn);
}
/*! SPECIES: get pressure tensor component YY defined on nodes */
double ****EMfields3D::getpYYsn() {
  return (pYYsn);
}
/*! SPECIES: get pressure tensor component YZ defined on nodes */
double ****EMfields3D::getpYZsn() {
  return (pYZsn);
}
/*! SPECIES: get pressure tensor component ZZ defined on nodes */
double ****EMfields3D::getpZZsn() {
  return (pZZsn);
}
/*! get current -Direction X */
double &EMfields3D::getJx(int indexX, int indexY, int indexZ) const {
  return (Jx[indexX][indexY][indexZ]);
}
/*! get current array X component * */
double ***EMfields3D::getJx() {
  return (Jx);
}
/*! get current -Direction Y */
double &EMfields3D::getJy(int indexX, int indexY, int indexZ) const {
  return (Jy[indexX][indexY][indexZ]);
}
/*! get current array Y component * */
double ***EMfields3D::getJy() {
  return (Jy);
}
/*! get current -Direction Z */
double &EMfields3D::getJz(int indexX, int indexY, int indexZ) const {
  return (Jz[indexX][indexY][indexZ]);
}
/*! get current array Z component * */
double ***EMfields3D::getJz() {
  return (Jz);
}
/*!SPECIES: get current array X component */
double ****EMfields3D::getJxs() {
  return (Jxs);
}
double ***& EMfields3D::getJxs(int is) {
  return (Jxs[is]);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJxs(int indexX, int indexY, int indexZ, int is) const {
  return (Jxs[is][indexX][indexY][indexZ]);
}
/*! get Magnetic Field component X array species is cell without the ghost cells */
double ***EMfields3D::getJxsc(int is) {
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Jxs[is][i][j][k]     +
			       Jxs[is][i + 1][j][k] +
			       Jxs[is][i][j + 1][k] +
			       Jxs[is][i][j][k + 1] +
			       Jxs[is][i + 1][j + 1][k] +
			       Jxs[is][i + 1][j][k + 1] +
			       Jxs[is][i][j + 1][k + 1] +
			       Jxs[is][i + 1][j + 1][k + 1]);

  return arr;
}
/*! SPECIES: get current array Y component */
double ****EMfields3D::getJys() {
  return (Jys);
}
double ***& EMfields3D::getJys(int is) {
  return (Jys[is]);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJys(int indexX, int indexY, int indexZ, int is) const {
  return (Jys[is][indexX][indexY][indexZ]);
}
/*! get current component Y array species is cell without the ghost cells */
double ***EMfields3D::getJysc(int is) {
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Jys[is][i][j][k]     + 
			       Jys[is][i + 1][j][k] +
			       Jys[is][i][j + 1][k] +
			       Jys[is][i][j][k + 1] +
			       Jys[is][i + 1][j + 1][k] +
			       Jys[is][i + 1][j][k + 1] +
			       Jys[is][i][j + 1][k + 1] +
			       Jys[is][i + 1][j + 1][k + 1]);

  return arr;
}
/*!SPECIES: get current array Z component */
double ****EMfields3D::getJzs() {
  return (Jzs);
}
double ***& EMfields3D::getJzs(int is) {
  return (Jzs[is]);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJzs(int indexX, int indexY, int indexZ, int is) const {
  return (Jzs[is][indexX][indexY][indexZ]);
}
/*! get current component Z array species is cell without the ghost cells */
double ***EMfields3D::getJzsc(int is) {
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
	arr[i][j][k] = .125 * (Jzs[is][i][j][k]     + 
			       Jzs[is][i + 1][j][k] +
			       Jzs[is][i][j + 1][k] +
			       Jzs[is][i][j][k + 1] +
			       Jzs[is][i + 1][j + 1][k] +
			       Jzs[is][i + 1][j][k + 1] +
			       Jzs[is][i][j + 1][k + 1] +
			       Jzs[is][i + 1][j + 1][k + 1]);

  return arr;
}

double EMfields3D::getRHOINIT(int is, int i, int j, int k){
  return (RHOINIT[is][i][j][k]);
}

/*! get the electric field energy */
/*! mlmd: i need the communicator also
  double EMfields3D::getEenergy(void) { */
double EMfields3D::getEenergy(MPI_Comm Comm) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++)
	localEenergy += .5 * dx * dy * dz * (Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI);

  MPI_Allreduce(&localEenergy, &totalEenergy, 1, MPI_DOUBLE, MPI_SUM, Comm);
  return (totalEenergy);

}
/*! get the magnetic field energy */
/*! mlmd: i need the communicator also
  double EMfields3D::getBenergy(void) { */
double EMfields3D::getBenergy(MPI_Comm Comm) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++){
	Bxt = Bxn[i][j][k]+Fext*Bx_ext[i][j][k];
	Byt = Byn[i][j][k]+Fext*By_ext[i][j][k];
	Bzt = Bzn[i][j][k]+Fext*Bz_ext[i][j][k];
	localBenergy += .5*dx*dy*dz*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, Comm);
  return (totalBenergy);
}


/*! Print info about electromagnetic field */
void EMfields3D::print(void) const {
}

/*! destructor: deallocate arrays */
EMfields3D::~EMfields3D() {
  // nodes
  delArr3(Ex, nxn, nyn);
  delArr3(Ey, nxn, nyn);
  delArr3(Ez, nxn, nyn);
  delArr3(Exth, nxn, nyn);
  delArr3(Eyth, nxn, nyn);
  delArr3(Ezth, nxn, nyn);
  delArr3(Bxn, nxn, nyn);
  delArr3(Byn, nxn, nyn);
  delArr3(Bzn, nxn, nyn);
  delArr3(rhon, nxn, nyn);
  delArr3(Jx, nxn, nyn);
  delArr3(Jy, nxn, nyn);
  delArr3(Jz, nxn, nyn);
  delArr3(Jxh, nxn, nyn);
  delArr3(Jyh, nxn, nyn);
  delArr3(Jzh, nxn, nyn);
  // nodes and species
  delArr4(rhons, ns, nxn, nyn);
  delArr4(Jxs, ns, nxn, nyn);
  delArr4(Jys, ns, nxn, nyn);
  delArr4(Jzs, ns, nxn, nyn);
  delArr4(pXXsn, ns, nxn, nyn);
  delArr4(pXYsn, ns, nxn, nyn);
  delArr4(pXZsn, ns, nxn, nyn);
  delArr4(pYYsn, ns, nxn, nyn);
  delArr4(pYZsn, ns, nxn, nyn);
  delArr4(pZZsn, ns, nxn, nyn);
  // central points
  delArr3(PHI, nxc, nyc);
  delArr3(Bxc, nxc, nyc);
  delArr3(Byc, nxc, nyc);
  delArr3(Bzc, nxc, nyc);
  delArr3(rhoc, nxc, nyc);
  delArr3(rhoh, nxc, nyc);
  delArr4(rhocs, ns, nxc, nyc);
  // various stuff needs to be deallocated too
  delArr3(tempXC, nxc, nyc);
  delArr3(tempYC, nxc, nyc);
  delArr3(tempZC, nxc, nyc);
  delArr3(tempXN, nxn, nyn);
  delArr3(tempYN, nxn, nyn);
  delArr3(tempZN, nxn, nyn);
  delArr3(tempC, nxc, nyc);
  delArr3(tempX, nxn, nyn);
  delArr3(tempY, nxn, nyn);
  delArr3(tempZ, nxn, nyn);
  delArr3(temp2X, nxn, nyn);
  delArr3(temp2Y, nxn, nyn);
  delArr3(temp2Z, nxn, nyn);
  delArr3(imageX, nxn, nyn);
  delArr3(imageY, nxn, nyn);
  delArr3(imageZ, nxn, nyn);
  delArr3(Dx, nxn, nyn);
  delArr3(Dy, nxn, nyn);
  delArr3(Dz, nxn, nyn);
  delArr3(vectX, nxn, nyn);
  delArr3(vectY, nxn, nyn);
  delArr3(vectZ, nxn, nyn);
  delArr3(divC, nxc, nyc);

  delete[]rhoINIT;
  delete[]rhoINJECT;
  delete[]DriftSpecies;

  delete[]dx_Ch;
  delete[]dy_Ch;
  delete[]dz_Ch;

  delArr4(RHOINIT, ns, nxc, nyc);
  
  delete injFieldsRight;
  delete injFieldsLeft;
  delete injFieldsTop;
  delete injFieldsBottom;
  delete injFieldsFront;
  delete injFieldsRear;

  delete []qom;

  delArr3(imageEX, nxc, nyc);
  delArr3(imageEY, nxc, nyc);
  delArr3(imageEZ, nxc, nyc);
  delArr3(imageBX, nxn, nyn);
  delArr3(imageBY, nxn, nyn);
  delArr3(imageBZ, nxn, nyn);

  delArr3(arr, nxn, nyn);
  delArr3(Lambda, nxn, nyn);

  if (CommToParent_InDel != MPI_COMM_NULL and MLMD_BC){
    delete[]RGBC_Info_Ghost;
    delete[]RGBC_Info_Active;

    if (MLMD_BCBufferArea){
      delete[]RGBC_Info_Buffer;
      delete[]DirBuffer;
    }

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      delete[]RGBC_Info_fix3B; 
    }

    delArr2(Ex_Active_BC, RG_numBCMessages_Active);
    delArr2(Ey_Active_BC, RG_numBCMessages_Active);
    delArr2(Ez_Active_BC, RG_numBCMessages_Active);

    delArr2(Ex_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Ey_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Ez_Ghost_BC, RG_numBCMessages_Ghost);

    delArr2(Exth_Active_BC, RG_numBCMessages_Active);
    delArr2(Eyth_Active_BC, RG_numBCMessages_Active);
    delArr2(Ezth_Active_BC, RG_numBCMessages_Active);

    delArr2(Exth_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Eyth_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Ezth_Ghost_BC, RG_numBCMessages_Ghost);

    delArr2(Bxn_Active_BC, RG_numBCMessages_Active);
    delArr2(Byn_Active_BC, RG_numBCMessages_Active);
    delArr2(Bzn_Active_BC, RG_numBCMessages_Active);

    delArr2(Bxn_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Byn_Ghost_BC, RG_numBCMessages_Ghost);
    delArr2(Bzn_Ghost_BC, RG_numBCMessages_Ghost);

    if (MLMD_BCBufferArea){
    
      delArr2(Ex_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Ey_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Ez_Buffer_BC, RG_numBCMessages_Buffer);
      
      delArr2(Exth_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Eyth_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Ezth_Buffer_BC, RG_numBCMessages_Buffer);
      
      delArr2(Bxn_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Byn_Buffer_BC, RG_numBCMessages_Buffer);
      delArr2(Bzn_Buffer_BC, RG_numBCMessages_Buffer);
    }// end if (MLMD_BCBufferArea){

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      delArr2(Bxc_fix3B_BC, RG_numBCMessages_fix3B);
      delArr2(Byc_fix3B_BC, RG_numBCMessages_fix3B);
      delArr2(Bzc_fix3B_BC, RG_numBCMessages_fix3B);
    } // end if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){

    delete[] RGMsg; 

    if (MLMD_BCBufferArea){
      delete[]RGMsgBuffer;
    }
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      delete[]RGMsgfix3B;
    }
    
  } // end if (CommToParent_InDel != MPI_COMM_NULL){
 
  if (numChildren >0 and MLMD_BC ){   

    delArr2(CG_Info_Ghost, numChildren);
    delete[]CG_numBCMessages_Ghost;
    delArr2(CG_Info_Active, numChildren);
    delete[]CG_numBCMessages_Active;
    delete[]CGMsg;   

    if (MLMD_BCBufferArea){
      delArr2(CG_Info_Buffer, numChildren);
      delete[]CG_numBCMessages_Buffer;  
      delete[]CGMsgBuffer; 
    }

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      delArr2(CG_Info_Fix3B, numChildren);
      delete[]CG_numBCMessages_Fix3B;
      delete[]CGMsgFix3B; 
    }

  } // end if (numChildren >0){

  if (BSNeeded){
      delArr3(Ex_BS, nxn, nyn);
      delArr3(Ey_BS, nxn, nyn);
      delArr3(Ez_BS, nxn, nyn);
    }

  /* deallocate datatypes */
  MPI_Type_free(&MPI_RGBC_struct);
  /* end deallocate datatypes */
}

/*! mlmd specific functions */

void EMfields3D::initWeightBC(Grid *grid, VirtualTopology3D *vct){

  bool VerboseCheck= false;
  

  /* what I am doing here
     phase 1: the refined grids calculates:
     -- how many points they need from the CGs
     -- from which cores in the CG
     -- side (bottom, top, left, right, front, bottom)
     -- CG coordinates (x,y,z) of the first point 


     phase 2: message with previous info is sent to the CG
     -- phase 2a: all cores in the RG send their message structure to core 0, local grid
     NB: send one message more with -1 in RG_core
     -- phase 2b: core 0, local grid, prepares and sends message for each core of the coarse grid, with info to build their weights
     -- phase 2c: all CG cores receive message and appropriately react building their weights
      
     phase 3: the coarse grid builds its weights, based on previous infos
  */
  /* Jun 2017: I include the possibility of setting E BC for a buffer length of RF -- 
     to do so, MLMD_BCBufferArea = true */
  MPI_Status status;

  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  // rank on the local grid
  int rank_local= vct->getCartesian_rank();

  // rank as a child in the parent-child communicator
  int rank_As_Child=-1;
  MPI_Comm CommToParent= vct->getCommToParent(); // if != MPI_COMM_NULL the grid is a child
  if (CommToParent != MPI_COMM_NULL){
    MPI_Comm_rank(CommToParent, &rank_As_Child);
  }

  /* rank as a parent in the parent-child communicator;
     since one grid may have more than one child, I have to calculate every time
     on the right communicatore */
  int rank_As_Parent;

  RG_numBCMessages_Ghost= 0;
  RG_numBCMessages_Active= 0;
  RG_numBCMessages_Buffer= 0;
  RG_numBCMessages_fix3B= 0;

  /* here, vectors where core 0 of the local child grid assembles the messages
     from all the local grid 
     de-allocated at the end of the function */
  RGBC_struct * RGBC_Info_Ghost_LevelWide;
  int RG_numBCMessages_Ghost_LevelWide=0;

  RGBC_struct * RGBC_Info_Active_LevelWide;
  int RG_numBCMessages_Active_LevelWide=0;

  RGBC_struct * RGBC_Info_Buffer_LevelWide;
  int RG_numBCMessages_Buffer_LevelWide=0;
  
  RGBC_struct * RGBC_Info_fix3B_LevelWide;
  int RG_numBCMessages_fix3B_LevelWide=0;

  /* phase 1 */
  // as a child
  if (CommToParent != MPI_COMM_NULL){

    RG_MaxMsgSize=0; // value is calculated in the initWeightBC_Phase1's
    /* instantiate ghost structure */
    RGBC_Info_Ghost = new RGBC_struct[MAX_RG_numBCMessages];
    RGBC_Info_Active = new RGBC_struct[MAX_RG_numBCMessages];

    RG_MaxMsgBufferSize= 0; // value is calculated in the initWeightBCBuffer_Phase1's 
    if (MLMD_BCBufferArea){
      RGBC_Info_Buffer = new RGBC_struct[MAX_RG_numBCMessages];
      DirBuffer = new char[MAX_RG_numBCMessages];
    }

    RG_Maxfix3BMsgSize= 0; // value is calculated in the initWeightBC_Phase1's   
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      RGBC_Info_fix3B = new RGBC_struct[MAX_RG_numBCMessages];
    }

    // -1 for ghost
    // 0 for active

   
    // ghost
    initWeightBC_Phase1(grid, vct, RGBC_Info_Ghost, &RG_numBCMessages_Ghost, -1);
    // active
    initWeightBC_Phase1(grid, vct, RGBC_Info_Active, &RG_numBCMessages_Active, 0);
    // buffer
    if (MLMD_BCBufferArea){
      initWeightBCBuffer_Phase1(grid, vct, RGBC_Info_Buffer, &RG_numBCMessages_Buffer, &RG_MaxMsgBufferSize);
    }
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      initWeightBCfix3B_Phase1(grid, vct, RGBC_Info_fix3B, &RG_numBCMessages_fix3B, &RG_Maxfix3BMsgSize);
    }
    
    // now I have all the info to initialise the BC vectors

    // rows: [0 - RG_numBCMessages_Active]                                                    
    // columns: [0 - RG_MaxMsgSize]                                                          
    Ex_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Ey_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Ez_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);

    Exth_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Eyth_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Ezth_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
     
    Bxn_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Byn_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);
    Bzn_Active_BC= newArr2(double, RG_numBCMessages_Active, RG_MaxMsgSize);

    // rows: [0 - RG_numBCMessages_Ghost]                                                                                                 
    // columns: [0 - RG_MaxMsgSize]                                                                            
    Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Ey_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Ez_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);

    Exth_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Eyth_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Ezth_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);

    Bxn_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Byn_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);
    Bzn_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);

    // instantiate this only if needed
    if (MLMD_BCBufferArea){
      // rows: [0 - RG_numBCMessages_Buffer] 
      // columns: [0 - RG_MaxMsgBufferSize]
      Ex_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Ey_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Ez_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      
      Exth_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Eyth_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Ezth_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      
      Bxn_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Byn_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      Bzn_Buffer_BC= newArr2(double, RG_numBCMessages_Buffer, RG_MaxMsgBufferSize);
      
    }

    // instantiate only if needed
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      Bxc_fix3B_BC= newArr2(double, RG_numBCMessages_fix3B, RG_Maxfix3BMsgSize);
      Byc_fix3B_BC= newArr2(double, RG_numBCMessages_fix3B, RG_Maxfix3BMsgSize);
      Bzc_fix3B_BC= newArr2(double, RG_numBCMessages_fix3B, RG_Maxfix3BMsgSize);
    }

    // this used when receiving
    // active and ghost
    RGMsg= new double[RG_MaxMsgSize *NumF];
    
    if (MLMD_BCBufferArea){
      RGMsgBuffer= new double[RG_MaxMsgBufferSize *NumF];
    }

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      RGMsgfix3B= new double[RG_Maxfix3BMsgSize *Numfix3B];
    }

    /**** check starts ****/
    // before proceeding, a check with the possibility of aborting if the check is failed
    int PG= vct->getParentGridNum();
    int localRank= vct->getCartesian_rank();
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double MsgLimsXMin, MsgLimsXMax, MsgLimsYMin, MsgLimsYMax, MsgLimsZMin, MsgLimsZMax;
    // check on ghost
    for (int i=0; i< RG_numBCMessages_Ghost; i++){
      int CG= RGBC_Info_Ghost[i].CG_core;
       
      MsgLimsXMin= RGBC_Info_Ghost[i].CG_x_first;
      MsgLimsXMax= RGBC_Info_Ghost[i].CG_x_first+ dx*(RGBC_Info_Ghost[i].np_x-1);

      MsgLimsYMin= RGBC_Info_Ghost[i].CG_y_first;
      MsgLimsYMax= RGBC_Info_Ghost[i].CG_y_first+ dy*(RGBC_Info_Ghost[i].np_y-1);

      MsgLimsZMin= RGBC_Info_Ghost[i].CG_z_first;
      MsgLimsZMax= RGBC_Info_Ghost[i].CG_z_first+ dz*(RGBC_Info_Ghost[i].np_z-1);

      grid->getParentLimits(vct, CG, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

      if (MsgLimsXMin < xmin or MsgLimsXMax > xmax or MsgLimsYMin < ymin or MsgLimsYMax > ymax or MsgLimsZMin < zmin or MsgLimsZMax > zmax){
	cout <<"G" <<numGrid <<"R" << localRank <<" Msg " << i << " we have a problem in initWeightBC, ghost BC, aborting ... " << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    } // for (int i=0; i< RG_numBCMessages_Ghost; i++)

    // check on active
    for (int i=0; i< RG_numBCMessages_Active; i++){
      int CG= RGBC_Info_Active[i].CG_core;
       
      MsgLimsXMin= RGBC_Info_Active[i].CG_x_first;
      MsgLimsXMax= RGBC_Info_Active[i].CG_x_first+ dx*(RGBC_Info_Active[i].np_x-1);

      MsgLimsYMin= RGBC_Info_Active[i].CG_y_first;
      MsgLimsYMax= RGBC_Info_Active[i].CG_y_first+ dy*(RGBC_Info_Active[i].np_y-1);

      MsgLimsZMin= RGBC_Info_Active[i].CG_z_first;
      MsgLimsZMax= RGBC_Info_Active[i].CG_z_first+ dz*(RGBC_Info_Active[i].np_z-1);

      grid->getParentLimits(vct, CG, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

      if (MsgLimsXMin < xmin or MsgLimsXMax > xmax or MsgLimsYMin < ymin or MsgLimsYMax > ymax or MsgLimsZMin < zmin or MsgLimsZMax > zmax){
	cout <<"G" <<numGrid <<"R" << localRank <<" Msg " << i << " we have a problem in initWeightBC, active BC, aborting ... " << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    // check on buffer
    if (MLMD_BCBufferArea){
      for (int i=0; i< RG_numBCMessages_Buffer; i++){
	int CG= RGBC_Info_Buffer[i].CG_core;
       
	MsgLimsXMin= RGBC_Info_Buffer[i].CG_x_first;
	MsgLimsXMax= RGBC_Info_Buffer[i].CG_x_first+ dx*(RGBC_Info_Buffer[i].np_x-1);

	MsgLimsYMin= RGBC_Info_Buffer[i].CG_y_first;
	MsgLimsYMax= RGBC_Info_Buffer[i].CG_y_first+ dy*(RGBC_Info_Buffer[i].np_y-1);
	
	MsgLimsZMin= RGBC_Info_Buffer[i].CG_z_first;
	MsgLimsZMax= RGBC_Info_Buffer[i].CG_z_first+ dz*(RGBC_Info_Buffer[i].np_z-1);

	grid->getParentLimits(vct, CG, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

	if (MsgLimsXMin < xmin or MsgLimsXMax > xmax or MsgLimsYMin < ymin or MsgLimsYMax > ymax or MsgLimsZMin < zmin or MsgLimsZMax > zmax){
	  cout <<"G" <<numGrid <<"R" << localRank <<" Msg " << i << " we have a problem in initWeightBC, buffer BC, aborting ... " << endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);
	}
      }
    }

    /**** check ends ****/
  } // end if (numGrid>0)
  /* end phase 1 */


  /*MPI_Barrier (vct->getComm());
  if (rank_local==0){
    cout << "Grid " << numGrid << " is after phase 1 of initWeight_BC" << endl;
    }*/

  /* phase 2, RG sends info to CG */
  // children grids
  if (CommToParent != MPI_COMM_NULL ){  

    // phase 2a: children assemble all the info in core 0 in the LOCAL child grid, 0-> XLEN*YLEN*ZLEN-1 
    if (rank_local > 0){
      /* send one message more; the last message has -1 in the RG_core
	 to signal end of 'valid' messages */
      /* -1 as tag for ghost, 0 as tag for active */

      // ghost
      MPI_Send(RGBC_Info_Ghost, RG_numBCMessages_Ghost+1, MPI_RGBC_struct, 0, TAG_BC_GHOST, vct->getComm());

      /*cout << "Rank local " << rank_local << ", grid " <<numGrid << "RG_numBCMessages_Ghost: " << RG_numBCMessages_Ghost << "+1" <<endl;

	for (int m=0; m< RG_numBCMessages_Ghost+1 ; m++){
	cout <<"Rank local " << rank_local << ", grid " <<numGrid << ", m: " << m <<", RGBC_Info_Ghost[m].RG_core: " << RGBC_Info_Ghost[m].RG_core << ", RGBC_Info_Ghost[m].CG_core: " << RGBC_Info_Ghost[m].CG_core<< endl;
	} */
      //active
      MPI_Send(RGBC_Info_Active, RG_numBCMessages_Active+1, MPI_RGBC_struct, 0, TAG_BC_ACTIVE, vct->getComm());

      /*cout << "ACTIVE: Rank local " << rank_local << ", grid " <<numGrid << "RG_numBCMessages_Active: " << RG_numBCMessages_Active << "+1" <<endl;

	for (int m=0; m< RG_numBCMessages_Active+1 ; m++){
	cout <<"ACTIVE: Rank local " << rank_local << ", grid " <<numGrid << ", m: " << m <<", RGBC_Info_Active[m].RG_core: " << RGBC_Info_Active[m].RG_core << ", RGBC_Info_Active[m].CG_core: " << RGBC_Info_Active[m].CG_core<< endl;
	}*/

      if (MLMD_BCBufferArea){
	MPI_Send(RGBC_Info_Buffer, RG_numBCMessages_Buffer+1, MPI_RGBC_struct, 0, TAG_BC_BUFFER, vct->getComm());
      }

      if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
	MPI_Send(RGBC_Info_fix3B, RG_numBCMessages_fix3B+1, MPI_RGBC_struct, 0, TAG_BC_FIX3B, vct->getComm());
      }
      
    }

    /*MPI_Barrier (vct->getComm());
    if (rank_local==0){
      cout << "Grid " << numGrid << " is after phase 2-prel of initWeight_BC (only for RGs)" << endl;
      }*/

   
    if  (rank_local==0){

      // only core 0 local needs to instantiate this 
      RGBC_Info_Ghost_LevelWide = new RGBC_struct[MAX_size_LevelWide];
      initWeightBC_Phase2a(grid, vct, RGBC_Info_Ghost_LevelWide, &RG_numBCMessages_Ghost_LevelWide, RGBC_Info_Ghost, RG_numBCMessages_Ghost, -1);

      RGBC_Info_Active_LevelWide = new RGBC_struct[MAX_size_LevelWide];
      initWeightBC_Phase2a(grid, vct, RGBC_Info_Active_LevelWide, &RG_numBCMessages_Active_LevelWide, RGBC_Info_Active, RG_numBCMessages_Active, 0);
       
      if (MLMD_BCBufferArea){
	RGBC_Info_Buffer_LevelWide = new RGBC_struct[MAX_size_LevelWide];
	initWeightBC_Phase2a(grid, vct, RGBC_Info_Buffer_LevelWide, &RG_numBCMessages_Buffer_LevelWide, RGBC_Info_Buffer, RG_numBCMessages_Buffer, 7);
      }

      if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
	RGBC_Info_fix3B_LevelWide = new RGBC_struct[MAX_size_LevelWide];
	initWeightBC_Phase2a(grid, vct, RGBC_Info_fix3B_LevelWide, &RG_numBCMessages_fix3B_LevelWide, RGBC_Info_fix3B, RG_numBCMessages_fix3B, 9);
      }
      
      /*cout << "Core 0 of grid " << numGrid <<": I have " << RG_numBCMessages_Ghost_LevelWide << " messages: " <<endl;               	
	for (int m=0; m< RG_numBCMessages_Ghost_LevelWide; m++){                                              
	cout << "Message " << m << " from RG core " << RGBC_Info_Ghost_LevelWide[m].RG_core << " to CG core" << RGBC_Info_Ghost_LevelWide[m].CG_core << " on Parent-Child communicator "<<endl;                         
	}

	cout << "ACTIVE Core 0 of grid " << numGrid <<": I have " << RG_numBCMessages_Active_LevelWide << " messages: " <<endl;
	for (int m=0; m< RG_numBCMessages_Active_LevelWide; m++){       
	cout << "ACTIVE message " << m << " from RG core " << RGBC_Info_Active_LevelWide[m].RG_core << " to CG core" << RGBC_Info_Active_LevelWide[m].CG_core << " on Parent-Child communicator "<<endl;                         
	}*/

      /*cout << "BUFFER Core 0 of grid " << numGrid <<": I have " << RG_numBCMessages_Buffer_LevelWide << " messages: " <<endl;
      for (int m=0; m< RG_numBCMessages_Buffer_LevelWide; m++){       
	cout << "BUFFER message " << m << " from RG core " << RGBC_Info_Buffer_LevelWide[m].RG_core << " to CG core" << RGBC_Info_Buffer_LevelWide[m].CG_core << " on Parent-Child communicator "<<endl;                         
	}*/

      
    } // end if (rank_local==0)
   
    /*MPI_Barrier (vct->getComm());
    if (rank_local==0){
      cout << "Grid " << numGrid << " is after phase 2a of initWeight_BC (only for RGs)" << endl;
      } */


    // phase 2b: core 0 of the child grid assembles & sends messages for all CG cores
    if (rank_local==0){
      //ghost
      initWeightBC_Phase2b(grid, vct, RGBC_Info_Ghost_LevelWide, RG_numBCMessages_Ghost_LevelWide, -1);
      // active
      initWeightBC_Phase2b(grid, vct, RGBC_Info_Active_LevelWide, RG_numBCMessages_Active_LevelWide, 0);
      //buffer
      if (MLMD_BCBufferArea){
	initWeightBC_Phase2b(grid, vct, RGBC_Info_Buffer_LevelWide, RG_numBCMessages_Buffer_LevelWide, 7);
      }
      if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
	initWeightBC_Phase2b(grid, vct, RGBC_Info_fix3B_LevelWide, RG_numBCMessages_fix3B_LevelWide, 9);
      }

    } // end if (rank_local==0), phase 2b
  } // end if on children grid
  

  /*MPI_Barrier(vct->getComm());
  if (rank_local==0){
    cout << "I am grid " << numGrid <<", I am after the barrier marking Phase2b" << endl;
    }*/
  
  /*MPI_Barrier(MPI_COMM_WORLD);
  if (vct->getSystemWide_rank()==0){
    cout << "After barrier phase 2" << endl;
    }*/

  // phase 2c, only for the parent
  // each CG core will receive a message per child grid
  if (numChildren > 0 ){

    // ghost & active
    CG_MaxSizeMsg=0;
    // buffer
    CG_MaxSizeBufferMsg=0;
    // fix3B
    CG_MaxSizeFix3BMsg=0;

    // ghost
    CG_Info_Ghost= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
    CG_numBCMessages_Ghost= new int[numChildren]; 

    //active
    CG_Info_Active= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
    CG_numBCMessages_Active= new int[numChildren]; 

    // buffer
    if (MLMD_BCBufferArea){
      CG_Info_Buffer= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
      CG_numBCMessages_Buffer= new int[numChildren]; 
    }

    // fix3B
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      CG_Info_Fix3B= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
      CG_numBCMessages_Fix3B= new int[numChildren]; 
    }

    for (int i=0; i<numChildren; i++){
      CG_numBCMessages_Ghost[i]=0;
      CG_numBCMessages_Active[i]=0;
      if (MLMD_BCBufferArea){
	CG_numBCMessages_Buffer[i]=0;
      }
      if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
	CG_numBCMessages_Fix3B[i]=0;
      }
    }
    
    initWeightBC_Phase2c(grid, vct, CG_Info_Ghost, CG_numBCMessages_Ghost, -1); 
    initWeightBC_Phase2c(grid, vct, CG_Info_Active, CG_numBCMessages_Active, 0);
    
    if (MLMD_BCBufferArea){
      initWeightBC_Phase2c(grid, vct, CG_Info_Buffer, CG_numBCMessages_Buffer, 7);
    }

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      initWeightBC_Phase2c(grid, vct, CG_Info_Fix3B, CG_numBCMessages_Fix3B, 9);
    }

    // same size for all children, used to build BC msg
    CGMsg = new double [CG_MaxSizeMsg * NumF];

    if (MLMD_BCBufferArea){
      CGMsgBuffer = new double [CG_MaxSizeBufferMsg * NumF];
    }

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      CGMsgFix3B = new double [CG_MaxSizeFix3BMsg * Numfix3B];
    }

    if (VerboseCheck== true and false){ 
      for (int ch=0; ch< numChildren; ch++){ // cycle on the children
	for (int m=0; m< CG_numBCMessages_Ghost[ch]; m++){ // cycle on the messages per child
	  // my rank as a parent on this particular parent- child communicator
	  MPI_Comm_rank(vct->getCommToChild(ch) , &rank_As_Parent);
	  cout << "GHOST: I am core " << rank_As_Parent << " (on the PC comm) " << " of grid " << numGrid <<  " and I will communicate with RG core " << CG_Info_Ghost[ch][m].RG_core << " of my " <<ch << "-th child, grid " << vct->getChildGridNum(ch) << " for ghost info" <<endl;  
	}
      }
      // active
      for (int ch=0; ch< numChildren; ch++){ // cycle on the children
	for (int m=0; m< CG_numBCMessages_Active[ch]; m++){ // cycle on the messages per child
	  // my rank as a parent on this particular parent- child communicator
	  MPI_Comm_rank(vct->getCommToChild(ch) , &rank_As_Parent);
	  cout << "ACTIVE: I am core " << rank_As_Parent << " (on the PC comm) " << " of grid " << numGrid <<  " and I will communicate with RG core " << CG_Info_Active[ch][m].RG_core << " of my " <<ch << "-th child, grid " << vct->getChildGridNum(ch) << " for active info" <<endl;
	}
      }
      // buffer
      for (int ch=0; ch< numChildren; ch++){ // cycle on the children
	for (int m=0; m< CG_numBCMessages_Buffer[ch]; m++){ // cycle on the messages per child
	  // my rank as a parent on this particular parent- child communicator
	  MPI_Comm_rank(vct->getCommToChild(ch) , &rank_As_Parent);
	  cout << "BUFFER: I am core " << rank_As_Parent << " (on the PC comm) " << " of grid " << numGrid <<  " and I will communicate with RG core " << CG_Info_Buffer[ch][m].RG_core << " of my " <<ch << "-th child, grid " << vct->getChildGridNum(ch) << " for buffer info" <<endl;
	}
      }

      
    }
    /* end checks */

  } // end  if (numChildren > 0 ), only for parents
  
  /*MPI_Barrier(MPI_COMM_WORLD);
    if (rank_local==0){
    cout << "I am grid " << numGrid <<", I am after the barrier marking Phase2c" << endl;
    // cout << "exiting now..."<< endl;
    }*/
   

  // if this core is involved in BC, instantiated the vectors to store En+1 (NOT En+theta) before smoothing
  
  if (numChildren >0){
    for (int ch=0; ch<numChildren; ch++){
      if (CG_numBCMessages_Ghost[ch]>0) {BSNeeded= true; break;}
      if (CG_numBCMessages_Active[ch]>0) {BSNeeded= true; break;}
    }
  }

  if (BSNeeded){
    Ex_BS = newArr3(double, nxn, nyn, nzn);
    Ey_BS = newArr3(double, nxn, nyn, nzn);
    Ez_BS = newArr3(double, nxn, nyn, nzn);
  }

  Trim_RGBC_Vectors(vct);
  //checks after trim
  if (CommToParent!= MPI_COMM_NULL and false){
    cout << "RG_numBCMessages_Active: " << RG_numBCMessages_Active << endl;
    /*for (int m=0; m< RG_numBCMessages_Active ; m++){                                                      
      cout <<"AFTER TRIM, grid " <<numGrid << ", m: " << m <<", RGBC_Info_Active[m].RG_core: " << RGBC_Info_Active[m].RG_core << ", RGBC_Info_Active[m].CG_core: " << RGBC_Info_Active[m].CG_core<< endl;
      }*/

    cout << "RG_numBCMessages_Ghost: " << RG_numBCMessages_Ghost << endl;
    /*for (int m=0; m< RG_numBCMessages_Ghost ; m++){                                                      
      cout <<"AFTER TRIM, grid " <<numGrid << ", m: " << m <<", RGBC_Info_Ghost[m].RG_core: " << RGBC_Info_Ghost[m].RG_core << ", RGBC_Info_Ghost[m].CG_core: " << RGBC_Info_Ghost[m].CG_core<< endl;
      }*/

    cout << "RG_numBCMessages_Buffer: " << RG_numBCMessages_Buffer << endl;
    /*for (int m=0; m< RG_numBCMessages_Buffer ; m++){                                                      
      cout <<"AFTER TRIM, grid " <<numGrid << ", m: " << m <<", RGBC_Info_Buffer[m].RG_core: " << RGBC_Info_Buffer[m].RG_core << ", RGBC_Info_Buffer[m].CG_core: " << RGBC_Info_Buffer[m].CG_core<< endl;
      }*/
    cout << "RG_numBCMessages_fix3B: " << RG_numBCMessages_fix3B << endl;
    for (int m=0; m< RG_numBCMessages_fix3B ; m++){                                                      
      cout <<"AFTER TRIM, grid " <<numGrid << ", m: " << m <<", RGBC_Info_fix3B[m].RG_core: " << RGBC_Info_fix3B[m].RG_core << ", RGBC_Info_fix3B[m].CG_core: " << RGBC_Info_fix3B[m].CG_core<< endl;
      }
    
  }
  if (numChildren > 0 and false){
    for (int ch=0; ch< numChildren; ch++){
       cout << "CG_numBCMessages_Active[ch]: " << CG_numBCMessages_Active[ch] <<endl;
       /*for (int d=0; d< CG_numBCMessages_Active[ch]; d++){
	cout << "AFTER TRIM ACTIVE, grid " << numGrid << ", child " << ch << ", CG core " << CG_Info_Active[ch][d].CG_core << " to RG core " << CG_Info_Active[ch][d].RG_core << endl;
	}
      cout << "CG_numBCMessages_Ghost[ch]: " << CG_numBCMessages_Ghost[ch] <<endl;
      //for (int d=0; d< CG_numBCMessages_Ghost[ch]; d++){
	cout << "AFTER TRIM GHOST, grid " << numGrid << ", child " << ch << ", CG core " << CG_Info_Ghost[ch][d].CG_core << " to RG core " << CG_Info_Ghost[ch][d].RG_core << endl;
	}*/
      if (MLMD_BCBufferArea ){
	cout << "CG_numBCMessages_Buffer[ch]: " << CG_numBCMessages_Buffer[ch] <<endl;
	for (int d=0; d< CG_numBCMessages_Buffer[ch]; d++){
	  cout << "AFTER TRIM BUFFER, grid " << numGrid << ", child " << ch << ", CG core " << CG_Info_Buffer[ch][d].CG_core << " to RG core " << CG_Info_Buffer[ch][d].RG_core << endl;
	}
      } // end  if (MLMD_BCBufferArea){

      if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
	cout << "CG_numBCMessages_Fix3B[ch]: " << CG_numBCMessages_Fix3B[ch] <<endl;
	for (int d=0; d< CG_numBCMessages_Fix3B[ch]; d++){
	  cout << "AFTER TRIM BUFFER, grid " << numGrid << ", child " << ch << ", CG core " << CG_Info_Fix3B[ch][d].CG_core << " to RG core " << CG_Info_Fix3B[ch][d].RG_core << endl;
	}
      } // end  if (MLMD_BCBufferArea){
      
      
    } // end for (int ch=0; ch< numChildren; ch++){
} // end if (numChildren>0)
  
  /* these deletes only for core 0 of RGs */  
  if (CommToParent != MPI_COMM_NULL && rank_local==0){
    delete[] RGBC_Info_Ghost_LevelWide;
    delete[] RGBC_Info_Active_LevelWide;
    if (MLMD_BCBufferArea){
      delete[] RGBC_Info_Buffer_LevelWide;
    }
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      delete[] RGBC_Info_fix3B_LevelWide;
    }
  }


  MPI_Barrier(vct->getComm());
  if (rank_local==0){
    cout << "I am grid " << numGrid <<", I am after the barrier at the end of initWeightBC" << endl;
  }
  
  /*if (MLMD_BCBufferArea){
    for (int m=0; m< RG_numBCMessages_Buffer; m++){
      cout << "R" << vct->getCartesian_rank() << " msg " << m << " of " << RG_numBCMessages_Buffer << endl;
      for (int j=0; j< RGBC_Info_Buffer[m].np_x* RGBC_Info_Buffer[m].np_y*RGBC_Info_Buffer[m].np_z; j++ ){
	cout << Ex_Buffer_BC[m][j] <<" " << Ey_Buffer_BC[m][j] <<" " << Ez_Buffer_BC[m][j]<<endl;
	cout << Exth_Buffer_BC[m][j] <<" " << Eyth_Buffer_BC[m][j] <<" " << Ezth_Buffer_BC[m][j]<<endl;
	cout << Bxn_Buffer_BC[m][j] <<" " << Byn_Buffer_BC[m][j] <<" " << Bzn_Buffer_BC[m][j]<<endl; 
      }
    }
    }*/


}


void EMfields3D::initWeightBC_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *RG_numBCMessages, int which){

  // this the size in the direction not used --- put 1 so i can cycle using <
  int DNS=1;

  bool DIR_0= true; // L/R
  bool DIR_1= true; // B/F
  bool DIR_2= true; // top/ bottom
  int countMSG=0; 

  if (! (which== -1 || which ==0 )){
    cout << "initWeightBC, phase 1, is receiving wrong inputs. Check the code. Aborting now..." << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  /*if (vct->getCartesian_rank()==0){
    cout << "Grid number " << numGrid << "entered initWeightBC_Phase1, which="  << which << endl;
    }*/
  int SW_rank=vct->getSystemWide_rank();
  
  int MS= nxn; if (nyn>MS) MS= nyn; if (nzn>MS) MS= nzn; 

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
  // DIR2: core rank in the CommToParent communicator
  int *Dir2_rank= new int [MS];  // this does not need to be this big, but anyway
  int *Dir2_IndexFirstPointperC= new int [MS];
  int Dir2_Ncores=0;
  /*******************************************************************/

  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  int rank_CTP= vct->getRank_CommToParent();
  int rank_G= vct->getSystemWide_rank();
  int rank_local= vct->getCartesian_rank();

  // NB: _s and _e are included!!!, so <= in the for
  int i_s, i_e;
  int j_s, j_e;
  int k_s, k_e;

  string FACE;

  // policy:
  // 1. explore the higher direction (e.g., y as opposed to x in bottom)
  // 2. do a for on the number of cores found there, and inside the for explore the lower direction (e.g., x in bottom)
  // 3. commit the message inside the for
  
    
  // this is the bottom face
  if (vct->getCoordinates(2)==0 && vct->getZleft_neighbor()== MPI_PROC_NULL && DIR_2){

    i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    // cycle on x and y
    if (which ==0){
      k_s=1; k_e=1; 
    }
    else if (which==-1){
      // i have to pass really everything because i won't be able to do a communicate
      k_s=0; k_e=0;
    }
    
    countMSG =  (*RG_numBCMessages);
    FACE= "BOTTOM";
    // Dir1, higher dimension: Y
    grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, k_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);
    
    for (int n=0; n<Dir1_Ncores; n++){ // it will find again the core in Dir 1, but it will also explore Dir 2
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e,  Dir1_IndexFirstPointperC[n], k_s, Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); // Dir2, lower dimension: X
      
      // build and commit each of these
      for (int NN=0; NN<Dir2_Ncores; NN++){
	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), Dir2_IndexFirstPointperC[NN], Dir1_IndexFirstPointperC[n], k_s, -1, Dir2_NPperC[NN], Dir1_NPperC[n], DNS, Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);
	
	//cout << "R" <<rank_G << " NN " << NN  << "BOTTOM: which " << which << ", init ["  <<  Dir2_IndexFirstPointperC[NN] <<"-" << Dir1_IndexFirstPointperC[n] <<"-" << k_s << " # points [" <<Dir2_NPperC[NN] <<"-" <<Dir1_NPperC[n] <<"-" <<DNS <<   "]"<<endl;
	
	(*RG_numBCMessages)++;
	
	int tmp= Dir2_NPperC[NN] *Dir1_NPperC[n] ;
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;
      
      
      }
      
    } // end for (int n=0; n<Dir1_Ncores; n++)

    
    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
    
  } // end bottom face
  
  // this is the top face
  if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor() == MPI_PROC_NULL && DIR_2){ 
    
    i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    // cycle on x and y
    if (which ==0){
      k_s=nzn-2; k_e=nzn-2; 
    }
    else if (which==-1){
      k_s=nzn-1; k_e=nzn-1;
    }

 
    countMSG =  (*RG_numBCMessages);
    FACE= "TOP";
    // Dir1, higher dimension: Y
    grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, k_s, Dir1_SPXperC, Dir1_SPYperC,Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);
    
    for (int n=0; n< Dir1_Ncores; n++){
      grid->RGBCExploreDirection(vct,FACE, 0, i_s, i_e, Dir1_IndexFirstPointperC[n], k_s, Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); // Dir2, lower dimension: X

      for (int NN=0; NN<Dir2_Ncores; NN++){
	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), Dir2_IndexFirstPointperC[NN], Dir1_IndexFirstPointperC[n], k_s, -1, Dir2_NPperC[NN], Dir1_NPperC[n], DNS, Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);


	(*RG_numBCMessages)++;

	int tmp= Dir2_NPperC[NN]* Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;
	

      }
    } // endfor (int n=0; n<Dir1_Ncores; n++)   

    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
  } // end top face
  
  // this is the left face
  if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor() == MPI_PROC_NULL && DIR_0){ 


    j_s=0; j_e= nyn-1;
    k_s=0; k_e= nzn-1;
    // cycle on y and z
    if (which ==0){
      i_s=1; i_e= 1;
    }
    else if (which==-1){
      i_s=0; i_e= 0;	       
      /*if (vct->getCoordinates(1)==0) j_s=0; else j_s=1;
	if (vct->getCoordinates(1)==YLEN-1) j_e=nyn-1; else j_e=nyn-2;
	if (vct->getCoordinates(2)==0) k_s=0; else k_s=1;
	if (vct->getCoordinates(2)==ZLEN-1) k_e=nzn-1; else k_e=nzn-2;*/
    }
    /* left face: 
       Dir 1 (higher dir) --> Z
       Dir 2 (lower dir) --> Y */
    
    countMSG =  (*RG_numBCMessages);
    FACE= "LEFT";

    // Dir1, higher dimensionality: Z
    grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

    for (int n=0; n< Dir1_Ncores; n++){
      // Dir2, Y
      grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n], Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); 

      for (int NN=0; NN< Dir2_Ncores; NN++){

	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), i_s, Dir2_IndexFirstPointperC[NN], Dir1_IndexFirstPointperC[n], -1, DNS, Dir2_NPperC[NN], Dir1_NPperC[n], Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);


	(*RG_numBCMessages)++;

	int tmp=  Dir2_NPperC[NN]*Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;
	
      }
    }// end for(int n=0; n< Dir1_Ncores; n++)
    
    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
  } // end left face
  
  // this is the right face
  if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor() == MPI_PROC_NULL && DIR_0){ 
    
    j_s=0; j_e= nyn-1;
    k_s=0; k_e= nzn-1;
    // cycle on y and z
    if (which ==0){
      i_s=nxn-2; i_e= nxn-2;
    }
    else if (which==-1){
      i_s=nxn-1; i_e= nxn-1;	       
      /*if (vct->getCoordinates(1)==0) j_s=0; else j_s=1;
	if (vct->getCoordinates(1)==YLEN-1) j_e=nyn-1; else j_e=nyn-2;
	if (vct->getCoordinates(2)==0) k_s=0; else k_s=1;
	if (vct->getCoordinates(2)==ZLEN-1) k_e=nzn-1; else k_e=nzn-2;*/
    }
    /* right face:
       Dir 1 (higher dir) --> Z
       Dir 2 (lower dir) --> Y */

    countMSG =  (*RG_numBCMessages);
    FACE= "RIGHT";

    // Dir1: higher dimensionality --> Z
    grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);
    

    for (int n=0; n< Dir1_Ncores; n++){
      // Dir2, Y
      grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n], Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC);

      for (int NN=0; NN< Dir2_Ncores; NN++){

	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), i_s, Dir2_IndexFirstPointperC[NN], Dir1_IndexFirstPointperC[n], -1, DNS, Dir2_NPperC[NN], Dir1_NPperC[n], Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);

	(*RG_numBCMessages)++;


	int tmp= Dir2_NPperC[NN]*Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;
	
      }
    }

    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
  } // end right face
  
  // this is the front face
  if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor() == MPI_PROC_NULL && DIR_1){

    i_s=0; i_e= nxn-1;
    k_s=0; k_e= nzn-1;
    // cycle on x and z
    if (which ==0){
      j_s=1; j_e= 1;
    }
    else if (which==-1){
      j_s=0; j_e= 0;	       
      /*if (vct->getCoordinates(0)==0) i_s=0; else i_s=1;
	if (vct->getCoordinates(0)==XLEN-1) i_e=nxn-1; else i_e=nxn-2;
	if (vct->getCoordinates(2)==0) k_s=0; else k_s=1;
	if (vct->getCoordinates(2)==ZLEN-1) k_e=nzn-1; else k_e=nzn-2;*/
    }
    /* front face:
       Dir 1 (higher dir) --> Z
       Dir 2 (lower dir) -->  X */

    countMSG =  (*RG_numBCMessages);
    FACE= "FRONT";
 
    // Dir 1, higher dimensionality: Z   
    grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

    for (int n=0; n< Dir1_Ncores; n++){
      // Dir2, X
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e, j_s, Dir1_IndexFirstPointperC[n], Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC);

      for (int NN=0; NN< Dir2_Ncores; NN++){
	
	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), Dir2_IndexFirstPointperC[NN], j_s, Dir1_IndexFirstPointperC[n], -1, Dir2_NPperC[NN], DNS, Dir1_NPperC[n], Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);

	(*RG_numBCMessages)++;

	int tmp= Dir2_NPperC[NN] * Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;
	
	
      }
    } // end for(int n=0; n< Dir1_Ncores; n++){

    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
  } // end front face
  
  // this is the back face
  if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor() == MPI_PROC_NULL && DIR_1){
    
    i_s=0; i_e= nxn-1;
    k_s=0; k_e= nzn-1;
    // cycle on x and z
    if (which ==0){
      j_s=nyn-2; j_e= nyn-2;
    }
    else if (which==-1){
      j_s=nyn-1; j_e= nyn-1;	       
      /*if (vct->getCoordinates(0)==0) i_s=0; else i_s=1;
	if (vct->getCoordinates(0)==XLEN-1) i_e=nxn-1; else i_e=nxn-2;
	if (vct->getCoordinates(2)==0) k_s=0; else k_s=1;
	if (vct->getCoordinates(2)==ZLEN-1) k_e=nzn-1; else k_e=nzn-2;*/
    }
    /* back face:
       Dir 1 (higher dimensionality) --> Z
       Dir 2 (lower dimensionality) --> X */
    
    countMSG =  (*RG_numBCMessages);    
    FACE= "BACK";

    // Dir 1, higher dimensionality: Z
    grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);
    
    for (int n=0; n< Dir1_Ncores; n++){
      // dir2, X
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e, j_s, Dir1_IndexFirstPointperC[n], Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC);
      
      for (int NN=0; NN< Dir2_Ncores; NN++){
	Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), Dir2_IndexFirstPointperC[NN], j_s, Dir1_IndexFirstPointperC[n], -1, Dir2_NPperC[NN], DNS, Dir1_NPperC[n], Dir2_SPXperC[NN], Dir2_SPYperC[NN], Dir2_SPZperC[NN], Dir2_rank[NN], rank_CTP, *RG_numBCMessages);

	(*RG_numBCMessages)++;

	int tmp= Dir2_NPperC[NN]* Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize) RG_MaxMsgSize= tmp;

      }
    } // end for(int n=0; n< Dir1_Ncores; n++){*/

    /*cout << "FACE " << FACE << ", Inside Phase 1, which " << which << endl;
      cout << "R" << SW_rank <<" send BC messages to " << (*RG_numBCMessages)-countMSG << "core(s)"  << endl;
      for (int i= countMSG; i< *RG_numBCMessages; i++){ 
      cout << " parent core " << RGBC_Info[i].CG_core <<" " << endl;
      }*/
  } // end back face
  

  // mega test: interpolating the entire grid (works for one core)
  /*
  if (which==-1){                                               
    cout << "MEGA TEST ENABLED-- DISABLE IN THE CODE" << endl;
   *RG_numBCMessages=0;
   Assign_RGBC_struct_Values(RGBC_Info + (*RG_numBCMessages), 1, 1, 1, -1, nxn-2, nyn-2, nzn-2, grid->getXN_P(1,1,1), grid->getYN_P(1,1,1), grid->getZN_P(1,1,1), 0, rank_CTP, *RG_numBCMessages);  
   (*RG_numBCMessages)++;                 
   RG_MaxMsgSize= nxn*nyn*nzn;                                                    
   }
  */
  /* */
  
  /* this only for testing communication; comment in the end */
  //TEST__Assign_RG_BC_Values(vct, RGBC_Info, RG_numBCMessages, which); // input with same # of cores per grid
  //TEST__Assign_RG_BC_Values_DNC(vct, RGBC_Info, RG_numBCMessages, which); // input with different # of cores per grid  
  /* END THIS IS A TEST; COMMENT IN THE END */
  
  /* for further use, i need to set the RG_core field of the first unused slot to -1 
     but DO NOT MODIFY THE NUMBER OF MSGs;
     I will just send a +1 */

  //cout << "R" <<SW_rank <<"RG_numBCMessages= " <<*RG_numBCMessages <<endl;
  RGBC_Info[*RG_numBCMessages].RG_core= -1;
  RGBC_Info[*RG_numBCMessages].CG_core= -1;
  //cout << "R" <<  rank_G <<":  END SONO QUI!!!" << endl;


  // all the deletes

  delete []Dir1_SPXperC;
  delete []Dir1_SPYperC;
  delete []Dir1_SPZperC;
  delete []Dir1_NPperC;
  delete []Dir1_rank;
  delete []Dir1_IndexFirstPointperC;

  delete []Dir2_SPXperC;
  delete []Dir2_SPYperC;
  delete []Dir2_SPZperC;
  delete []Dir2_NPperC;
  delete []Dir2_rank;
  delete []Dir2_IndexFirstPointperC;
}

void EMfields3D::initWeightBCBuffer_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg){

  // this the size in the direction not used --- put 1 so i can cycle using <
  int DNS=1;

  bool DIR_0= true; // L/R
  bool DIR_1= true; // B/F
  bool DIR_2= true; // top/ bottom

  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  int rank_CTP= vct->getRank_CommToParent();
  int rank_G= vct->getSystemWide_rank();
  int rank_local= vct->getCartesian_rank();

  // NB: _s and _e are included!!!, so <= in the for
  int i_s, i_e;
  int j_s, j_e;
  int k_s, k_e;

  string FACE;

  // policy:
  // 1. explore the higher direction (e.g., y as opposed to x in bottom)
  // 2. do a for on the number of cores found there, and inside the for explore the lower direction (e.g., x in bottom)
  // 3. commit the message inside the for

  char DIR;
    
  // this is the bottom face
  if (vct->getCoordinates(2)==0 && vct->getZleft_neighbor()== MPI_PROC_NULL && DIR_2){

    i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    k_s=2; k_e=k_s+ BufZ; 
    
    if(k_e > nzn-1) {
      cout << "WARNING: the BufZ value you selected is too high w.r.t. the number of cells; I am reducing it so that k_s+BufZ= nzn-1 in initWeightBCBuffer_Phase1" << endl;
      BufZ= nzn-1-k_s;
      k_e= nzn-1;
      cout << "BufZ is now " << BufZ << endl;
    }


    DIR= 'B';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);

  } // end bottom face
  


  // this is the top face
  if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor() == MPI_PROC_NULL && DIR_2){ 
    
    i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    k_e= nzn-3; k_s= k_e- BufZ;
    if(k_s <0) 
      {
	cout << "WARNING: the BufZ value you selected is too high w.r.t. the number of cells; I am reducing it so that k_s=k_e-BufZ= 0 in initWeightBCBuffer_Phase1" << endl;
	BufZ= k_e;
	k_s= 0;
	cout << "BufZ is now " <<BufZ <<endl;
      }
    DIR= 'T';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
  } // end top face
  
  // this is the left face
  if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor() == MPI_PROC_NULL && DIR_0){ 

    j_s=0; j_e= nyn-1;
    k_s=0; k_e= nzn-1;
    i_s=2; i_e= i_s+ BufX;
    if(i_e > nxn-1) {
      cout << "WARNING: the BufX value you selected is too high w.r.t. the number of cells; I am reducing it so that i_s+BufX= nxn-1 in initWeightBCBuffer_Phase1" << endl;
      BufX= nxn-1-i_s;
      i_e= nxn-1;
      cout << "BufX is now " <<BufX <<endl;
    }


    DIR= 'L';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);

  } // end left face
  
  // this is the right face
  if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor() == MPI_PROC_NULL && DIR_0){ 
    
    j_s=0; j_e= nyn-1;
    k_s=0; k_e= nzn-1;
    i_e=nxn-3; i_s= i_e- BufX;
    if(i_s <0)
      {
	cout << "WARNING: the BufX value you selected is too high w.r.t. the number of cells; I am reducing it so that i_e-BufX= 0 in initWeightBCBuffer_Phase1" << endl;
	BufX= i_e;
	i_s= 0;
	cout << "BufX is now " <<BufX <<endl;
      }

    DIR= 'R';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);

  } // end right face
  
  // this is the front face
  if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor() == MPI_PROC_NULL && DIR_1){

    i_s=0; i_e= nxn-1;
    k_s=0; k_e= nzn-1;
    j_s=2; j_e= j_s+ BufY;
    if(j_e > nyn-1) {
      cout << "WARNING: the BufY value you selected is too high w.r.t. the number of cells; I am reducing it so that j_s+BufY= nyn-1 in initWeightBCBuffer_Phase1" << endl;
      BufY= nyn-1-j_s;
      j_e= nyn-1;
      cout << "BufY is now " <<BufY <<endl;
    }
    DIR= 'F';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);

  } // end front face
  
  // this is the back face
  if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor() == MPI_PROC_NULL && DIR_1){
    
    i_s=0; i_e= nxn-1;
    k_s=0; k_e= nzn-1;
    j_e=nyn-3; j_s= j_e- BufY;
    if(j_s <0) 
    {
      cout << "WARNING: the BufY value you selected is too high w.r.t. the number of cells; I am reducing it so that j_e-BufY= 0 in initWeightBCBuffer_Phase1" << endl;
      BufY= j_e;
      j_s= 0;
      cout << "BufY is now " <<BufY <<endl;
    }

    DIR= 'b';
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);

  } // end back face
  

  RGBC_Info[*numMsg].RG_core= -1;
  RGBC_Info[*numMsg].CG_core= -1;

}

/** fix3B **/
/** NB: I am reusing the same structure for InterpolateOldBCell; 
    anyway, if InterpolateOldBCell, fix3B is included **/
void EMfields3D::initWeightBCfix3B_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg){

  // this the size in the direction not used --- put 1 so i can cycle using <
  int DNS=1;

  bool DIR_0= true; // L/R
  bool DIR_1= true; // B/F
  bool DIR_2= true; // top/ bottom

  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  int rank_CTP= vct->getRank_CommToParent();
  int rank_G= vct->getSystemWide_rank();
  int rank_local= vct->getCartesian_rank();

  // NB: _s and _e are included!!!, so <= in the for
  int i_s, i_e;
  int j_s, j_e;
  int k_s, k_e;

  string FACE;

  // policy:
  // 1. explore the higher direction (e.g., y as opposed to x in bottom)
  // 2. do a for on the number of cores found there, and inside the for explore the lower direction (e.g., x in bottom)
  // 3. commit the message inside the for

  char DIR;
    
  if (MLMD_fixBCenters){

    // this is the bottom face
    if (vct->getCoordinates(2)==0 && vct->getZleft_neighbor()== MPI_PROC_NULL && DIR_2){
      
      i_s=0; i_e= nxc-1;
      j_s=0; j_e= nyc-1;
      k_s=0; k_e=2;  
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) k_e= k_e+BufZ-1;
      
      DIR= 'B';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end bottom face
    
    
    // this is the top face
    if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor() == MPI_PROC_NULL && DIR_2){ 
      
      i_s=0; i_e= nxc-1;
      j_s=0; j_e= nyc-1;
      k_e= nzc-1; k_s= nzc-3;
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) k_s= k_s-(BufZ-1);
      
      DIR= 'T';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end top face
    
    // this is the left face
    if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor() == MPI_PROC_NULL && DIR_0){ 
      
      j_s=0; j_e= nyc-1;
      k_s=0; k_e= nzc-1;
      i_s=0; i_e= 2;
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) i_e= i_e+BufX-1;
      
      DIR= 'L';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end left face
    
    // this is the right face
    if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor() == MPI_PROC_NULL && DIR_0){ 
      
      j_s=0; j_e= nyc-1;
      k_s=0; k_e= nzc-1;
      i_e=nxc-1; i_s= nxc-3;
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) i_s= i_s-(BufX-1);
      
      DIR= 'R';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end right face
    
    // this is the front face
    if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor() == MPI_PROC_NULL && DIR_1){
      
      i_s=0; i_e= nxc-1;
      k_s=0; k_e= nzc-1;
      j_s=0; j_e= 2;
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) j_e= j_e+BufY-1;
      
      DIR= 'F';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end front face
    
    // this is the back face
    if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor() == MPI_PROC_NULL && DIR_1){
      
      i_s=0; i_e= nxc-1;
      k_s=0; k_e= nzc-1;
      j_e=nyc-1; j_s= nyc-3;
      
      // the particles, repopulated before moving them, are moved in fields completely interpolated; here, (+ Buf-1)
      // considering that Buf nodes are ON TOP OF ghost + active, and all cells are included
      if (RepopulateBeforeMover) j_s= j_s-(BufY-1);
      
      DIR= 'b';
      grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
      
    } // end back face
  } // end if (MLMD_fixBCenters)

  if (MLMD_InterpolateOldBCell){
    
    if (vct->getCartesian_rank() == 0) {
      cout << "ATTENTION: I will interpolate values of ALL B^n cells" << endl;
    }
    
    /* interpolating EVERYTHING */
    i_s=0; i_e= nxc-1;
    j_s=0; j_e= nyc-1;
    k_s=0; k_e= nzc-1;  
    
    DIR= 'B';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
  }

  RGBC_Info[*numMsg].RG_core= -1;
  RGBC_Info[*numMsg].CG_core= -1;

}

/** end fix3B **/


/* phase 2a of initWeightBC:
   core 0 of each child grid receives all the messages to send to the corresponding coarse grid 
   a level-wide message structure is built */
void EMfields3D::initWeightBC_Phase2a(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info_LevelWide,int * RG_numBCMessages_LevelWide, RGBC_struct *RGBC_Info, int RG_numBCMessages, int which){

  // this sizes are local to the grid
  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  MPI_Status status;

  int TAG;
  if (which==-1) TAG=TAG_BC_GHOST;
  if (which==0) TAG=TAG_BC_ACTIVE;
  if (which==3) TAG=TAG_PROJ;
  if (which==-10) TAG=TAG_II;
  if (which == 7) TAG= TAG_BC_BUFFER;
  if (which == 9) TAG= TAG_BC_FIX3B;

  /* here, core 0 local receives one message per other RG core */
  RGBC_struct * buffer_rcv_Core0;
  buffer_rcv_Core0 =  new RGBC_struct[MAX_RG_numBCMessages];
  
  //cout << "MAX number of msg for level 0: " << MAX_size_LevelWide << endl;
  
  //*RG_numBCMessages_LevelWide= 0;
  
  // first, core0 copies its own messages into the level-wide structure
  
  for (int m=0; m<RG_numBCMessages; m++){
    RGBC_Info_LevelWide[m]= RGBC_Info[m];
  }
  
  *RG_numBCMessages_LevelWide= RG_numBCMessages;
  

  string WW;
  if (which==0) WW= "ACTIVE";
  if (which==1) WW= "GHOST";
  if (which==3) WW= "PROJ";
  if (which==-10) WW= "II";
  if (which== 7) WW= "BUFFER";
  if (which== 9) WW= "FIX3B";
  /*cout <<WW << " grid " << numGrid << " after core 0 has copies his messages, RG_numBCMessages_LevelWide= " << * RG_numBCMessages_LevelWide <<endl;*/
  // end copying, starts receiving
  
  
  for (int c=0; c<XLEN*YLEN*ZLEN-1; c++){
    int rcvd=0;
    MPI_Recv(buffer_rcv_Core0, MAX_RG_numBCMessages, MPI_RGBC_struct, MPI_ANY_SOURCE, TAG, vct->getComm(), &status);
    while( buffer_rcv_Core0[rcvd].RG_core != -1  ){
      
      
      RGBC_Info_LevelWide[*RG_numBCMessages_LevelWide]= buffer_rcv_Core0[rcvd]; 
      (*RG_numBCMessages_LevelWide)++;
      
      rcvd++;
      if (*RG_numBCMessages_LevelWide==MAX_size_LevelWide){
	cout << "initWeightBC_Phase2a: The number of msgs that you plan on sending (>" << *RG_numBCMessages_LevelWide <<") exceeds capabilities;" << endl;
	cout << "increase buffer size; aborting now " << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
      
    }
  } // end for (int c=0; c<XLEN*YLEN*ZLEN-1; c++){
  
  delete[] buffer_rcv_Core0;
}

/*core 0 of the child grid assembles & sends messages for all CG cores*/
void EMfields3D::initWeightBC_Phase2b(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info_LevelWide,int RG_numBCMessages_LevelWide, int which){

  /* here, vectors where core 0 of the local child grid assembles the messages                               
     to all the coarse grid core in the parent-child communicator         
     de-allocated at the end of the function */
  /* one message per CG core -- hence, the sizes are the sizes of the CG */

  RGBC_struct ** RGBC_Info_ToCGCore;
  int * RG_numBCMessages_ToCGCore;

  // # of cores in the parent, differnt from the # of cores in the child
  int parentGrid= vct->getParentGridNum();
  int ParentCoreNum= vct->getXLEN(parentGrid) * vct->getYLEN(parentGrid) *vct->getZLEN(parentGrid);

  int TAG;
  if (which==-1) TAG=TAG_BC_GHOST;
  if (which==0) TAG=TAG_BC_ACTIVE;
  if (which==3) TAG=TAG_PROJ;
  if (which==-10) TAG=TAG_II;
  if (which==7) TAG=TAG_BC_BUFFER;
  if (which == 9) TAG= TAG_BC_FIX3B;

  MPI_Comm CommToParent= vct->getCommToParent();

  RGBC_Info_ToCGCore =  newArr2(RGBC_struct, ParentCoreNum, MAX_RG_numBCMessages);
  RG_numBCMessages_ToCGCore= new int[ParentCoreNum];
  
  for (int c=0; c< ParentCoreNum; c++) {
    RG_numBCMessages_ToCGCore[c]=0;
  } // end for (int c=0; c< XLEN*YLEN*ZLEN; c++)
  
  for (int m=0; m< RG_numBCMessages_LevelWide; m++) {
    int where= RGBC_Info_LevelWide[m].CG_core;
    RGBC_Info_ToCGCore[where][RG_numBCMessages_ToCGCore[where]]=RGBC_Info_LevelWide[m] ;
    RG_numBCMessages_ToCGCore[where]++;

    if (RG_numBCMessages_ToCGCore[where] == MAX_RG_numBCMessages){
      cout << "initWeightBC_Phase2b: The number of msgs that you plan on sending (>" << MAX_RG_numBCMessages <<") exceeds capabilities;" << endl;
      cout << "increase buffer size; aborting now " << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  
  // this to stop the read when receiving; need to send one more
  for (int cg=0; cg< ParentCoreNum; cg++){
    RGBC_Info_ToCGCore[cg][RG_numBCMessages_ToCGCore[cg]].RG_core= -1;
  }
    
  // core 0 of local child grid send a message to all CG cores on the parent-child communicator
  // need to send a message more because the RG_core field of the extra msg contains the -1 flag

  // NB: this is the size of the parent, not of the child
  for (int cg=0; cg <ParentCoreNum; cg++){
    //cout << "grid " << numGrid << " sends first, RG_Core " << RGBC_Info_ToCGCore[cg][0].RG_core << endl;
    //cout << "ACTIVE: grid " <<numGrid << ", msg to CG core " << cg << ": "<< RG_numBCMessages_ToCGCore[cg] << "+1" << endl; // ok up to here

    MPI_Send(&(RGBC_Info_ToCGCore[cg][0]),RG_numBCMessages_ToCGCore[cg] +1, MPI_RGBC_struct, cg, TAG, CommToParent  );
    
  }

  delArr2(RGBC_Info_ToCGCore, ParentCoreNum);
  delete[] RG_numBCMessages_ToCGCore;
}

/* phase 2c, only for the parent                                                                          
   each CG core will receive a message per child grid  */
void EMfields3D::initWeightBC_Phase2c(Grid *grid, VirtualTopology3D *vct, RGBC_struct ** CG_Info, int* CG_numBCMessages, int which){

  MPI_Status status;

  int tmp;
  int tmp_proj;
  int tmp_PPP;
  /* here, each CG core receive messages for core 0 (local) form RG */
  RGBC_struct * CG_buffer;
  CG_buffer= new RGBC_struct[MAX_RG_numBCMessages];

  int TAG;
  if (which==-1) TAG=TAG_BC_GHOST;
  if (which==0) TAG=TAG_BC_ACTIVE;
  if (which==3) TAG= TAG_PROJ;
  if (which==-10) TAG= TAG_II;
  if (which==7) TAG= TAG_BC_BUFFER;
  if (which == 9) TAG= TAG_BC_FIX3B;

  int rank_As_Parent;
  
  for (int ch=0; ch < numChildren; ch ++){
    // receive 

    MPI_Comm_rank(vct->getCommToChild(ch), &rank_As_Parent );

    MPI_Recv(CG_buffer, MAX_RG_numBCMessages, MPI_RGBC_struct, MPI_ANY_SOURCE, TAG, vct->getCommToChild(ch), &status);
    
    // process
    while (CG_buffer[CG_numBCMessages[ch]].RG_core != -1 ){
      CG_Info[ch][CG_numBCMessages[ch]]= CG_buffer[CG_numBCMessages[ch]];
      
      // the size of this msg; consider extrems are INCLUDED
      // the max value will be used as max size for BC msg
      
      if (which == -1 or which==0){ 
	tmp=( CG_buffer[CG_numBCMessages[ch]].np_x ) * ( CG_buffer[CG_numBCMessages[ch]].np_y ) * ( CG_buffer[CG_numBCMessages[ch]].np_z);
	if ( tmp> CG_MaxSizeMsg  ) CG_MaxSizeMsg = tmp;
      }

      else if (which ==3){
	tmp_proj=( CG_buffer[CG_numBCMessages[ch]].np_x ) * ( CG_buffer[CG_numBCMessages[ch]].np_y ) * ( CG_buffer[CG_numBCMessages[ch]].np_z);
	if ( tmp_proj> size_CG_ProjMsg  ) size_CG_ProjMsg = tmp_proj;
      }

      else if (which ==-10){
	tmp_proj=( CG_buffer[CG_numBCMessages[ch]].np_x ) * ( CG_buffer[CG_numBCMessages[ch]].np_y ) * ( CG_buffer[CG_numBCMessages[ch]].np_z);
	if ( tmp_proj> CG_MaxSizeMsg_II) CG_MaxSizeMsg_II = tmp_proj;
      }
      
      else if (which ==7){
	tmp_proj=( CG_buffer[CG_numBCMessages[ch]].np_x ) * ( CG_buffer[CG_numBCMessages[ch]].np_y ) * ( CG_buffer[CG_numBCMessages[ch]].np_z);
	if ( tmp_proj> CG_MaxSizeBufferMsg) CG_MaxSizeBufferMsg = tmp_proj;
      }
      else if (which ==9){
	tmp_proj=( CG_buffer[CG_numBCMessages[ch]].np_x ) * ( CG_buffer[CG_numBCMessages[ch]].np_y ) * ( CG_buffer[CG_numBCMessages[ch]].np_z);
	if ( tmp_proj> CG_MaxSizeFix3BMsg) CG_MaxSizeFix3BMsg = tmp_proj;
      }
      else{
	cout << "Wrong tag, aborting..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
	
      /* some sanity check ---
	 the CG_core field of the msg should correspond to the rank of this core 
	 on the parent-child communicator: check */
      
      if (rank_As_Parent != (CG_Info[ch][CG_numBCMessages[ch]]).CG_core ){
	cout << "initWeightBC_Phase2c has detected an anomaly in msg from RG grid ..." << endl;
	cout << "aborting now ..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      } //else {cout << "Sanity check in initWeightBC_Phase2c passed ..." << endl;}
      /* end sanity check */
      
      CG_numBCMessages[ch]++;
    } // end while

    if (which==3){ // to instantiate the vector with the weight
      tmp_PPP= CG_numBCMessages[ch];
      if (tmp_PPP> Max_CG_numProjMessages) {Max_CG_numProjMessages= tmp_PPP;}
    }

  } // end for (int ch=0; ch < numChildren; ch ++){

  delete[]CG_buffer;
}

/* Trim functions: RGBC_Info_** and CG_Info_** are initially allocated very large for                     
   internal reasons;                                                                                      
   once initWeightBC has finished, we want to free to unused memory */
/* this for RGBC_Info_*, vector with size MAX_RG_numBCMessages                                            
   as desired_Size, pass the # of slots used  */
void EMfields3D::Trim_RGBC_Vectors(VirtualTopology3D *vct){

  int dim;
  MPI_Comm CommToParent= vct->getCommToParent();
  
  if (CommToParent != MPI_COMM_NULL){
    // ghost
    dim= RG_numBCMessages_Ghost+1;
    RGBC_struct *tmp= new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      tmp[i]= RGBC_Info_Ghost[i];
    }
    delete[] RGBC_Info_Ghost;
    RGBC_Info_Ghost = new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      RGBC_Info_Ghost[i]= tmp[i];
    }
    delete[]tmp;

    // active
    dim= RG_numBCMessages_Active+1;
    tmp= new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      tmp[i]= RGBC_Info_Active[i];
    }
    delete[] RGBC_Info_Active;
    RGBC_Info_Active = new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      RGBC_Info_Active[i]= tmp[i];
    }
    delete[]tmp;

    // buffer
    if (MLMD_BCBufferArea){
      dim= RG_numBCMessages_Buffer+1;
      tmp= new RGBC_struct[dim];
      
      for (int i=0; i<dim; i++){
	tmp[i]= RGBC_Info_Buffer[i];
      }
      delete[] RGBC_Info_Buffer;
      RGBC_Info_Buffer = new RGBC_struct[dim];
      
      for (int i=0; i<dim; i++){
	RGBC_Info_Buffer[i]= tmp[i];
      }
      delete[]tmp;

      char* TMP;
      TMP= new char[dim];
      
      for (int i=0; i<dim; i++){
	TMP[i] = DirBuffer[i];
      }
      delete[] DirBuffer;
      DirBuffer= new char[dim];

      for (int i=0; i<dim; i++){
	DirBuffer[i]= TMP[i];
      }
      delete[] TMP;
      
    } // end if (MLMD_BCBufferArea){

    /** fix3B **/
    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      dim= RG_numBCMessages_fix3B+1;
      tmp= new RGBC_struct[dim];
      
      for (int i=0; i<dim; i++){
	tmp[i]= RGBC_Info_fix3B[i];
      }
      delete[] RGBC_Info_fix3B;
      RGBC_Info_fix3B = new RGBC_struct[dim];
      
      for (int i=0; i<dim; i++){
	RGBC_Info_fix3B[i]= tmp[i];
      }
      delete[]tmp;

    } // end if (MLMD_fixBCenters){
    /** end fix3B **/
  }
  
  if (numChildren >0){
    // ghost
    dim=0;
    for (int ch=0; ch< numChildren; ch++){
      if (CG_numBCMessages_Ghost[ch]>dim) {dim= CG_numBCMessages_Ghost[ch];}
    }
    dim= dim+1;
    RGBC_struct **tmp2= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	tmp2[ch][d]= CG_Info_Ghost[ch][d];
      }
    }
    delArr2(CG_Info_Ghost, numChildren);
    CG_Info_Ghost= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	CG_Info_Ghost[ch][d]= tmp2[ch][d];
      }
    }
   
    delArr2(tmp2, numChildren);

    // active
    dim=0;
    for (int ch=0; ch< numChildren; ch++){
      if (CG_numBCMessages_Active[ch]>dim) {dim= CG_numBCMessages_Active[ch];}
    }
    dim=dim+1;
    tmp2= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	tmp2[ch][d]= CG_Info_Active[ch][d];
      }
    }
    delArr2(CG_Info_Active, numChildren);
    CG_Info_Active= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	CG_Info_Active[ch][d]= tmp2[ch][d]; //SF not here
      }
    }
    delArr2(tmp2, numChildren);


    if (MLMD_BCBufferArea){
      dim=0;
      for (int ch=0; ch< numChildren; ch++){
	if (CG_numBCMessages_Buffer[ch]>dim) {dim= CG_numBCMessages_Buffer[ch];}
      }
      dim=dim+1;
      tmp2= newArr2(RGBC_struct, numChildren, dim);
      for (int ch=0; ch<numChildren; ch++){
	for (int d=0; d<dim; d++){
	  tmp2[ch][d]= CG_Info_Buffer[ch][d];
	}
      }
      delArr2(CG_Info_Buffer, numChildren);
      CG_Info_Buffer= newArr2(RGBC_struct, numChildren, dim);
      for (int ch=0; ch<numChildren; ch++){
	for (int d=0; d<dim; d++){
	  CG_Info_Buffer[ch][d]= tmp2[ch][d]; //SF not here
	}
      }
      delArr2(tmp2, numChildren);
     
    } // end if (MLMD_BCBufferArea){

    if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
      dim=0;
      for (int ch=0; ch< numChildren; ch++){
	if (CG_numBCMessages_Fix3B[ch]>dim) {dim= CG_numBCMessages_Fix3B[ch];}
      }
      dim=dim+1;
      tmp2= newArr2(RGBC_struct, numChildren, dim);
      for (int ch=0; ch<numChildren; ch++){
	for (int d=0; d<dim; d++){
	  tmp2[ch][d]= CG_Info_Fix3B[ch][d];
	}
      }
      delArr2(CG_Info_Fix3B, numChildren);
      CG_Info_Fix3B= newArr2(RGBC_struct, numChildren, dim);
      for (int ch=0; ch<numChildren; ch++){
	for (int d=0; d<dim; d++){
	  CG_Info_Fix3B[ch][d]= tmp2[ch][d]; //SF not here
	}
      }
      delArr2(tmp2, numChildren);
      
    } // end if (MLMD_fixBCenters){

  } // end if(numChildren >0)
}

void EMfields3D::Trim_RGBC_Vectors_II(VirtualTopology3D *vct){

  int dim;
  MPI_Comm CommToParent= vct->getCommToParent();
  RGBC_struct *tmp;
  RGBC_struct **tmp2;

  if (CommToParent != MPI_COMM_NULL){

    // II
    dim= RG_numBCMessages_II+1;
    tmp= new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      tmp[i]= RGBC_Info_II[i];
    }
    delete[] RGBC_Info_II;
    RGBC_Info_II = new RGBC_struct[dim];
    
    for (int i=0; i<dim; i++){
      RGBC_Info_II[i]= tmp[i];
    }
    delete[]tmp;
    
  }
  
  if (numChildren >0){

    dim=0;
    for (int ch=0; ch< numChildren; ch++){
      if (CG_numBCMessages_II[ch]>dim) {dim= CG_numBCMessages_II[ch];}
    }
    dim=dim+1;
    tmp2= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	tmp2[ch][d]= CG_Info_II[ch][d];
      }
    }
    delete[] CG_Info_II;
    CG_Info_II= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch<numChildren; ch++){
      for (int d=0; d<dim; d++){
	CG_Info_II[ch][d]= tmp2[ch][d]; //SF not here
      }
    }
    delete[] tmp2;

  } // end if(numChildren >0)
}

/* BC related operations */
/* sendBC: coarse grids sends BCs to the refined grids */
void EMfields3D::sendBC(Grid *grid, VirtualTopology3D *vct){

  // this is only for coarses grids 
  if (numChildren < 1){ return; }
  
  int msg; // this is the fake msg that will be sent to test the communication
  int rank_As_Parent;
  int dest; // destination on the parent-child communicator
  MPI_Comm PC_Comm; // the parent-child communicator
  int tag;

  // IMPORTANT: before sending BC, make sure ghost nodes are fixed
  //   This is enforced here, even if it should be already ok from the Field Solved  
  //   but DO NOT USE the communicateNode BC because you don't need to mess up with the BC 
  
  /*communicateNode(nxn, nyn, nzn, Ex, vct);
    communicateNode(nxn, nyn, nzn, Ey, vct);
    communicateNode(nxn, nyn, nzn, Ez, vct);
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);*/
  
  
  // RGBC_struct ** CG_Info_Ghost;
  //  int *CG_numBCMessages_Ghost;
  
  //  RGBC_struct ** CG_Info_Active;
  //  int *CG_numBCMessages_Active; 

  // I am sending ghost E, ghost Eth, ghost B
  // active E, active Eth, active B
  // at the moment i do not use active B
  

  for (int ch=0; ch< numChildren; ch ++){
    // active
    for (int m=0; m<CG_numBCMessages_Active[ch]; m++){
      sendOneBC(vct, grid, CG_Info_Active[ch][m], ch, 0);
    }

  }

  // THIS BARRIER HAS TO STAY!!!
  // to make sure all active msgs are sent before the ghost start being sent (otherwise the child, waiting for an active, may actually receive the ghost first)
  MPI_Barrier(vct->getComm()); 
  // END THIS BARRIER HAS TO STAY!!!

  for (int ch=0; ch< numChildren; ch ++){

    for (int m=0; m<CG_numBCMessages_Ghost[ch]; m++){
      sendOneBC(vct, grid, CG_Info_Ghost[ch][m], ch, -1);
    }
  } // end cycle on child

  if (MLMD_BCBufferArea ){
    MPI_Barrier(vct->getComm());
    // this barrier has to stay!!!
    for (int ch=0; ch< numChildren; ch ++){
      for (int m=0; m<CG_numBCMessages_Buffer[ch]; m++){
	sendOneBC(vct, grid, CG_Info_Buffer[ch][m], ch, 7);
	cout << "I am sending Buffer BC " << endl;
      }
    }  
  } // end if (MLMD_BCBufferArea ){

  // NB: for MLMD_InterpolateOldBCell I am re-using the MLMD_fixBCenters infrastructure
  if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){
    MPI_Barrier(vct->getComm());
    // this barrier has to stay!!!
    for (int ch=0; ch< numChildren; ch ++){
      for (int m=0; m<CG_numBCMessages_Fix3B[ch]; m++){
	sendOneBC(vct, grid, CG_Info_Fix3B[ch][m], ch, 9);
      }
    }  
  } // end if (MLMD_fixBCenters ){



  /*MPI_Barrier(vct->getComm());
  if (vct->getCartesian_rank()==0){
    cout << "grid " << numGrid << " Finished sending BC" << endl;
    }*/


}
/* receiveBC: refined grids receive BCs from the coarse grids */
void EMfields3D::receiveBC(Grid *grid, VirtualTopology3D *vct){

  MPI_Comm CommToParent= vct->getCommToParent();
  MPI_Comm CommToParent_BCGhost= vct->getCommToParent_BCGhost();
  MPI_Comm CommToParent_BCBuffer= vct->getCommToParent_BCBuffer();
  MPI_Comm CommToParent_BCFix3B= vct->getCommToParent_BCFix3B();
  if (CommToParent == MPI_COMM_NULL) {return;}  // if you are not a refined grid, no need to be here

  /* RGBC_struct * RGBC_Info_Ghost;
     int RG_numBCMessages_Ghost;

     RGBC_struct * RGBC_Info_Active;
     int RG_numBCMessages_Active; */
  //cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " trying to receive msg " <<endl;

  MPI_Status status;
  bool found;
  int countExp;
  int count;

  
  bool Testing= true; // put to false during production
  // NB: Msg instantiated together with E*_***_BC
  // active
  for (int m=0; m< RG_numBCMessages_Active; m++){

    MPI_Recv(RGMsg, RG_MaxMsgSize *NumF, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count );
    

    found= false;
    // check with stored msgs
    for (int i=0; i< RG_numBCMessages_Active; i++){

      if (RGBC_Info_Active[i].MsgID == status.MPI_TAG and RGBC_Info_Active[i].CG_core == status.MPI_SOURCE){ // NB: the tag gives you the msg position in the BC vector
	
	//cout << "R" << vct->getSystemWide_rank() <<" R" << vct->getRank_CommToParent() <<" on PC communicator, has received and matchedActive msg from " << status.MPI_SOURCE <<" with size " <<count << " and tag " <<status.MPI_TAG <<endl;

	found= true;
	  
	countExp= RGBC_Info_Active[i].np_x*RGBC_Info_Active[i].np_y*RGBC_Info_Active[i].np_z  ;

	if (Testing){
	  if (countExp *NumF  != count){
	    cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but Active size does not check: received: " << count << " expected " << countExp*NumF << ", aborting ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	}
	  
	// put BC here
	// the tag gives you the position in the BC vector
	
	// E
	memcpy(Ex_Active_BC[status.MPI_TAG], RGMsg, sizeof(double)*countExp );
	memcpy(Ey_Active_BC[status.MPI_TAG], RGMsg +   countExp, sizeof(double)*countExp );
	memcpy(Ez_Active_BC[status.MPI_TAG], RGMsg + 2*countExp, sizeof(double)*countExp );
	
	// Eth
	memcpy(Exth_Active_BC[status.MPI_TAG], RGMsg+ 3*countExp, sizeof(double)*countExp );
	memcpy(Eyth_Active_BC[status.MPI_TAG], RGMsg+ 4*countExp, sizeof(double)*countExp );
	memcpy(Ezth_Active_BC[status.MPI_TAG], RGMsg+ 5*countExp, sizeof(double)*countExp );
	
	// B
	memcpy(Bxn_Active_BC[status.MPI_TAG], RGMsg +6*countExp, sizeof(double)*countExp );
	memcpy(Byn_Active_BC[status.MPI_TAG], RGMsg +7*countExp, sizeof(double)*countExp );
	memcpy(Bzn_Active_BC[status.MPI_TAG], RGMsg +8*countExp, sizeof(double)*countExp );


	//
	if (false){
	  cout << "R" << vct->getSystemWide_rank() << " ACTIVE TAG " << status.MPI_TAG << " starts from [" << RGBC_Info_Active[i].ix_first <<"-" << RGBC_Info_Active[i].iy_first <<"-" << RGBC_Info_Active[i].iz_first <<"], Exth: " <<endl;

	  for (int j= 0 ; j< countExp; j++) cout << Ex_Active_BC[status.MPI_TAG][j] << " " ;
	  
	  cout << endl;

	} // end if true

	break;

      } // end if (RGBC_Info_Active[i].MsgID == status.MPI_TAG )
              
    } //for (int i=0; i< RG_numBCMessages_Active; i++){ // here end search of the msg

    if (Testing){
      if (found == false){
	cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received an active msg I cannot match with my record, aborting..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    
  } // here is when i receive it 

  MPI_Barrier(vct->getComm());
 

  // ghosts
  for (int m=0; m< RG_numBCMessages_Ghost; m++){

    MPI_Recv(RGMsg, RG_MaxMsgSize *NumF, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent_BCGhost, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count );
    
    found= false;
    // check with stored msgs
    for (int i=0; i< RG_numBCMessages_Ghost; i++){

      if (RGBC_Info_Ghost[i].MsgID == status.MPI_TAG and (RGBC_Info_Ghost[i].CG_core == status.MPI_SOURCE)){ // NB: the tag gives you the msg position in the BC vector
	
	found= true;

	//cout << "R" << vct->getSystemWide_rank() <<" R" << vct->getRank_CommToParent() <<" on PC communicator, has received and matched Ghost msg from " << status.MPI_SOURCE <<" with size " <<count << " and tag " <<status.MPI_TAG <<endl;
	  
	countExp= RGBC_Info_Ghost[i].np_x*RGBC_Info_Ghost[i].np_y*RGBC_Info_Ghost[i].np_z  ;

	if (Testing){
	  if (countExp *NumF  != count ){
	    cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but Ghost size does not check: received: " << count << " expected " << countExp*NumF << ", aborting ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	}
	  
	// put BC here
	// the tag gives you the position in the BC vector
	
	// E
	memcpy(Ex_Ghost_BC[status.MPI_TAG], RGMsg, sizeof(double)*countExp );
	memcpy(Ey_Ghost_BC[status.MPI_TAG], RGMsg +   countExp, sizeof(double)*countExp );
	memcpy(Ez_Ghost_BC[status.MPI_TAG], RGMsg + 2*countExp, sizeof(double)*countExp );
	
	// Eth
	memcpy(Exth_Ghost_BC[status.MPI_TAG], RGMsg+ 3*countExp, sizeof(double)*countExp );
	memcpy(Eyth_Ghost_BC[status.MPI_TAG], RGMsg+ 4*countExp, sizeof(double)*countExp );
	memcpy(Ezth_Ghost_BC[status.MPI_TAG], RGMsg+ 5*countExp, sizeof(double)*countExp );
	
	// B
	memcpy(Bxn_Ghost_BC[status.MPI_TAG], RGMsg +6*countExp, sizeof(double)*countExp );
	memcpy(Byn_Ghost_BC[status.MPI_TAG], RGMsg +7*countExp, sizeof(double)*countExp );
	memcpy(Bzn_Ghost_BC[status.MPI_TAG], RGMsg +8*countExp, sizeof(double)*countExp );


	if (false){
	  cout << "R" << vct->getSystemWide_rank() << " GHOST TAG " << status.MPI_TAG << " starts from [" << RGBC_Info_Ghost[i].ix_first <<"-" << RGBC_Info_Ghost[i].iy_first <<"-" << RGBC_Info_Ghost[i].iz_first<< "], Exth: " <<endl;

	  for (int j= 0 ; j< countExp; j++) cout << Ex_Ghost_BC[status.MPI_TAG][j] << " " ;
	  
	  cout << endl;

	} // end if true
	break;
      } // end if (RGBC_Info_Ghost[i].MsgID == status.MPI_TAG )
           
    } // end search msg

    if (Testing){
      if (found == false){
	cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received a ghostmsg I cannot match with my record, aborting..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
  } // end receive msg

  

  // buffer
  if (MLMD_BCBufferArea ){
    cout << "I am receiving Buffer BC" << endl;
    MPI_Barrier(vct->getComm());
    for (int m=0; m< RG_numBCMessages_Buffer; m++){
      
      MPI_Recv(RGMsgBuffer, RG_MaxMsgBufferSize *NumF, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent_BCBuffer, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &count );
      
      found= false;
      // check with stored msgs
      for (int i=0; i< RG_numBCMessages_Buffer; i++){

	if (RGBC_Info_Buffer[i].MsgID == status.MPI_TAG and (RGBC_Info_Buffer[i].CG_core == status.MPI_SOURCE)){ // NB: the tag gives you the msg position in the BC vector
	  
	  found= true;
	  
	  //cout << "R" << vct->getSystemWide_rank() <<" R" << vct->getRank_CommToParent() <<" on PC communicator, has received and matched buffer msg from " << status.MPI_SOURCE <<" with size " <<count << " and tag " <<status.MPI_TAG <<endl;
	  
	  countExp= RGBC_Info_Buffer[i].np_x*RGBC_Info_Buffer[i].np_y*RGBC_Info_Buffer[i].np_z  ;

	  if (Testing){
	    if (countExp *NumF  != count ){
	      cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but Buffer size does not check: received: " << count << " expected " << countExp*NumF << ", aborting ..." << endl;
	      MPI_Abort(MPI_COMM_WORLD, -1);
	    }
	  }
	  
	  // put BC here
	  // the tag gives you the position in the BC vector
	  
	  // E
	  memcpy(Ex_Buffer_BC[status.MPI_TAG], RGMsgBuffer, sizeof(double)*countExp );
	  memcpy(Ey_Buffer_BC[status.MPI_TAG], RGMsgBuffer +   countExp, sizeof(double)*countExp );
	  memcpy(Ez_Buffer_BC[status.MPI_TAG], RGMsgBuffer + 2*countExp, sizeof(double)*countExp );
	  
	  // Eth
	  memcpy(Exth_Buffer_BC[status.MPI_TAG], RGMsgBuffer+ 3*countExp, sizeof(double)*countExp );
	  memcpy(Eyth_Buffer_BC[status.MPI_TAG], RGMsgBuffer+ 4*countExp, sizeof(double)*countExp );
	  memcpy(Ezth_Buffer_BC[status.MPI_TAG], RGMsgBuffer+ 5*countExp, sizeof(double)*countExp );
	  
	  // B
	  memcpy(Bxn_Buffer_BC[status.MPI_TAG], RGMsgBuffer +6*countExp, sizeof(double)*countExp );
	  memcpy(Byn_Buffer_BC[status.MPI_TAG], RGMsgBuffer +7*countExp, sizeof(double)*countExp );
	  memcpy(Bzn_Buffer_BC[status.MPI_TAG], RGMsgBuffer +8*countExp, sizeof(double)*countExp );
	  
	  
	  if (false){
	    cout << "R" << vct->getSystemWide_rank() << " Buffer TAG " << status.MPI_TAG << " starts from [" << RGBC_Info_Buffer[i].ix_first <<"-" << RGBC_Info_Buffer[i].iy_first <<"-" << RGBC_Info_Buffer[i].iz_first<< "], Exth: " <<endl;
	    
	    for (int j= 0 ; j< countExp; j++) cout << Ex_Buffer_BC[status.MPI_TAG][j] << " " ;
	    
	    cout << endl;
	    
	  } // end if true
	  break;
	} // end if (RGBC_Info_Ghost[i].MsgID == status.MPI_TAG )
	
      } // end search msg
      
      if (Testing){
	if (found == false){
	  cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received a buffer msg I cannot match with my record, aborting..." << endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);
	}
      }
    } // end receive msg
    
  } // end if (MLMD_BCBufferArea)

  /** fix3B **/
  // i am reusing MLMD_fixBCenters infrastructure for MLMD_InterpolateOldBCell
  if (MLMD_fixBCenters or MLMD_InterpolateOldBCell){

    MPI_Barrier(vct->getComm());
    for (int m=0; m< RG_numBCMessages_fix3B; m++){
      
      MPI_Recv(RGMsgfix3B, RG_Maxfix3BMsgSize *Numfix3B, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent_BCFix3B, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &count );
      
      found= false;
      // check with stored msgs
      for (int i=0; i< RG_numBCMessages_fix3B; i++){

	if (RGBC_Info_fix3B[i].MsgID == status.MPI_TAG and (RGBC_Info_fix3B[i].CG_core == status.MPI_SOURCE)){ // NB: the tag gives you the msg position in the BC vector
	  
	  found= true;
	  
	  //cout << "R" << vct->getSystemWide_rank() <<" R" << vct->getRank_CommToParent() <<" on PC communicator, has received and matched buffer msg from " << status.MPI_SOURCE <<" with size " <<count << " and tag " <<status.MPI_TAG <<endl;
	  
	  countExp= RGBC_Info_fix3B[i].np_x* RGBC_Info_fix3B[i].np_y* RGBC_Info_fix3B[i].np_z  ;

	  if (Testing){
	    if (countExp *Numfix3B  != count ){
	      cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but Fix3B size does not check: received: " << count << " expected " << countExp*Numfix3B << ", aborting ..." << endl;
	      MPI_Abort(MPI_COMM_WORLD, -1);
	    }
	  }
	  
	  // put BC here
	  // the tag gives you the position in the BC vector
	  
	  // B Centers
	  memcpy(Bxc_fix3B_BC[status.MPI_TAG], RGMsgfix3B, sizeof(double)*countExp );
	  memcpy(Byc_fix3B_BC[status.MPI_TAG], RGMsgfix3B +   countExp, sizeof(double)*countExp );
	  memcpy(Bzc_fix3B_BC[status.MPI_TAG], RGMsgfix3B + 2*countExp, sizeof(double)*countExp );
	  
	  if (false){
	    cout << "R" << vct->getSystemWide_rank() << " Fic3B TAG " << status.MPI_TAG << " starts from [" << RGBC_Info_fix3B[i].ix_first <<"-" << RGBC_Info_fix3B[i].iy_first <<"-" << RGBC_Info_fix3B[i].iz_first<< "], Exth: " <<endl;
	    
	    
	    cout << endl;
	    
	  } // end if true
	  break;
	} // end if (RGBC_Info_Fix3B[i].MsgID == status.MPI_TAG )
	
      } // end search msg
      
      if (Testing){
	if (found == false){
	  cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received a buffer msg I cannot match with my record, aborting..." << endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);
	}
      }
    } // end receive msg
    
  } // end if (MLMD_fixBCenters)

  /** end fix3B **/



}

void EMfields3D::receiveInitialInterpolation(Grid *grid, VirtualTopology3D *vct){

  MPI_Comm CommToParent= vct->getCommToParent();

  if (CommToParent == MPI_COMM_NULL) {return;}  // if you are not a refined grid, no need to be here

  MPI_Status status;
  bool found;
  int countExp;
  int count;

  
  bool Testing= true; // put to false during production
  // NB: Msg instantiated together with E*_***_BC
  // active
  for (int m=0; m< RG_numBCMessages_II; m++){

    MPI_Recv(RGMsg_II, RG_MaxMsgSize_II *NumF_II, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count );
    
    found= false;
    // check with stored msgs
    for (int i=0; i< RG_numBCMessages_II; i++){

      if (RGBC_Info_II[i].MsgID == status.MPI_TAG and RGBC_Info_II[i].CG_core == status.MPI_SOURCE){ // NB: the tag gives you the msg position in the BC vector

	found= true;
	  
	countExp= RGBC_Info_II[i].np_x*RGBC_Info_II[i].np_y*RGBC_Info_II[i].np_z  ;

	if (Testing){
	  if (countExp *NumF_II  != count){
	    cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but II size does not check: received: " << count << " expected " << countExp*NumF_II << ", aborting ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	}
	  
	//cout << "RG_numBCMessages_II: " << RG_numBCMessages_II << " countExp: " <<countExp <<endl;

	/*for (int ii=0; ii< countExp; ii++){
	  cout << "In receiveII, RGMsg_II[ii]: " << RGMsg_II[ii] << endl;
	  }*/

	for (int cc=0; cc< countExp; cc++){
	  Ex_II[status.MPI_TAG][cc]= RGMsg_II[cc];
	  Ey_II[status.MPI_TAG][cc]= RGMsg_II[countExp + cc];
	  Ez_II[status.MPI_TAG][cc]= RGMsg_II[2*countExp + cc];

	  Bxn_II[status.MPI_TAG][cc]= RGMsg_II[3*countExp + cc];
	  Byn_II[status.MPI_TAG][cc]= RGMsg_II[4*countExp + cc];
	  Bzn_II[status.MPI_TAG][cc]= RGMsg_II[5*countExp + cc];

	  if (ns>0)
	    rho0_II[status.MPI_TAG][cc]= RGMsg_II[6*countExp + cc];
	  if (ns>1)
	    rho1_II[status.MPI_TAG][cc]= RGMsg_II[7*countExp + cc];
	  if (ns>2)
	    rho2_II[status.MPI_TAG][cc]= RGMsg_II[8*countExp + cc];
	  if (ns>3)
	    rho3_II[status.MPI_TAG][cc]= RGMsg_II[9*countExp + cc];
	}

	break;

      } // end if (RGBC_Info_Active[i].MsgID == status.MPI_TAG )
    
    } //for (int i=0; i< RG_numBCMessages_Active; i++){ // here end search of the msg

    if (Testing){
      if (found == false){
	cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received an II msg I cannot match with my record, aborting..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    
  } // here is when i receive it 

  MPI_Barrier(vct->getComm());

}

void EMfields3D::buildBCMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo, int Size, double *Msg ){
  
  // NB: if i send more stuff, do all the interpolation here
  // i use the three indexes, knowing that one will be zero

  // position of the initial RG point in the CG grid
  // points are displacement over this

  //cout << "Building msg starting from " << RGInfo.CG_x_first <<", " << RGInfo.CG_y_first << ", " <<RGInfo.CG_z_first <<" size "<< RGInfo.np_x << ", "  << RGInfo.np_y <<", " << RGInfo.np_z << endl;

  double x0= RGInfo.CG_x_first;
  double y0= RGInfo.CG_y_first;
  double z0= RGInfo.CG_z_first;

  double inv_dx= 1./dx;
  double inv_dy= 1./dy;
  double inv_dz= 1./dz;

  double xp, yp, zp;
  int count =0;

  //cout << "inside building msg RGInfo.np_x " << RGInfo.np_x << " RGInfo.np_y " << RGInfo.np_y << " RGInfo.np_z " << RGInfo.np_z <<endl;

  //cout <<"building msg, starts @ " << x0 << "; " << y0 << "; " << z0 << endl;

  for (int i= 0; i<RGInfo.np_x; i++){
    for (int j=0; j<RGInfo.np_y; j++){
      for (int k=0; k<RGInfo.np_z; k++){
	
	// ok, now each RG point is treated as a particle here
	// and i need the E at the point position

	xp= x0 + i*dx_Ch[ch];
	yp= y0 + j*dy_Ch[ch];
	zp= z0 + k*dz_Ch[ch]; 
	
	//cout << "building msg, point " << xp << "; " <<yp <<"; " << zp << endl;

	// this is copied from the mover
	const double ixd = floor((xp - xStart) * inv_dx);
	const double iyd = floor((yp - yStart) * inv_dy);
	const double izd = floor((zp - zStart) * inv_dz);
	int ix = 2 + int (ixd);
	int iy = 2 + int (iyd);
	int iz = 2 + int (izd);
	if (ix < 1)
	  ix = 1;
	if (iy < 1)
	  iy = 1;
	if (iz < 1)
	  iz = 1;
	if (ix > nxn - 1)
	  ix = nxn - 1;
	if (iy > nyn - 1)
	  iy = nyn - 1;
	if (iz > nzn - 1)
	  iz = nzn - 1;

	double xi  [2];
	double eta [2];
	double zeta[2];

	xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
	eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
	zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
	xi  [1] = grid->getXN(ix,iy,iz) - xp;
	eta [1] = grid->getYN(ix,iy,iz) - yp;
	zeta[1] = grid->getZN(ix,iy,iz) - zp;

	/*cout << "point: " << xp <<", " << yp <<", " << zp << endl;
	cout << "xi[0] " << xi[0] << " xi[1] " << xi[1] << " eta[0] " << eta[0] << " eta[1] " << eta[1] << " zeta[0] " << zeta[0] << " zeta[1] " << zeta[1];
	cout << "invVOL: " << invVOL << endl;*/

	double Exl = 0.0;
	double Eyl = 0.0;
	double Ezl = 0.0;

	double Exthl = 0.0;
	double Eythl = 0.0;
	double Ezthl = 0.0;

	double Bxl = 0.0;
	double Byl = 0.0;
	double Bzl = 0.0;

	const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
	const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
	const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
	const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
	const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
	const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
	const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
	const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
	
	//                                                                                                                                                       
	Exl += weight000 * (Ex[ix][iy][iz]             );
	Exl += weight001 * (Ex[ix][iy][iz - 1]         );
	Exl += weight010 * (Ex[ix][iy - 1][iz]         );
	Exl += weight011 * (Ex[ix][iy - 1][iz - 1]     );
	Exl += weight100 * (Ex[ix - 1][iy][iz]         );
	Exl += weight101 * (Ex[ix - 1][iy][iz - 1]     );
	Exl += weight110 * (Ex[ix - 1][iy - 1][iz]     );
	Exl += weight111 * (Ex[ix - 1][iy - 1][iz - 1] );

	Eyl += weight000 * (Ey[ix][iy][iz]             );
	Eyl += weight001 * (Ey[ix][iy][iz - 1]         );
	Eyl += weight010 * (Ey[ix][iy - 1][iz]         );
	Eyl += weight011 * (Ey[ix][iy - 1][iz - 1]     );
	Eyl += weight100 * (Ey[ix - 1][iy][iz]         );
	Eyl += weight101 * (Ey[ix - 1][iy][iz - 1]     );
	Eyl += weight110 * (Ey[ix - 1][iy - 1][iz]     );
	Eyl += weight111 * (Ey[ix - 1][iy - 1][iz - 1] );

	Ezl += weight000 * (Ez[ix][iy][iz]             );
	Ezl += weight001 * (Ez[ix][iy][iz - 1]         );
	Ezl += weight010 * (Ez[ix][iy - 1][iz]         );
	Ezl += weight011 * (Ez[ix][iy - 1][iz - 1]     );
	Ezl += weight100 * (Ez[ix - 1][iy][iz]         );
	Ezl += weight101 * (Ez[ix - 1][iy][iz - 1]     );
	Ezl += weight110 * (Ez[ix - 1][iy - 1][iz]     );
	Ezl += weight111 * (Ez[ix - 1][iy - 1][iz - 1] ); 

	/*Exl += weight000 * (Ex_BS[ix][iy][iz]             );
        Exl += weight001 * (Ex_BS[ix][iy][iz - 1]         );
        Exl += weight010 * (Ex_BS[ix][iy - 1][iz]         );
        Exl += weight011 * (Ex_BS[ix][iy - 1][iz - 1]     );
        Exl += weight100 * (Ex_BS[ix - 1][iy][iz]         );
        Exl += weight101 * (Ex_BS[ix - 1][iy][iz - 1]     );
        Exl += weight110 * (Ex_BS[ix - 1][iy - 1][iz]     );
        Exl += weight111 * (Ex_BS[ix - 1][iy - 1][iz - 1] );

        Eyl += weight000 * (Ey_BS[ix][iy][iz]             );
        Eyl += weight001 * (Ey_BS[ix][iy][iz - 1]         );
        Eyl += weight010 * (Ey_BS[ix][iy - 1][iz]         );
        Eyl += weight011 * (Ey_BS[ix][iy - 1][iz - 1]     );
        Eyl += weight100 * (Ey_BS[ix - 1][iy][iz]         );
        Eyl += weight101 * (Ey_BS[ix - 1][iy][iz - 1]     );
        Eyl += weight110 * (Ey_BS[ix - 1][iy - 1][iz]     );
        Eyl += weight111 * (Ey_BS[ix - 1][iy - 1][iz - 1] );

        Ezl += weight000 * (Ez_BS[ix][iy][iz]             );
        Ezl += weight001 * (Ez_BS[ix][iy][iz - 1]         );
        Ezl += weight010 * (Ez_BS[ix][iy - 1][iz]         );
        Ezl += weight011 * (Ez_BS[ix][iy - 1][iz - 1]     );
        Ezl += weight100 * (Ez_BS[ix - 1][iy][iz]         );
        Ezl += weight101 * (Ez_BS[ix - 1][iy][iz - 1]     );
        Ezl += weight110 * (Ez_BS[ix - 1][iy - 1][iz]     );
        Ezl += weight111 * (Ez_BS[ix - 1][iy - 1][iz - 1] ); */


	// TH
	Exthl += weight000 * (Exth[ix][iy][iz]             );
	Exthl += weight001 * (Exth[ix][iy][iz - 1]         );
	Exthl += weight010 * (Exth[ix][iy - 1][iz]         );
	Exthl += weight011 * (Exth[ix][iy - 1][iz - 1]     );
	Exthl += weight100 * (Exth[ix - 1][iy][iz]         );
	Exthl += weight101 * (Exth[ix - 1][iy][iz - 1]     );
	Exthl += weight110 * (Exth[ix - 1][iy - 1][iz]     );
	Exthl += weight111 * (Exth[ix - 1][iy - 1][iz - 1] );

	Eythl += weight000 * (Eyth[ix][iy][iz]             );
	Eythl += weight001 * (Eyth[ix][iy][iz - 1]         );
	Eythl += weight010 * (Eyth[ix][iy - 1][iz]         );
	Eythl += weight011 * (Eyth[ix][iy - 1][iz - 1]     );
	Eythl += weight100 * (Eyth[ix - 1][iy][iz]         );
	Eythl += weight101 * (Eyth[ix - 1][iy][iz - 1]     );
	Eythl += weight110 * (Eyth[ix - 1][iy - 1][iz]     );
	Eythl += weight111 * (Eyth[ix - 1][iy - 1][iz - 1] );

	Ezthl += weight000 * (Ezth[ix][iy][iz]             );
	Ezthl += weight001 * (Ezth[ix][iy][iz - 1]         );
	Ezthl += weight010 * (Ezth[ix][iy - 1][iz]         );
	Ezthl += weight011 * (Ezth[ix][iy - 1][iz - 1]     );
	Ezthl += weight100 * (Ezth[ix - 1][iy][iz]         );
	Ezthl += weight101 * (Ezth[ix - 1][iy][iz - 1]     );
	Ezthl += weight110 * (Ezth[ix - 1][iy - 1][iz]     );
	Ezthl += weight111 * (Ezth[ix - 1][iy - 1][iz - 1] );
	// end TH

	// B
	Bxl += weight000 * (Bxn[ix][iy][iz]             );
	Bxl += weight001 * (Bxn[ix][iy][iz - 1]         );
	Bxl += weight010 * (Bxn[ix][iy - 1][iz]         );
	Bxl += weight011 * (Bxn[ix][iy - 1][iz - 1]     );
	Bxl += weight100 * (Bxn[ix - 1][iy][iz]         );
	Bxl += weight101 * (Bxn[ix - 1][iy][iz - 1]     );
	Bxl += weight110 * (Bxn[ix - 1][iy - 1][iz]     );
	Bxl += weight111 * (Bxn[ix - 1][iy - 1][iz - 1] );

	Byl += weight000 * (Byn[ix][iy][iz]             );
	Byl += weight001 * (Byn[ix][iy][iz - 1]         );
	Byl += weight010 * (Byn[ix][iy - 1][iz]         );
	Byl += weight011 * (Byn[ix][iy - 1][iz - 1]     );
	Byl += weight100 * (Byn[ix - 1][iy][iz]         );
	Byl += weight101 * (Byn[ix - 1][iy][iz - 1]     );
	Byl += weight110 * (Byn[ix - 1][iy - 1][iz]     );
	Byl += weight111 * (Byn[ix - 1][iy - 1][iz - 1] );

	Bzl += weight000 * (Bzn[ix][iy][iz]             );
	Bzl += weight001 * (Bzn[ix][iy][iz - 1]         );
	Bzl += weight010 * (Bzn[ix][iy - 1][iz]         );
	Bzl += weight011 * (Bzn[ix][iy - 1][iz - 1]     );
	Bzl += weight100 * (Bzn[ix - 1][iy][iz]         );
	Bzl += weight101 * (Bzn[ix - 1][iy][iz - 1]     );
	Bzl += weight110 * (Bzn[ix - 1][iy - 1][iz]     );
	Bzl += weight111 * (Bzn[ix - 1][iy - 1][iz - 1] );
	// end B
	
	Msg[0*Size +count]= Exl;
	Msg[1*Size +count]= Eyl;
	Msg[2*Size +count]= Ezl;

	Msg[3*Size +count]= Exthl;
	Msg[4*Size +count]= Eythl;
	Msg[5*Size +count]= Ezthl;

	Msg[6*Size +count]= Bxl;
	Msg[7*Size +count]= Byl;
	Msg[8*Size +count]= Bzl;

	count ++;

      }
    }
  }
  //cout <<"R" << vct->getSystemWide_rank() << ", count " << count << "for child " << ch << endl;
}

/* fix3B */
void EMfields3D::buildFix3BMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo, int Size, double *Msg ){
  
  // NB: if i send more stuff, do all the interpolation here
  // i use the three indexes, knowing that one will be zero

  // position of the initial RG point in the CG grid
  // points are displacement over this

  //cout << "Building msg starting from " << RGInfo.CG_x_first <<", " << RGInfo.CG_y_first << ", " <<RGInfo.CG_z_first <<" size "<< RGInfo.np_x << ", "  << RGInfo.np_y <<", " << RGInfo.np_z << endl;

  double x0= RGInfo.CG_x_first;
  double y0= RGInfo.CG_y_first;
  double z0= RGInfo.CG_z_first;

  double inv_dx= 1./dx;
  double inv_dy= 1./dy;
  double inv_dz= 1./dz;

  double xp, yp, zp;
  int count =0;

  //cout << "inside building msg RGInfo.np_x " << RGInfo.np_x << " RGInfo.np_y " << RGInfo.np_y << " RGInfo.np_z " << RGInfo.np_z <<endl;

  //cout <<"building msg, starts @ " << x0 << "; " << y0 << "; " << z0 << endl;

  for (int i= 0; i<RGInfo.np_x; i++){
    for (int j=0; j<RGInfo.np_y; j++){
      for (int k=0; k<RGInfo.np_z; k++){
	
	// ok, now each RG point is treated as a particle here
	// and i need the E at the point position

	xp= x0 + i*dx_Ch[ch];
	yp= y0 + j*dy_Ch[ch];
	zp= z0 + k*dz_Ch[ch]; 
	

	// this is done centers to centers
	/*// here, xp, yp, zp is the position of the CENTERS

	//cout << "building msg, point " << xp << "; " <<yp <<"; " << zp << endl;

	// this is copied from the mover
	//const double ixd = floor((xp - xStart) * inv_dx);
	//const double iyd = floor((yp - yStart) * inv_dy);
	//const double izd = floor((zp - zStart) * inv_dz);
	const double ixd = floor((xp - (xStart +dx )) * inv_dx);
	const double iyd = floor((yp - (yStart +dy )) * inv_dy);
	const double izd = floor((zp - (zStart +dz )) * inv_dz);
	int ix = 2 + int (ixd);
	int iy = 2 + int (iyd);
	int iz = 2 + int (izd);
	if (ix < 1)
	  ix = 1;
	if (iy < 1)
	  iy = 1;
	if (iz < 1)
	  iz = 1;
	if (ix > nxc - 1)
	  ix = nxc - 1;
	if (iy > nyc - 1)
	  iy = nyc - 1;
	if (iz > nzc - 1)
	  iz = nzc - 1;

	double xi  [2];
	double eta [2];
	double zeta[2];

	xi  [0] = xp - grid->getXC(ix-1,iy  ,iz  );
	eta [0] = yp - grid->getYC(ix  ,iy-1,iz  );
	zeta[0] = zp - grid->getZC(ix  ,iy  ,iz-1);
	xi  [1] = grid->getXC(ix,iy,iz) - xp;
	eta [1] = grid->getYC(ix,iy,iz) - yp;
	zeta[1] = grid->getZC(ix,iy,iz) - zp;

	//cout << "point: " << xp <<", " << yp <<", " << zp << endl;
	//cout << "xi[0] " << xi[0] << " xi[1] " << xi[1] << " eta[0] " << eta[0] << " eta[1] " << eta[1] << " zeta[0] " << zeta[0] << " zeta[1] " << zeta[1];
	//cout << "invVOL: " << invVOL << endl;
	
	double Bxl = 0.0;
	double Byl = 0.0;
	double Bzl = 0.0;

	const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
	const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
	const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
	const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
	const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
	const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
	const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
	const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
	
	// B
	Bxl += weight000 * (Bxc[ix][iy][iz]             );
	Bxl += weight001 * (Bxc[ix][iy][iz - 1]         );
	Bxl += weight010 * (Bxc[ix][iy - 1][iz]         );
	Bxl += weight011 * (Bxc[ix][iy - 1][iz - 1]     );
	Bxl += weight100 * (Bxc[ix - 1][iy][iz]         );
	Bxl += weight101 * (Bxc[ix - 1][iy][iz - 1]     );
	Bxl += weight110 * (Bxc[ix - 1][iy - 1][iz]     );
	Bxl += weight111 * (Bxc[ix - 1][iy - 1][iz - 1] );

	Byl += weight000 * (Byc[ix][iy][iz]             );
	Byl += weight001 * (Byc[ix][iy][iz - 1]         );
	Byl += weight010 * (Byc[ix][iy - 1][iz]         );
	Byl += weight011 * (Byc[ix][iy - 1][iz - 1]     );
	Byl += weight100 * (Byc[ix - 1][iy][iz]         );
	Byl += weight101 * (Byc[ix - 1][iy][iz - 1]     );
	Byl += weight110 * (Byc[ix - 1][iy - 1][iz]     );
	Byl += weight111 * (Byc[ix - 1][iy - 1][iz - 1] );

	Bzl += weight000 * (Bzc[ix][iy][iz]             );
	Bzl += weight001 * (Bzc[ix][iy][iz - 1]         );
	Bzl += weight010 * (Bzc[ix][iy - 1][iz]         );
	Bzl += weight011 * (Bzc[ix][iy - 1][iz - 1]     );
	Bzl += weight100 * (Bzc[ix - 1][iy][iz]         );
	Bzl += weight101 * (Bzc[ix - 1][iy][iz - 1]     );
	Bzl += weight110 * (Bzc[ix - 1][iy - 1][iz]     );
	Bzl += weight111 * (Bzc[ix - 1][iy - 1][iz - 1] );
	// end B*/

	// this on nodes, still to test
	// this is copied from the mover                                          
        const double ixd = floor((xp - xStart) * inv_dx);
        const double iyd = floor((yp - yStart) * inv_dy);
        const double izd = floor((zp - zStart) * inv_dz);
        int ix = 2 + int (ixd);
        int iy = 2 + int (iyd);
        int iz = 2 + int (izd);
        if (ix < 1)
          ix = 1;
        if (iy < 1)
          iy = 1;
        if (iz < 1)
          iz = 1;
        if (ix > nxn - 1)
          ix = nxn - 1;
        if (iy > nyn - 1)
          iy = nyn - 1;
        if (iz > nzn - 1)
          iz = nzn - 1;

        double xi  [2];
        double eta [2];
        double zeta[2];

	xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
        eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
        zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
        xi  [1] = grid->getXN(ix,iy,iz) - xp;
        eta [1] = grid->getYN(ix,iy,iz) - yp;
        zeta[1] = grid->getZN(ix,iy,iz) - zp;

        /*cout << "point: " << xp <<", " << yp <<", " << zp << endl;              
	  cout << "xi[0] " << xi[0] << " xi[1] " << xi[1] << " eta[0] " << eta[0] <\
	  < " eta[1] " << eta[1] << " zeta[0] " << zeta[0] << " zeta[1] " << zeta[1];       
	  cout << "invVOL: " << invVOL << endl;*/

        double Bxl = 0.0;
        double Byl = 0.0;
        double Bzl = 0.0;

	const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
        const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
        const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
        const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
        const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
        const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
        const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
        const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;

	// B                                                                      
        Bxl += weight000 * (Bxn[ix][iy][iz]             );
        Bxl += weight001 * (Bxn[ix][iy][iz - 1]         );
        Bxl += weight010 * (Bxn[ix][iy - 1][iz]         );
        Bxl += weight011 * (Bxn[ix][iy - 1][iz - 1]     );
        Bxl += weight100 * (Bxn[ix - 1][iy][iz]         );
        Bxl += weight101 * (Bxn[ix - 1][iy][iz - 1]     );
        Bxl += weight110 * (Bxn[ix - 1][iy - 1][iz]     );
        Bxl += weight111 * (Bxn[ix - 1][iy - 1][iz - 1] );

        Byl += weight000 * (Byn[ix][iy][iz]             );
        Byl += weight001 * (Byn[ix][iy][iz - 1]         );
        Byl += weight010 * (Byn[ix][iy - 1][iz]         );
        Byl += weight011 * (Byn[ix][iy - 1][iz - 1]     );
        Byl += weight100 * (Byn[ix - 1][iy][iz]         );
        Byl += weight101 * (Byn[ix - 1][iy][iz - 1]     );
        Byl += weight110 * (Byn[ix - 1][iy - 1][iz]     );
        Byl += weight111 * (Byn[ix - 1][iy - 1][iz - 1] );

        Bzl += weight000 * (Bzn[ix][iy][iz]             );
        Bzl += weight001 * (Bzn[ix][iy][iz - 1]         );
        Bzl += weight010 * (Bzn[ix][iy - 1][iz]         );
        Bzl += weight011 * (Bzn[ix][iy - 1][iz - 1]     );
        Bzl += weight100 * (Bzn[ix - 1][iy][iz]         );
        Bzl += weight101 * (Bzn[ix - 1][iy][iz - 1]     );
        Bzl += weight110 * (Bzn[ix - 1][iy - 1][iz]     );
        Bzl += weight111 * (Bzn[ix - 1][iy - 1][iz - 1] );
	
	Msg[0*Size +count]= Bxl;
	Msg[1*Size +count]= Byl;
	Msg[2*Size +count]= Bzl;

	count ++;

      }
    }
  }
  //cout <<"R" << vct->getSystemWide_rank() << ", count " << count << "for child " << ch << endl;

  if (MLMD_InterpolateOldBCell)
    cout << "I have successfully used fix3B structure to build interpolateB msg " << endl;

}

/* end fix3B */



/* the RG sets the received BC                                                                                              
       Fx, Fy, Fz: the field where BC are applied                      
       Fx_BC, Fy_BC, Fz_BC: the BCs              
       RGBC_Info: the RGBC_Info struct (Active or Ghost)                   
       RG_numBCMessages: number of messages (Active or Ghost) */
void EMfields3D::setBC_Nodes(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){

  //  sumEzBC=0.0;

  for (int m=0; m<RG_numBCMessages; m++){
    // points are included: <=
    // remember one index is 0
    // check here actually

    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){
	  
	  //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize); 	  
	  Fx[ix_f + i][iy_f + j][iz_f + k]= Fx_BC[m][count];
	  Fy[ix_f + i][iy_f + j][iz_f + k]= Fy_BC[m][count];
	  Fz[ix_f + i][iy_f + j][iz_f + k]= Fz_BC[m][count];
	  
	  /*sumEzBC+= fabs(Fz[ix_f + i][iy_f + j][iz_f + k]);
	  //cout << "E th BC: " << Fz[ix_f + i][iy_f + j][iz_f + k] << endl;
	  if (sumEzBC > DBL_EPSILON) {
	    MPI_Abort(MPI_COMM_WORLD, -1);
	    }*/

	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  
  /*if (sumEzBC>0) {
    cout << "Ez BC for light wave should be zero, fatal error, exithing";
    MPI_Abort(MPI_COMM_WORLD, -1);
    }*/
  //cout <<"sumEzBC= " << sumEzBC << endl;

  return;
}

void EMfields3D::setBC_Nodes_TwoLess(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){
  for (int m=0; m<RG_numBCMessages; m++){
    // points are included: <=
    // remember one index is 0
    // check here actually

    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  if (DirBuffer[m]=='L' and i>2 + BufX -2) continue; // left
	  if (DirBuffer[m]=='B' and k>2 + BufZ -2) continue; // bottom
	  if (DirBuffer[m]=='F' and j>2 + BufY -2) continue; // front
	  if (DirBuffer[m]=='R' and i< nxn-3-BufX +2) continue; // right
	  if (DirBuffer[m]=='T' and k< nzn-3-BufZ +2) continue; // top
	  if (DirBuffer[m]=='b' and j< nyn-3-BufY +2) continue; // back

	  //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize); 	  
	  Fx[ix_f + i][iy_f + j][iz_f + k]= Fx_BC[m][count];
	  Fy[ix_f + i][iy_f + j][iz_f + k]= Fy_BC[m][count];
	  Fz[ix_f + i][iy_f + j][iz_f + k]= Fz_BC[m][count];
	  
	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

/* to average RG and CG interpolated solution
       Fx, Fy, Fz: the field where BC are applied                      
       Fx_BC, Fy_BC, Fz_BC: the BCs              
       RGBC_Info: the RGBC_Info struct (Active or Ghost)                   
       RG_numBCMessages: number of messages (Active or Ghost) */
void EMfields3D::AverageBC_Nodes(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){

  double Fac;

  double *** CopyFx= newArr3(double, nxn, nyn, nzn);
  double *** CopyFy= newArr3(double, nxn, nyn, nzn);
  double *** CopyFz= newArr3(double, nxn, nyn, nzn);
  
  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++){
	CopyFx[i][j][k]= Fx[i][j][k];
	CopyFy[i][j][k]= Fy[i][j][k];
	CopyFz[i][j][k]= Fz[i][j][k];
      }
  
  for (int m=0; m<RG_numBCMessages; m++){
    // points are included: <=
    // remember one index is 0
    // check here actually

  
    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  Fac= BufferFactor(DirBuffer[m], ix_f + i, iy_f + j, iz_f + k, BufX, BufY, BufZ);
	  	  
	  //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize); 	  
	  Fx[ix_f + i][iy_f + j][iz_f + k]= Fac*Fx_BC[m][count] + (1.-Fac)* CopyFx[ix_f + i][iy_f + j][iz_f + k];
	  Fy[ix_f + i][iy_f + j][iz_f + k]= Fac*Fy_BC[m][count] + (1.-Fac)* CopyFy[ix_f + i][iy_f + j][iz_f + k];
	  Fz[ix_f + i][iy_f + j][iz_f + k]= Fac*Fz_BC[m][count] + (1.-Fac)* CopyFz[ix_f + i][iy_f + j][iz_f + k];
	  
	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  delArr3(CopyFx, nxn, nyn);
  delArr3(CopyFy, nxn, nyn);
  delArr3(CopyFz, nxn, nyn);

  return;
}

void EMfields3D::setBC_Nodes(VirtualTopology3D * vct, double **Fx_BC, double **Fy_BC, double **Fz_BC, double **Fa_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){
  for (int m=0; m<RG_numBCMessages; m++){
    // points are included: <=
    // remember one index is 0
    // check here actually

    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){
	  
	  if (ns >0)
	    rhons[0][ix_f + i][iy_f + j][iz_f + k]= Fx_BC[m][count];
	  if (ns >1)
	    rhons[1][ix_f + i][iy_f + j][iz_f + k]= Fy_BC[m][count];
	  if (ns >2)
	    rhons[2][ix_f + i][iy_f + j][iz_f + k]= Fz_BC[m][count];
	  if (ns >3)
	    rhons[3][ix_f + i][iy_f + j][iz_f + k]= Fa_BC[m][count];
	  
	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

void EMfields3D::setBC_Nodes_RENORM(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){
  for (int m=0; m<RG_numBCMessages; m++){
    // points are included: <=
    // remember one index is 0
    // check here actually

    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;
    
    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){
	  
	  //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize); 	  
	  Fx[ix_f + i][iy_f + j][iz_f + k]= Fx_BC[m][count]*(dx*dx);
	  Fy[ix_f + i][iy_f + j][iz_f + k]= Fy_BC[m][count]*(dx*dx);
	  Fz[ix_f + i][iy_f + j][iz_f + k]= Fz_BC[m][count]*(dx*dx);
	  
	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

// QUI
/* the RG sets the received BC                                                                                              
       Fx, Fy, Fz: the image
       vectX, vectY, vectZ: the intermediate solution in the iteration
       Fx_BC, Fy_BC, Fz_BC: the BCs              
       RGBC_Info: the RGBC_Info struct (Active or Ghost)                   
       RG_numBCMessages: number of messages (Active or Ghost) */
void EMfields3D::setBC_NodesImage(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double ***vectX, double ***vectY, double ***vectZ, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){

  for (int m=0; m<RG_numBCMessages; m++){


    // one direction has to be 0
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_NodesImage N0 !=2 (may be ok)" << endl;}

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  if (ix_f + i< 0 or ix_f + i>nxn-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }
	  if (iy_f + j< 0 or iy_f + j>nyn-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }
	  if (iz_f + k< 0 or iz_f + k>nzn-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }

	  // I am putting here Image = intermediate solution;
	  // in the source I will put the BC
	  Fx[ix_f + i][iy_f + j][iz_f + k]= vectX[ix_f + i][iy_f + j][iz_f + k] ;
	  Fy[ix_f + i][iy_f + j][iz_f + k]= vectY[ix_f + i][iy_f + j][iz_f + k] ;
	  Fz[ix_f + i][iy_f + j][iz_f + k]= vectZ[ix_f + i][iy_f + j][iz_f + k] ;
	  
	  count ++;


	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

/* this function is exactly the same as the previous setBC_NodesImage
 only differences: unnecessary inputs deleted & vector size included
 (this just for the safety check, to be deleted in later versions)*/
void EMfields3D::setBC_NodesImage(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double ***vectX, double ***vectY, double ***vectZ, RGBC_struct * RGBC_Info, int RG_numBCMessages, int nx, int ny, int nz){

  for (int m=0; m<RG_numBCMessages; m++){


    // one direction has to be 0
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_NodesImage N0 !=2 (may be ok)" << endl;}

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  if (ix_f + i< 0 or ix_f + i>nx-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }
	  if (iy_f + j< 0 or iy_f + j>ny-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }
	  if (iz_f + k< 0 or iz_f + k>nz-1){ cout <<"Problem is here, aborting..." << endl; MPI_Abort(MPI_COMM_WORLD, -1);  }

	  // I am putting here Image = intermediate solution;
	  // in the source I will put the BC
	  Fx[ix_f + i][iy_f + j][iz_f + k]= vectX[ix_f + i][iy_f + j][iz_f + k] ;
	  Fy[ix_f + i][iy_f + j][iz_f + k]= vectY[ix_f + i][iy_f + j][iz_f + k] ;
	  Fz[ix_f + i][iy_f + j][iz_f + k]= vectZ[ix_f + i][iy_f + j][iz_f + k] ;
	  
	  count ++;


	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

void EMfields3D::setBC_NodesImage_RENORM(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double ***vectX, double ***vectY, double ***vectZ, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){


  for (int m=0; m<RG_numBCMessages; m++){

    //cout << "R" << vct->getSystemWide_rank() << "I am putting BC image " << m <<"/ " << RG_numBCMessages <<" starting from [" << RGBC_Info[m].ix_first<<"-" << RGBC_Info[m].iy_first  <<"-" << RGBC_Info[m].iz_first<<"] for nodes [" << RGBC_Info[m].np_x <<"-" << RGBC_Info[m].np_y << "-" << RGBC_Info[m].np_z <<"]"<< endl;
    // one direction has to be 0
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_NodesImage N0 !=2 (may be ok)" << endl;}

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  // I am putting here Image = intermediate solution;
	  // in the source I will put the BC
	  Fx[ix_f + i][iy_f + j][iz_f + k]= vectX[ix_f + i][iy_f + j][iz_f + k]*(dx*dx) ;
	  Fy[ix_f + i][iy_f + j][iz_f + k]= vectY[ix_f + i][iy_f + j][iz_f + k]*(dx*dx) ;
	  Fz[ix_f + i][iy_f + j][iz_f + k]= vectZ[ix_f + i][iy_f + j][iz_f + k]*(dx*dx) ;

	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;
}

void EMfields3D::sendOneBC(VirtualTopology3D * vct, Grid * grid,  RGBC_struct CG_Info, int ch, int which){
  
  // 0: active
  // -1: ghost [but only for debug]
  int dest, tag;
  string TYPE;
  if (which== 0) TYPE= "ACTIVE";
  if (which== -1) TYPE= "GHOST";
  if (which== -10) TYPE= "II";
  if (which== 7) TYPE= "BUFFER";
  if (which== 9) TYPE= "FIX3B";

  MPI_Comm CommToCh;
  if (which==0 or which== -10)
    CommToCh=vct->getCommToChild(ch);
  else if (which== -1)
    CommToCh=vct->getCommToChild_BCGhost(ch);
  else if (which== 7)
    CommToCh=vct->getCommToChild_BCBuffer(ch);
  else if (which== 9)
    CommToCh=vct->getCommToChild_BCFix3B(ch);
  
  int rankAsParent;
  int rankAsParent_BCGhost;

  /*if (vct->getCommToChild(ch) != MPI_COMM_NULL and vct->getCommToChild_BCGhost(ch)== MPI_COMM_NULL){
    cout << "FATAL ERROR: duplicate failed  "<<endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // prova
  MPI_Comm_rank(vct->getCommToChild(ch), &rankAsParent ); 
  MPI_Comm_rank(vct->getCommToChild_BCGhost(ch), &rankAsParent_BCGhost);

  if (rankAsParent != rankAsParent_BCGhost){
    cout << "FATAL ERROR: ranks are not the same on communicator and copy  "<<endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
    }*/
  
  int HW=0;
  
  if (! ( CG_Info.np_x==1) ) { HW++; }
  if (! ( CG_Info.np_y==1) ) { HW++; }
  if (! ( CG_Info.np_z==1) ) { HW++; }

  int Size= CG_Info.np_x* CG_Info.np_y*CG_Info.np_z;

  /*if (HW != 2) {
    cout <<"R" << vct->getSystemWide_rank() << " WARNING: in sendOneBC, HW= " << HW << endl;
    }*/
      
  if (which == 0 or which == -1){
    buildBCMsg(vct, grid, ch, CG_Info, Size, CGMsg );
  }
  if (which == -10){
    buildIIMsg(vct, grid, ch, CG_Info, Size );
  }
  if (which == 7){
    buildBCMsg(vct, grid, ch, CG_Info, Size, CGMsgBuffer );
  }
  if (which == 9){
    buildFix3BMsg(vct, grid, ch, CG_Info, Size, CGMsgFix3B );
  }
  

  //cout << "R" << vct->getSystemWide_rank()  << " after building msg for ch " << ch  << endl;      

  dest= CG_Info.RG_core;  // on the parent-child communicator
  tag= CG_Info.MsgID;
    

  MPI_Request request;
  MPI_Status status;
  
  if (which == 0 or which == -1){
    MPI_Isend(CGMsg, Size*NumF, MPI_DOUBLE, dest, tag, CommToCh, &request);
    MPI_Wait(&request, &status);

  }
  if (which== -10){
    MPI_Isend(CGMsg_II, Size*NumF_II, MPI_DOUBLE, dest, tag, CommToCh, &request);
    MPI_Wait(&request, &status);

  }
  if (which== 7){
    MPI_Isend(CGMsgBuffer, Size*NumF, MPI_DOUBLE, dest, tag, CommToCh, &request);
    MPI_Wait(&request, &status);
  }
  if (which== 9){
    MPI_Isend(CGMsgFix3B, Size*Numfix3B, MPI_DOUBLE, dest, tag, CommToCh, &request);
    MPI_Wait(&request, &status);
  }
  

  //cout << "R" << vct->getSystemWide_rank() <<", R" << rankAsParent <<" on PC communicator, has sent this " << TYPE <<" msg to core " << dest << " in PC communicator, with size " << Size <<" and tag " << tag<< endl;

  return;

}

void EMfields3D::MPI_Barrier_ParentChild(VirtualTopology3D* vct){

  int RS= vct->getSystemWide_rank();
  int RR= vct->getCartesian_rank();
  int TOT= vct->getNprocs();

  /*if (RR==0){
    cout << "Grid " <<numGrid << " before MPI_Barrier_ParentChild- fields" << endl;
    }*/

  if (vct->getCommToParent() != MPI_COMM_NULL)  MPI_Barrier(vct->getCommToParent());

  for (int ch=0; ch< numChildren; ch++)
    if (vct->getCommToChild(ch) != MPI_COMM_NULL)
      MPI_Barrier(vct->getCommToChild(ch));

  /*if (RR==0){
    cout << "Grid " <<numGrid << " at the end of MPI_Barrier_ParentChild- fields" << endl;
    }*/

  if (RS==0)
    cout << "Everybody at the end of EMfields3D::MPI_Barrier_ParentChild" << endl;
}

void EMfields3D::initDoubleGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize DOUBLE GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
	  /*globalx= grid->getXN(i, j, k) + coarsedx + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + coarsedy + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + coarsedz + grid->getOz_SW();*/

	  globalx= grid->getXN(i, j, k) + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + grid->getOz_SW();
         
          double yB = globaly - .25 * Ly;
	  double yT = globaly - .75 * Ly;
          double yBd = yB / delta;
          double yTd = yT / delta;

	  double xB = globalx - .25 * Lx;
	  double xT = globalx - .75 * Lx;
	  double xBd = xB / delta;
	  double xTd = xT / delta;
	    
	  double xpert;
	  double ypert;
	    
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              double sech_yBd = 1. / cosh(yBd);
              double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) + tanh(-yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          
	  xpert = globalx- Lx/4;
	  ypert = globaly- Ly/4;
	  double deltax= Lx/2.0;
	  double deltay= Ly/2.0;
	  if (xpert < Lx/2 and ypert < Ly/2)
	    {
	      Bxn[i][j][k] +=(B0x*pertGEM)*(M_PI/deltay)*cos(2*M_PI*xpert/deltax)*sin(M_PI*ypert/deltay  );
	      Byn[i][j][k] = B0y -(B0x*pertGEM)*(2*M_PI/deltax)*sin(2*M_PI*xpert/deltax)*cos(M_PI*ypert/deltay);
	    }
	  // add the second initial GEM perturbation                                                                     
	  xpert = globalx- 3*Lx/4;
	  ypert = globaly- 3*Ly/4;
	  if (xpert > Lx/2 and ypert > Ly/2)
	    {
	      Bxn[i][j][k] +=(B0x*pertGEM)*(M_PI/deltay)*cos(2*M_PI*xpert/deltax)*sin(M_PI*ypert/deltay  );
	      Byn[i][j][k] = B0y -(B0x*pertGEM)*(2*M_PI/deltax)*sin(2*M_PI*xpert/deltax)*cos(M_PI*ypert/deltay);
	    }
	  double exp_pert;
	  // add the initial X perturbation                                                                              
	  xpert = globalx- Lx/4;
	  ypert = globaly- Ly/4;
	  exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

	  Bxn[i][j][k] +=(B0x*pertX)*exp_pert*(
					       -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
					       -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
					       );

	  Byn[i][j][k] +=(B0x*pertX)*exp_pert*(
					       cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
					       +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
					       );
	  // add the second initial X perturbation                                                                       
	  xpert = globalx- 3*Lx/4;
	  ypert = globaly- 3*Ly/4;
	  exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

	  Bxn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						-cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
						-cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
						);

	  Byn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
						+sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
						);
	  // guide field                                                                                                 
	  Bzn[i][j][k] = B0z;
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

    // Lambda
    
    /*for (int i=0; i < nxn; i++)                   
      for (int j=0; j < nyn; j++)                               
	for (int k=0; k < nzn; k++){                          

	  double yC1=  (grid->getOy_SW() +grid->getYN(i, j, k) - 1./4.*Ly)/ (10*delta); //Lambda[i][j][k]=2.0 * M_PI / dy* fabs(tanh(yC));
	  double yC2=  (grid->getOy_SW() +grid->getYN(i, j, k) - 3./4.*Ly)/ (10*delta);
	  Lambda[i][j][k]=2.0 * M_PI / dy* (-1 + fabs(tanh(yC1)) + fabs(tanh(yC2)) );
	    
	} 
    */
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }

}


// *** //
void EMfields3D::initMAX_Show_RG_BC(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X

  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize MAX_Show_RG_BC" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
         }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
	    
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0y;
          Byn[i][j][k] = B0y;
	  Bzn[i][j][k] = B0z;
          
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }

}

// *** //

void EMfields3D::initDoubleGEM_CentralPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize DOUBLE GEM Challenge with Central Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
	  globalx= grid->getXN(i, j, k) + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + grid->getOz_SW();

          double yB = globaly - .25 * Ly;
	  double yT = globaly - .75 * Ly;
          double yBd = yB / delta;
          double yTd = yT / delta;

	  double xB = globalx - .25 * Lx;
	  double xT = globalx - .75 * Lx;
	  double xBd = xB / delta;
	  double xTd = xT / delta;
	    
	  double xpert;
	  double ypert;
	    
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              double sech_yBd = 1. / cosh(yBd);
              double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] +=rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) + tanh(-yTd));
	  //cout << "Bxn[i][j][k] "<< Bxn[i][j][k] <<endl;
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          
	  // add the second initial GEM perturbation                                                                     
	  xpert = globalx- Lx/2;
	  ypert = globaly- 3*Ly/4;
	  if (xpert > Lx/4 and xpert <3*Lx/4 and ypert > Ly/2)
	    {
	      Bxn[i][j][k] +=(B0x*pertGEM)*(M_PI/deltay)*cos(2*M_PI*xpert/deltax)*sin(M_PI*ypert/deltay  );
	      Byn[i][j][k] = B0y -(B0x*pertGEM)*(2*M_PI/deltax)*sin(2*M_PI*xpert/deltax)*cos(M_PI*ypert/deltay);
	    }
	  double exp_pert;

	  // add the second initial X perturbation                                                                       
	  xpert = globalx- Lx/2;
	  ypert = globaly- 3*Ly/4;
	  exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

	  Bxn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						-cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
						-cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
						);

	  Byn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
						+sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
						);
	  // guide field                                                                                                 
	  Bzn[i][j][k] = B0z;
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

    // Lambda
    
    /*for (int i=0; i < nxn; i++)                   
      for (int j=0; j < nyn; j++)                               
	for (int k=0; k < nzn; k++){                          

	  double yC1=  (grid->getOy_SW() +grid->getYN(i, j, k) - 1./4.*Ly)/ (10*delta); //Lambda[i][j][k]=2.0 * M_PI / dy* fabs(tanh(yC));
	  double yC2=  (grid->getOy_SW() +grid->getYN(i, j, k) - 3./4.*Ly)/ (10*delta);
	  Lambda[i][j][k]=2.0 * M_PI / dy* (-1 + fabs(tanh(yC1)) + fabs(tanh(yC2)) );
	    
	} 
    */
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initWeightProj(Grid *grid, VirtualTopology3D *vct){

  // copy as much as possible form initWeightBC

  int rank_local=vct->getCartesian_rank();
  MPI_Comm CommToParent= vct->getCommToParent(); // if != MPI_COMM_NULL the grid is a child
 
  size_RG_ProjMsg=0;
  size_CG_ProjMsg=0;
  RG_numProjMessages=0;

  /* here, vectors where core 0 of the local child grid assembles the messages
      from all the local grid 
      de-allocated at the end of the function */
  RGBC_struct * RGProj_Info_LevelWide;
  int RG_numProjMessage_LevelWide=0;

  /* phase 1, as a child */
  if (CommToParent!= MPI_COMM_NULL){

    RGProj_Info= new RGBC_struct[MAX_RG_numBCMessages];

    initWeightProj_Phase1(grid, vct);

    // *3 for the three componenets of E
    RG_ProjMsg= new double[(size_RG_ProjMsg+1)*3];

    // here, in the BC versions there are several checks which i did NOT implement
  }

  /* end phase 1 */

  /* phase 2, RG sends info to CG */
  // children grids
  if (CommToParent != MPI_COMM_NULL) {

    // phase 2a: children assemble all the info in core 0 in the LOCAL child grid

    if (rank_local >0){
    /* send one message more; the last message has -1 in the RG_core
       to signal end of 'valid' messages */

      MPI_Send(RGProj_Info, RG_numProjMessages+1, MPI_RGBC_struct, 0, TAG_PROJ, vct->getComm());
    }// end if (rank_local >0){
    if (rank_local ==0){

      // only rank_local==0 needs to instantiate this
      RGProj_Info_LevelWide = new RGBC_struct[MAX_size_LevelWide];
      
      // recycling initWeightBC_Phase2a, maybe change name 
      initWeightBC_Phase2a(grid, vct, RGProj_Info_LevelWide, &RG_numProjMessage_LevelWide, RGProj_Info, RG_numProjMessages, 3);

    } // end if (rank_local ==0){

    
    // phase 2b: core 0 of the child grid assembles & sends messages for all CG cores
    if (rank_local==0){
      // recycling initWeightBC_Phase2b, maybe change name  
      initWeightBC_Phase2b(grid, vct, RGProj_Info_LevelWide, RG_numProjMessage_LevelWide, 3);

    } // end if (rank_local==0), phase 2b

  } // end if (CommToParent != MPI_COMM_NULL) {
  
  // phase 2c, only for the parent
  // each CG core will receive a message per child grid
  if (numChildren > 0 ){

    size_CG_ProjMsg=0;
    Max_CG_numProjMessages=0;

    CGProj_Info= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
    CG_numProjMessages = new int [numChildren];

    for (int i=0; i< numChildren; i++){
      CG_numProjMessages[i]=0;
    }

    // recycling initWeightBC_Phase2c, maybe change name 
    initWeightBC_Phase2c(grid, vct, CGProj_Info, CG_numProjMessages, 3);

    // do i need to instantiate stuff??
    /*CGProjVectors_Needed= false;

    for (int ch=0; ch<numChildren; ch++){
      if (CG_numProjMessages[ch]>0) {CGProjVectors_Needed= true; break;}
      }*/
    
    if (true){

      // NOW, instantiate vector with the weights and fill them
      
      double inv_dx= 1./dx;
      double inv_dy= 1./dy;
      double inv_dz= 1./dz;
      
      // the vector where you will actually receive the projection message
      CG_ProjMsg= new double[3*size_CG_ProjMsg +3];

      ProjWeight000= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight001= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight010= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight011= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight100= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight101= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight110= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjWeight111= newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      
      ProjIX=newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjIY=newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      ProjIZ=newArr3(double, numChildren, Max_CG_numProjMessages, size_CG_ProjMsg);
      
      for (int ch=0; ch< numChildren; ch++){
	for (int m=0; m< CG_numProjMessages[ch]; m++){
	  RGBC_struct msg= CGProj_Info[ch][m];
	  
	  double x0= msg.CG_x_first;
	  double y0= msg.CG_y_first;
	  double z0= msg.CG_z_first;
	  
	  double xp, yp, zp;
	  int count=0;
	  
	  for (int i= 0; i<msg.np_x; i++){
	    for (int j=0; j<msg.np_y; j++){
	      for (int k=0; k<msg.np_z; k++){
		
		// ok, now each RG point is treated as a particle here
		// and i need the E at the point position
		
		xp= x0 + i*dx_Ch[ch];
		yp= y0 + j*dy_Ch[ch];
		zp= z0 + k*dz_Ch[ch]; 
		
		// this is copied from the mover
		const double ixd = floor((xp - xStart) * inv_dx);
		const double iyd = floor((yp - yStart) * inv_dy);
		const double izd = floor((zp - zStart) * inv_dz);
		int ix = 2 + int (ixd);
		int iy = 2 + int (iyd);
		int iz = 2 + int (izd);

		if (ix <1 or iy< 1 or iz<1){
		  cout << "Grid " << numGrid << " R " << vct->getCartesian_rank() << "RN position on CG: " << xp << " " <<yp <<" " << zp << " xstart " << xStart <<" " << yStart <<" " << zStart << " Ix" << ix << " " <<iy <<" " <<iz <<endl;
		  cout << "Aborting "<< endl;
		  MPI_Abort(MPI_COMM_WORLD, -1);
		  
		}

		if (ix < 1)
		  ix = 1;
		if (iy < 1)
		  iy = 1;
		if (iz < 1)
		  iz = 1;
		if (ix > nxn - 1)
		  ix = nxn - 1;
		if (iy > nyn - 1)
		  iy = nyn - 1;
		if (iz > nzn - 1)
		  iz = nzn - 1;
		
		double xi  [2];
		double eta [2];
		double zeta[2];
		
		xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
		eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
		zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
		xi  [1] = grid->getXN(ix,iy,iz) - xp;
		eta [1] = grid->getYN(ix,iy,iz) - yp;
		zeta[1] = grid->getZN(ix,iy,iz) - zp;
		
		ProjWeight000[ch][m][count] = xi[0] * eta[0] * zeta[0] * invVOL;
		ProjWeight001[ch][m][count] = xi[0] * eta[0] * zeta[1] * invVOL;
		ProjWeight010[ch][m][count] = xi[0] * eta[1] * zeta[0] * invVOL;
		ProjWeight011[ch][m][count] = xi[0] * eta[1] * zeta[1] * invVOL;
		ProjWeight100[ch][m][count] = xi[1] * eta[0] * zeta[0] * invVOL;
		ProjWeight101[ch][m][count] = xi[1] * eta[0] * zeta[1] * invVOL;
		ProjWeight110[ch][m][count] = xi[1] * eta[1] * zeta[0] * invVOL;
		ProjWeight111[ch][m][count] = xi[1] * eta[1] * zeta[1] * invVOL;
		
		ProjIX[ch][m][count]= ix;
		ProjIY[ch][m][count]= iy;
		ProjIZ[ch][m][count]= iz;
		
		//cout << "ProjWeight000[ch][m][count]: " << ProjWeight000[ch][m][count] << " xi[0] " << xi[0] << " xi[1] " << xi[1]<<" eta[0] " << eta[0] << " zeta[0] " << zeta[0] <<" invVOL" << invVOL << " xp " << xp <<" yp " << yp << " zp " << zp << endl; 
		count++;
	      }
	    }
	  }
	} // end for (int m=0; m< CG_numProjMessages[m]; m++){
      } // end for(int i=0; i< numChildren; i++){

      // now instantantiate the (very wasteful) vectors for projection
      ExthProjSt= newArr4(double, numChildren, nxn, nyn, nzn);
      EythProjSt= newArr4(double, numChildren, nxn, nyn, nzn);
      EzthProjSt= newArr4(double, numChildren, nxn, nyn, nzn);
      DenProjSt= newArr4(double, numChildren, nxn, nyn, nzn);

      // DenProjSt can be filled right now

      eqValue(0.0, DenProjSt, numChildren, nxn, nyn, nzn);

      for (int ch=0; ch< numChildren; ch++){
	for (int m=0; m< CG_numProjMessages[ch]; m++){
	  int CC= CGProj_Info[ch][m].np_x* CGProj_Info[ch][m].np_y * CGProj_Info[ch][m].np_z;
	  for (int c=0; c<CC; c++){
	    
	    int IX= ProjIX[ch][m][c];
	    int IY= ProjIY[ch][m][c];
	    int IZ= ProjIZ[ch][m][c];
	    
	    DenProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][m][c] ;
            DenProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][m][c] ; 
            DenProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][m][c] ;
            DenProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][m][c]; 
            DenProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][m][c] ;
            DenProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][m][c]; 
            DenProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][m][c] ;
            DenProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][m][c]; 
	    
	  } // int CC= CGProj_Info[ch][i].np_x* CGProj_Info[ch][i].np_y * CGProj_Info[ch][i].np_z;
	} // end  for (int i=0; i< CG_numProjMessages[ch]; i++){

	communicateNode_Proj(nxn, nyn, nzn, DenProjSt[ch], 0, 0, 0, 0, 0, 0, vct) ;
      } // end for (int ch=0; ch< numChildren; ch++){
      
    } // end if(Needed) --> became end if (true)
  } // end  if (numChildren > 0 ){
  
  Trim_Proj_Vectors(vct); 

  // here, to decide wether you have to enter ApplyProjection
  ApplyProjection= false;
  for (int ch=0; ch< numChildren; ch++)
    for (int i=0; i<nxn; i++)
      for (int j=0; j<nyn; j++)
	for (int k=0; k<nzn; k++){
	  // I am checking DenProjSt[ch][i][j][k] because you may have a != 0 after the communicateNode_Proj
	  // not necessarily only after the messaging thing
 	  if (DenProjSt[ch][i][j][k] >0.0){
	    ApplyProjection= true;
	    break;
	  }
	}
  if (ApplyProjection){
    Ex_n = newArr3(double, nxn, nyn, nzn);
    Ey_n = newArr3(double, nxn, nyn, nzn);
    Ez_n = newArr3(double, nxn, nyn, nzn);
  }
  // end here, to decide wether you have to enter ApplyProjection

  if (CommToParent != MPI_COMM_NULL && rank_local==0){
    delete[] RGProj_Info_LevelWide;
  }

  /*
  int RL= vct->getCartesian_rank();
  MPI_Barrier(MPI_COMM_WORLD);
  if (numChildren >0){
    for (int ch=0; ch < numChildren; ch++){
      int CG_PC= vct->getRank_CommToChildren(ch);

      int waste=0;
      for (int ii=0; ii< 10000* RL; ii++){
	waste= waste+1;
      }
      cout << waste << endl;
      cout << "CG core " << CG_PC << " expects " << CG_numProjMessages[ch] << " projections messages:" << endl;
      for (int j=0; j <CG_numProjMessages[ch]; j++){
	RGBC_struct msg = CGProj_Info[ch][j];
	if (CG_PC != msg.CG_core) {cout << "FATAL ERROR IN RANK, ABORTING" << endl; MPI_Abort(MPI_COMM_WORLD, -1);}
	cout << "CG core " << msg.CG_core << " RG core " << msg.RG_core << " Id " << msg.MsgID << " count "  << msg.np_x*msg.np_y*msg.np_z<< endl; 

	
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (CommToParent != MPI_COMM_NULL){
    int RG_PC= vct->getRank_CommToParent();
    
    int waste=0;
    for (int ii=0; ii< 100000000* RL; ii++){
      waste= waste+1;
    }
    cout << waste << endl;

    cout << "RG core " << RG_PC << " sends " << RG_numProjMessages << " projections messages:" << endl;
    for (int j=0; j< RG_numProjMessages; j++){
      RGBC_struct msg = RGProj_Info[j];
      if (RG_PC != msg.RG_core) {cout << "FATAL ERROR IN RANK, RG, ABORTING" << endl; MPI_Abort(MPI_COMM_WORLD, -1);}
      cout << "RG core " << msg.RG_core << " CG core " << msg.CG_core << " Id " << msg.MsgID << " count "  << msg.np_x*msg.np_y*msg.np_z << endl;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Abort(MPI_COMM_WORLD, -1);

  MPI_Barrier(MPI_COMM_WORLD);*/

  MPI_Barrier(vct->getComm());
  if (rank_local==0){
    cout << "I am grid " << numGrid <<", I am after the barrier at the end of initWeightProj" << endl;
  }
  
}



void EMfields3D::initWeightProj_Phase1(Grid *grid, VirtualTopology3D *vct){

  // RG cores build the map for proj info

  // i do not want to send back the BC
  int i_s=2, i_e= nxn-1-2;
  int j_s=2, j_e= nyn-1-2;
  int k_s=2, k_e= nzn-1-2;
  
  // the right core sends the shared node
  if (vct->getXleft_neighbor()!= MPI_PROC_NULL and vct->getXLEN()>1) i_s=1; // so i don't send the shared node twice
  if (vct->getYleft_neighbor()!= MPI_PROC_NULL and vct->getYLEN()>1) j_s=1; // so i don't send the shared node twice  
  if (vct->getZleft_neighbor()!= MPI_PROC_NULL and vct->getZLEN()>1) k_s=1; // so i don't send the shared node twice  

  
  char nn= 'p';
  Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RGProj_Info, &RG_numProjMessages, &size_RG_ProjMsg, vct, nn);

}

void EMfields3D::Explore3DAndCommit(Grid *grid, int i_s, int i_e, int j_s, int j_e, int k_s, int k_e, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg, VirtualTopology3D * vct , char  dir){
  // here, the same as Explore3DAndCommit for particles
  
  // policy:
  // explore Z dir
  // for on the number of cores found there: explore Y dir
  // for on the number of cores found there: explore X dir
  // finally, commit  NB: all faces should have the same c

  int MS= nxn; if (nyn>MS) MS= nyn; if (nzn>MS) MS= nzn;
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
  grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

  for (int n=0; n<Dir1_Ncores; n++){ // it will find again the core in Dir 1, but it will also explore Dir 2 
    // Y dir / Dir 2
    grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n],  Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); 
    
    for (int m=0; m< Dir2_Ncores; m++){ //it will find again the core in Dir 1, but it will also explore Dir 2 
      // X dir / Dir 3
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e, Dir2_IndexFirstPointperC[m],  Dir1_IndexFirstPointperC[n], Dir3_SPXperC, Dir3_SPYperC, Dir3_SPZperC, Dir3_NPperC, Dir3_rank, &Dir3_Ncores, Dir3_IndexFirstPointperC);

      for (int NN=0; NN< Dir3_Ncores; NN++){
	// using the function written for BCs, it's the same

	//void Assign_RGBC_struct_Values(RGBC_struct *s, int ix_first_tmp, int iy_first_tmp, int iz_first_tmp, int BCside_tmp, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp, int ID);

	Assign_RGBC_struct_Values(RGBC_Info + (*numMsg), Dir3_IndexFirstPointperC[NN], Dir2_IndexFirstPointperC[m], Dir1_IndexFirstPointperC[n], -1, Dir3_NPperC[NN], Dir2_NPperC[m], Dir1_NPperC[n], Dir3_SPXperC[NN], Dir3_SPYperC[NN], Dir3_SPZperC[NN], Dir3_rank[NN], rank_CTP, *numMsg);
	
	
	if (dir != 'p'){
	  DirBuffer[*numMsg]= dir;
	}

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


void EMfields3D::Trim_Proj_Vectors(VirtualTopology3D *vct){
  int dim;
  
  MPI_Comm CommToParent= vct->getCommToParent();

  if (CommToParent != MPI_COMM_NULL){ // then trim RGProj_Info
    
    dim= RG_numProjMessages+1;
    RGBC_struct *tmp= new RGBC_struct[dim];
    
    for (int i=0; i< dim; i++)
      tmp[i]= RGProj_Info[i];

    delete[]RGProj_Info;
    
    RGProj_Info= new RGBC_struct[dim];

    for (int i=0; i< dim; i++)
      RGProj_Info[i]= tmp[i];
    
    delete[]tmp;
  }
  
  if (numChildren>0){
    dim=0;
    for (int ch=0; ch < numChildren; ch++){
      if (CG_numProjMessages[ch]> dim) {dim= CG_numProjMessages[ch];}
    }

    dim= dim+1;

    RGBC_struct **tmp2= newArr2(RGBC_struct, numChildren, dim);
    
    for (int ch=0; ch < numChildren; ch++)
      for (int j=0; j< CG_numProjMessages[ch]; j++)
	tmp2[ch][j]= CGProj_Info[ch][j];
    
    delArr2(CGProj_Info, numChildren);

    CGProj_Info= newArr2(RGBC_struct, numChildren, dim);
    for (int ch=0; ch < numChildren; ch++)
      for (int j=0; j< CG_numProjMessages[ch]; j++)
	CGProj_Info[ch][j]= tmp2[ch][j];
    
    delArr2(tmp2, numChildren);
  } // end if (numChildren>0){

}

void EMfields3D::sendProjection(Grid *grid, VirtualTopology3D *vct){
  
  if (vct->getCommToParent()==MPI_COMM_NULL) return; // only children do this

  cout << "numGrid " << numGrid << " R " << vct->getCartesian_rank() << " is sending proj" << endl;
  
  int dest;
  int tag;
  MPI_Request request;
  MPI_Status status;

  for (int m=0; m< RG_numProjMessages; m++){
    
    // build and send each RG Proj msg
    RGBC_struct msg= RGProj_Info[m];
    int count=0;

    int Size= msg.np_x*msg.np_y*msg.np_z;

    for (int i= 0; i< msg.np_x; i++)
      for (int j= 0; j< msg.np_y; j++)
	for (int k= 0; k< msg.np_z; k++){
	  int xC= msg.ix_first+i;
	  int yC= msg.iy_first+j;
	  int zC= msg.iz_first+k;
	  
	  RG_ProjMsg[0*Size + count]=Exth[xC][yC][zC];
	  RG_ProjMsg[1*Size + count]=Eyth[xC][yC][zC];
	  RG_ProjMsg[2*Size + count]=Ezth[xC][yC][zC];

	  count ++;
	}
      
    dest= msg.CG_core;
    tag= msg.MsgID;

    /* i did not even try with MPI_Send, just went directly for the non-blocking
       to let RG go immediately after sending*/
    MPI_Isend(RG_ProjMsg, Size*3, MPI_DOUBLE, dest, tag, vct->getCommToParent_Proj(), &request);
    MPI_Wait(&request, &status);
    

  } // end for (int m=0; m< RG_numProjMessages; m++){

}

void EMfields3D::receiveProjection(Grid *grid, VirtualTopology3D *vct){
  
  int RR= vct->getCartesian_rank();

  if (numChildren <1) return; // this is only for parents

  //if (CGProjVectors_Needed== false) return; // if you are not involved in proj operators, exit

  cout << "Grid " << numGrid << " R " << RR << " has started receiveProjection" <<endl;

  /* set staging vectors (ExthProjSt, EythProjSt, EzthProjSt, DenProjSt) to 0
     before starting accumulatimg */
  eqValue(0.0, ExthProjSt, numChildren, nxn, nyn, nzn);
  eqValue(0.0, EythProjSt, numChildren, nxn, nyn, nzn);
  eqValue(0.0, EzthProjSt, numChildren, nxn, nyn, nzn);


  /* I have to receive child by child because different children send on different communicators
     this will introduce delays;
     should not produce DD */
  MPI_Status status;
  int count;
  bool found;
  int countExp;
  bool Testing= true;

  for (int ch=0; ch< numChildren; ch++){ // check projection from all children

    int R_PC= vct->getRank_CommToChildren(ch);
    for (int m=0; m< CG_numProjMessages[ch]; m++){ // check all the messages to this core
      
      MPI_Recv(CG_ProjMsg, (size_CG_ProjMsg+1)*3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, vct->getCommToChild_Proj(ch), &status);
      MPI_Get_count(&status, MPI_DOUBLE, &count );
    
      // store them immediately after receiving them, hoping it will minimize waits

      // 1. try to match the msg
      found= false;
      
      for (int i=0; i< CG_numProjMessages[ch]; i++){ // for to match the msg
	// I have to match BOTH id and RG_Core, because the id's are unique to the RG cores, not to the CG cores
	
	if (CGProj_Info[ch][i].MsgID == status.MPI_TAG and CGProj_Info[ch][i].RG_core == status.MPI_SOURCE){
	  found = true;
	  
	  countExp= CGProj_Info[ch][i].np_x* CGProj_Info[ch][i].np_y* CGProj_Info[ch][i].np_z;

	  //cout << "Grid "<< numGrid << " R " <<vct->getRank_CommToChildren(ch)<< " found msg: count: " << countExp << " S " << status.MPI_SOURCE << " ID " << status.MPI_TAG << endl;

	  if   (Testing){
	    if ( ! (countExp*3 == count)){
	      cout << "Grid " << numGrid << " R " << RR << ": fatal error in receiveProjection, I expect a msg with size " << countExp*3 << ", I got " << count ;
	      cout << "Aborting ..." << endl;
	      MPI_Abort(MPI_COMM_WORLD, -1);
	    }
	  } // if (Testing){ 

	  // STORE IT
	  int CC= CGProj_Info[ch][i].np_x* CGProj_Info[ch][i].np_y * CGProj_Info[ch][i].np_z;
	  for (int c=0; c<CC; c++){
	    
	    int IX= ProjIX[ch][i][c];
	    int IY= ProjIY[ch][i][c];
	    int IZ= ProjIZ[ch][i][c];

	    //cout << "Grid " << numGrid << " R " << vct->getCartesian_rank() << "IX " <<IX << " IY " << IY <<" IZ " <<IZ <<endl;
	    // NB: the denominator is written only once, in initWeightProj
	    
	    ExthProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * CG_ProjMsg[0*CC+c];
	    ExthProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * CG_ProjMsg[0*CC+c];

	    EythProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * CG_ProjMsg[1*CC+c];
	    EythProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * CG_ProjMsg[1*CC+c];

	    EzthProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * CG_ProjMsg[2*CC+c];
	    EzthProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * CG_ProjMsg[2*CC+c];

	    /* this is to test the projection operator, see instruction in TestProjection 

	    ExthProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * 1;
	    ExthProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * 1;
	    ExthProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * 1;
	    ExthProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * 1;
	    ExthProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * 1;
	    ExthProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * 1;
	    ExthProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * 1;
	    ExthProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * 1;

	    EythProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * 1;
	    EythProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * 1;
	    EythProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * 1;
	    EythProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * 1;
	    EythProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * 1;
	    EythProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * 1;
	    EythProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * 1;
	    EythProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * 1;

	    EzthProjSt[ch][IX][IY][IZ] += ProjWeight000[ch][i][c] * 1;
	    EzthProjSt[ch][IX][IY][IZ-1] += ProjWeight001[ch][i][c] * 1;
	    EzthProjSt[ch][IX][IY-1][IZ] += ProjWeight010[ch][i][c] * 1;
	    EzthProjSt[ch][IX][IY-1][IZ-1] += ProjWeight011[ch][i][c] * 1;
	    EzthProjSt[ch][IX-1][IY][IZ] += ProjWeight100[ch][i][c] * 1;
	    EzthProjSt[ch][IX-1][IY][IZ-1] += ProjWeight101[ch][i][c] * 1;
	    EzthProjSt[ch][IX-1][IY-1][IZ] += ProjWeight110[ch][i][c] * 1;
	    EzthProjSt[ch][IX-1][IY-1][IZ-1] += ProjWeight111[ch][i][c] * 1;*/
	    
	    //cout << "Grid " << numGrid << " R " << vct->getCartesian_rank() << " I have just finished dealing with a msg " << ExthProjSt[ch][IX][IY][IZ] << " " << ProjWeight000[ch][i][c] <<" " <<CG_ProjMsg[0*CC+c] << endl;
	  }
	  
	  break;
	}   //  if (CGProj_Info[ch][i].MsgID == status.MPI_TAG and ...
      
	
      } // end for (int i=0; i< CG_numProjMessages[ch]; i++){ // for to match it
      
      if (found== false){ // I could not match the msg
	  cout << "Grid " << numGrid << " R (PC rank) " << R_PC<< ": fatal error in receiveProjection, I could not match a msg..." ;
	  cout << "Msg is from core (PC comm) " << status.MPI_SOURCE << " size:  " << count << " ID:  " << status.MPI_TAG << " size " << count << endl;
	  cout << "Aborting ..." << endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);
      } // end if (found== false){ 
      
    } // end for (int m=0; m< CG_numProjMessages[m]; m++){ // check all the messages to this core 
  

  } // end for (int ch=0; ch< numChildren; ch++){ // check projection from all children 

  for (int ch=0; ch < numChildren; ch++){
    // this is if the projection ends up in shared nodes
    communicateNode_Proj(nxn, nyn, nzn, ExthProjSt[ch], 0, 0, 0, 0, 0, 0, vct) ;
    communicateNode_Proj(nxn, nyn, nzn, EythProjSt[ch], 0, 0, 0, 0, 0, 0, vct) ;
    communicateNode_Proj(nxn, nyn, nzn, EzthProjSt[ch], 0, 0, 0, 0, 0, 0, vct) ;
  }
  //cout << "Grid " << numGrid <<" R " << " has finished receiveProjection" <<endl;

}

void EMfields3D::TestProjection(Grid *grid, VirtualTopology3D *vct){

  /* instructions for the test:
     1. this method has to be positioned here:
        if (MLMD_PROJECTION ){
        EMf->receiveProjection(grid,vct);
        EMf->TestProjection(grid, vct);
       }
     2. in receiveProjection, do not sum 'received value'* weight, but 1*weight
     3. in the inputfile, enable also BC and particle repopulation, otherwise RG messages
        that should be received at different cycles may be erroneously received at the same cycke
     4. Ex =1 says that you received everything and you are appropriately dividing
        Ey[i][j][k]= ExthProjSt[ch][i][j][k]
	Ez[i][j][k]= DenProjSt [ch][i][j][k] counts the contribution for each grid point
     5. to avoid mess, disable mover and solver
     GOOD LUCK!
  */
     

  if (numChildren <1) return; // this is only for parents   

  // I have to set E to 0 otherwise it will keep summing
  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++){
	Ex[i][j][k]=0;
	Ey[i][j][k]=0;
	Ez[i][j][k]=0;
      }


  if (true){

    for (int ch=0; ch< numChildren; ch++){
      for (int i=0; i<nxn; i++)
	for (int j=0; j<nyn; j++)
	  for (int k=0; k<nzn; k++){
	    if (fabs(DenProjSt [ch][i][j][k]) >0){
	      
	      Ex[i][j][k]= ExthProjSt[ch][i][j][k] /DenProjSt [ch][i][j][k];
	      Ey[i][j][k]= ExthProjSt[ch][i][j][k] ;
	      Ez[i][j][k]= DenProjSt [ch][i][j][k];
	      
	    }
	  }
    }
  }

}

void EMfields3D::initTestProjection(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize initTestProjection " << endl;
      cout << "------------------------------------------" << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
	  /*globalx= grid->getXN(i, j, k) + coarsedx + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + coarsedy + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + coarsedz + grid->getOz_SW();*/
         
	  globalx= grid->getXN(i, j, k) + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + grid->getOz_SW();

          // initialize the density for species
          for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
	  if (numGrid ==0){
	    Ex[i][j][k] = 0.0;
	    Ey[i][j][k] = 0.0;
	    Ez[i][j][k] = 0.0;}
	  else{
	    Ex[i][j][k] = 1*globalx;
	    Ey[i][j][k] = 2*globaly;
	    Ez[i][j][k] = 3*globalz;

	    Exth[i][j][k] = 1*globalx;
	    Eyth[i][j][k] = 2*globaly;
	    Ezth[i][j][k] = 3*globalz;
	    //cout << "Ex[i][j][k] " << Ex[i][j][k]  << "Ey[i][j][k] " << Ey[i][j][k]  << "Ez[i][j][k] " << Ez[i][j][k] << endl;
	  }
          // Magnetic field
          Bxn[i][j][k] = B0x;
          Byn[i][j][k] = B0y;
	  Byn[i][j][k] = B0z;
          // add the initial X perturbation
          
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


void EMfields3D::initTestBC(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize initTestBC " << endl;
      cout << "------------------------------------------" << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
	  /*globalx= grid->getXN(i, j, k) + coarsedx + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + coarsedy + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + coarsedz + grid->getOz_SW();*/

	  globalx= grid->getXN(i, j, k) +grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) +grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) +grid->getOz_SW();
         
	    
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
	  if (numGrid ==0){
	    Ex[i][j][k] = globalx;
	    Ey[i][j][k] = 10+ globaly;
	    Ez[i][j][k] = 20 +globalz;

	    Bxn[i][j][k] = 30 +globalx;
	    Byn[i][j][k] = 40 +globaly;
	    Bzn[i][j][k] = 50 +globalz;

	  }
	  else{
	    Ex[i][j][k] = 0.0;
	    Ey[i][j][k] = 0.0;
	    Ez[i][j][k] = 0.0;
	    
	    Bxn[i][j][k] = 0.0;
	    Byn[i][j][k] = 0.0;
	    Byn[i][j][k] = 0.0;

	    //cout << "Ex[i][j][k] " << Ex[i][j][k]  << "Ey[i][j][k] " << Ey[i][j][k]  << "Ez[i][j][k] " << Ez[i][j][k] << endl;
	  }
          
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

/** fix3B **/
void EMfields3D::initTestFix3B(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize initTestFix3B " << endl;
      cout << "------------------------------------------" << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {

	  globalx= grid->getXN(i, j, k) +grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) +grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) +grid->getOz_SW();
         
	    
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
	    rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	   

	  if (numGrid ==0){
	    Bxn[i][j][k] = 10 +globalx;
	    Byn[i][j][k] = 100 +globaly;
	    Bzn[i][j][k] = 1000 +globalz;
	  }
	  else{
	    Bxn[i][j][k] = 0.0;
	    Byn[i][j][k] = 0.0;
	    Byn[i][j][k] = 0.0;

	    //cout << "Ex[i][j][k] " << Ex[i][j][k]  << "Ey[i][j][k] " << Ey[i][j][k]  << "Ez[i][j][k] " << Ez[i][j][k] << endl;
	  }
          
        }
    // communicate ghost
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);

  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}
/** end fix3B **/

void EMfields3D::applyProjection(Grid *grid, VirtualTopology3D *vct, Collective * col){
  
  // NB: at this stage, I assume the refined grids at the same levels do not overlap

  if (numChildren <1) return; // this is only for parents

  if (ApplyProjection == true ){
  
    cout << "Applying projection" << endl;

    for (int ch=0; ch < numChildren; ch++){

      double RFx= dx/ dx_Ch[ch];
      double RFy= dy/ dy_Ch[ch];
      double RFz= dz/ dz_Ch[ch];
      
      for (int i=0; i< nxn; i++)
	for (int j=0; j< nyn; j++)
	  for (int k=0; k< nzn; k++){
	    //if (DenProjSt[ch][i][j][k]>0){
	    // to exclude boundary nodes
	    if (DenProjSt[ch][i][j][k]> RFx*RFy*RFz*0.99){
	      // average 
	      cout << "Exth[i][j][k] before: " << Exth[i][j][k] << endl;
	      Exth[i][j][k]= (Exth[i][j][k] + ExthProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k])/2;
	      Eyth[i][j][k]= (Eyth[i][j][k] + EythProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k])/2;
	      Ezth[i][j][k]= (Ezth[i][j][k] + EzthProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k])/2;
	      cout << "Exth[i][j][k] after: " << Exth[i][j][k] << " ExthProjSt[ch][i][j][k] " << ExthProjSt[ch][i][j][k] << " DenProjSt[ch][i][j][k] " << DenProjSt[ch][i][j][k] << endl;
	      // substitution
	      /*Exth[i][j][k]= ExthProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k];
		Eyth[i][j][k]= EythProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k];
		Ezth[i][j][k]= EzthProjSt[ch][i][j][k] / DenProjSt[ch][i][j][k];*/
	      
	    } // end if (DenProjSt[ch][i][j][k]>0){
	  } // end for (int k=0; k< nzn; k++)
    } // end cycle on children
    // now regenerate E n+1

    addscale(1 / th, -(1.0 - th) / th, Ex, Ex_n, Exth, nxn, nyn, nzn);  
    addscale(1 / th, -(1.0 - th) / th, Ey, Ey_n, Eyth, nxn, nyn, nzn);
    addscale(1 / th, -(1.0 - th) / th, Ez, Ez_n, Ezth, nxn, nyn, nzn);
    
    smoothE_NoComm(Smooth, vct, col);
  }

}


void EMfields3D::initWeightBC_InitialInterpolation(Grid *grid, VirtualTopology3D *vct){
  // this to do the inital field interpolatin, G2R

  NumF_II= 6+ns;
  bool VerboseCheck= false;

  MPI_Status status;

  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  // rank on the local grid
  int rank_local= vct->getCartesian_rank();

  // rank as a child in the parent-child communicator
  int rank_As_Child=-1;
  MPI_Comm CommToParent= vct->getCommToParent(); // if != MPI_COMM_NULL the grid is a child
  if (CommToParent != MPI_COMM_NULL){
    MPI_Comm_rank(CommToParent, &rank_As_Child);
  }

  /* rank as a parent in the parent-child communicator;
     since one grid may have more than one child, I have to calculate every time
     on the right communicatore */
  int rank_As_Parent;
  
  RG_numBCMessages_II= 0;

  /* here, vectors where core 0 of the local child grid assembles the messages
     from all the local grid 
     de-allocated at the end of the function */

  RGBC_struct * RGBC_Info_II_LevelWide;
  int RG_numBCMessages_II_LevelWide=0;

  /* phase 1 */
  // as a child
  if (CommToParent != MPI_COMM_NULL){

    RG_MaxMsgSize_II=0; // values is calulcated in the initWeightBC_Phase1's
    /* instantiate ghost structure */
    RGBC_Info_II = new RGBC_struct[MAX_RG_numBCMessages];

    for (int i=0; i< MAX_RG_numBCMessages; i++){
      RGBC_Info_II[i].ix_first= 0;
    }

    // -1 for ghost
    // 0 for active
    // -10 for II

    // active
    initWeightBC_InitialInterpolation_Phase1(grid, vct);

    // now I have all the info to initialise the BC vectors

    cout << "Before instantiating RG_numBCMessages_II: " << RG_numBCMessages_II << " RG_MaxMsgSize " << RG_MaxMsgSize << endl;
    // rows: [0 - RG_numBCMessages_II]                                                    
    // columns: [0 - RG_MaxMsgSize]                                                          
    Ex_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    Ey_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    Ez_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
         
    Bxn_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    Byn_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    Bzn_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);

    rho0_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    rho1_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    rho2_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    rho3_II= newArr2(double, RG_numBCMessages_II, RG_MaxMsgSize_II);
    
    if (ns > 4){
      cout << "Inside initWeightBC_InitialInterpolation: it is written only for 4 species max, but here i have more... Aborting..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    
    // this used when receiving
    RGMsg_II= new double[RG_MaxMsgSize_II *NumF_II];

    /**** check starts ****/
    // before proceeding, a check with the possibility of aborting if the check is failed
    int PG= vct->getParentGridNum();
    int localRank= vct->getCartesian_rank();
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double MsgLimsXMin, MsgLimsXMax, MsgLimsYMin, MsgLimsYMax, MsgLimsZMin, MsgLimsZMax;

    // check on active
    for (int i=0; i< RG_numBCMessages_II; i++){
      int CG= RGBC_Info_II[i].CG_core;
       
      MsgLimsXMin= RGBC_Info_II[i].CG_x_first;
      MsgLimsXMax= RGBC_Info_II[i].CG_x_first+ dx*(RGBC_Info_II[i].np_x-1);

      MsgLimsYMin= RGBC_Info_II[i].CG_y_first;
      MsgLimsYMax= RGBC_Info_II[i].CG_y_first+ dy*(RGBC_Info_II[i].np_y-1);

      MsgLimsZMin= RGBC_Info_II[i].CG_z_first;
      MsgLimsZMax= RGBC_Info_II[i].CG_z_first+ dz*(RGBC_Info_II[i].np_z-1);

      grid->getParentLimits(vct, CG, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

      if (MsgLimsXMin < xmin or MsgLimsXMax > xmax or MsgLimsYMin < ymin or MsgLimsYMax > ymax or MsgLimsZMin < zmin or MsgLimsZMax > zmax){
	cout <<"G" <<numGrid <<"R" << localRank <<" Msg " << i << " we have a problem in initWeightBC, active BC, aborting ... " << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }


    /**** check ends ****/
  } // end if (numGrid>0)
  /* end phase 1 */

  /* phase 2, RG sends info to CG */
  // children grids
  if (CommToParent != MPI_COMM_NULL ){  

    // phase 2a: children assemble all the info in core 0 in the LOCAL child grid, 0-> XLEN*YLEN*ZLEN-1 
    if (rank_local > 0){
      /* send one message more; the last message has -1 in the RG_core
	 to signal end of 'valid' messages */
      /* -1 as tag for ghost, 0 as tag for active */

      // II
      MPI_Send(RGBC_Info_II, RG_numBCMessages_II+1, MPI_RGBC_struct, 0, TAG_II, vct->getComm());
    }

    if  (rank_local==0){

      RGBC_Info_II_LevelWide = new RGBC_struct[MAX_size_LevelWide];
      initWeightBC_Phase2a(grid, vct, RGBC_Info_II_LevelWide, &RG_numBCMessages_II_LevelWide, RGBC_Info_II, RG_numBCMessages_II, -10);
       
    } // end if (rank_local==0)
   
 
    // phase 2b: core 0 of the child grid assembles & sends messages for all CG cores
    if (rank_local==0){
      // II
      initWeightBC_Phase2b(grid, vct, RGBC_Info_II_LevelWide, RG_numBCMessages_II_LevelWide, -10);
    } // end if (rank_local==0), phase 2b
  } // end if on children grid
  

  MPI_Barrier(MPI_COMM_WORLD);
  if (vct->getSystemWide_rank()==0){
    cout << "After barrier phase 2" << endl;
  }

  // phase 2c, only for the parent
  // each CG core will receive a message per child grid
  if (numChildren > 0 ){

    CG_MaxSizeMsg_II=0;

    //active
    CG_Info_II= newArr2(RGBC_struct, numChildren, MAX_RG_numBCMessages);
    CG_numBCMessages_II= new int[numChildren]; 

    for (int i=0; i<numChildren; i++){
      CG_numBCMessages_II[i]=0;
    }

    initWeightBC_Phase2c(grid, vct, CG_Info_II, CG_numBCMessages_II, -10);

    // same size for all children, used to build BC msg
    CGMsg_II = new double [CG_MaxSizeMsg_II * NumF_II];

  } // end  if (numChildren > 0 ), only for parents
  
  Trim_RGBC_Vectors_II(vct);

  /* these deletes only for core 0 of RGs */  
  if (CommToParent != MPI_COMM_NULL && rank_local==0){
    delete[] RGBC_Info_II_LevelWide;

  }
  
  if (numGrid >0){
    cout << "Grid " <<numGrid << " R " << vct->getCartesian_rank() <<" receives " <<RG_numBCMessages_II << " II msgs" << endl;
    for (int i=0; i< RG_numBCMessages_II; i++){
      cout <<"msg " << i << " of " << RG_numBCMessages_II <<" from " << RGBC_Info_II[i].CG_core << " size: " << RGBC_Info_II[i].np_x*RGBC_Info_II[i].np_y*RGBC_Info_II[i].np_z << " RGBC_Info_II[i].MsgID  " << RGBC_Info_II[i].MsgID << " RGBC_Info_II[i].CG_core " << RGBC_Info_II[i].CG_core << " RGBC_Info_II[i].RG_core " << RGBC_Info_II[i].RG_core <<endl;
    }
  }
  for (int ch=0; ch < numChildren; ch++){
    cout <<"Grid "<< numGrid << " R " << vct->getCartesian_rank()<< " sends " << CG_numBCMessages_II[ch] <<" II msgs" << endl;
    for (int m=0; m< CG_numBCMessages_II[ch]; m++){
      cout <<"msg "<< m <<" of " <<CG_numBCMessages_II[ch] << " to " << CG_Info_II[ch][m].RG_core << " size: " <<CG_Info_II[ch][m].np_x*CG_Info_II[ch][m].np_y*CG_Info_II[ch][m].np_z << " CG_Info_II[ch][m].MsgID " <<CG_Info_II[ch][m].MsgID << " CG_Info_II[ch][m].CG_core " <<CG_Info_II[ch][m].CG_core << " CG_Info_II[ch][m].RG_core " <<CG_Info_II[ch][m].RG_core  << endl;
    }
  }
  
}

void EMfields3D::initWeightBC_InitialInterpolation_Phase1(Grid *grid, VirtualTopology3D *vct){

  // RG cores build the map for proj info

  // send everything
  int i_s=0, i_e= nxn-1;
  int j_s=0, j_e= nyn-1;
  int k_s=0, k_e= nzn-1;
    
  // policy:
  // explore Z dir
  // for on the number of cores found there: explore Y dir
  // for on the number of cores found there: explore X dir
  // finally, commit  NB: all faces should have the same c

  int MS= nxn; if (nyn>MS) MS= nyn; if (nzn>MS) MS= nzn;
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
  grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

  for (int n=0; n<Dir1_Ncores; n++){ // it will find again the core in Dir 1, but it will also explore Dir 2 
    // Y dir / Dir 2
    grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n],  Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); 
    
    for (int m=0; m< Dir2_Ncores; m++){ //it will find again the core in Dir 1, but it will also explore Dir 2 
      // X dir / Dir 3
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e, Dir2_IndexFirstPointperC[m],  Dir1_IndexFirstPointperC[n], Dir3_SPXperC, Dir3_SPYperC, Dir3_SPZperC, Dir3_NPperC, Dir3_rank, &Dir3_Ncores, Dir3_IndexFirstPointperC);

      for (int NN=0; NN< Dir3_Ncores; NN++){
	// using the function written for BCs, it's the same

	//void Assign_RGBC_struct_Values(RGBC_struct *s, int ix_first_tmp, int iy_first_tmp, int iz_first_tmp, int BCside_tmp, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp, int ID);

	Assign_RGBC_struct_Values(RGBC_Info_II + RG_numBCMessages_II, Dir3_IndexFirstPointperC[NN], Dir2_IndexFirstPointperC[m], Dir1_IndexFirstPointperC[n], -1, Dir3_NPperC[NN], Dir2_NPperC[m], Dir1_NPperC[n], Dir3_SPXperC[NN], Dir3_SPYperC[NN], Dir3_SPZperC[NN], Dir3_rank[NN], rank_CTP, RG_numBCMessages_II);
	
	RG_numBCMessages_II++;
	
	int tmp= Dir3_NPperC[NN]*Dir2_NPperC[m]*Dir1_NPperC[n];
	if (tmp > RG_MaxMsgSize_II) RG_MaxMsgSize_II= tmp; 
	
      } // end for (int NN=0; NN< Dir3_Ncores; NN++){
    } // end  for (int m=0; m< Dir2_Ncores; m++){
  } // end for (int n=0; n<Dir1_Ncores; n++){
  // end here, the same as Explore3DAndCommit for particles


  // the msg signaling the end
  RGBC_Info_II[RG_numBCMessages_II].RG_core= -1;
  RGBC_Info_II[RG_numBCMessages_II].CG_core= -1;
  

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

void EMfields3D::sendInitialInterpolation(Grid *grid, VirtualTopology3D *vct){
  // this is only for coarses grids 
  if (numChildren < 1){ return; }
  
  int msg; // this is the fake msg that will be sent to test the communication
  int rank_As_Parent;
  int dest; // destination on the parent-child communicator
  MPI_Comm PC_Comm; // the parent-child communicator
  int tag;

  for (int ch=0; ch< numChildren; ch ++){

    // active
    for (int m=0; m<CG_numBCMessages_II[ch]; m++){
      sendOneBC(vct, grid, CG_Info_II[ch][m], ch, -10);
    }

  }

  MPI_Barrier(vct->getComm());
  /*if (vct->getCartesian_rank()==0){
    cout << "grid " << numGrid << " Finished sending II" << endl;
    }*/
}

void EMfields3D::buildIIMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo, int Size ){
  
  // NB: if i send more stuff, do all the interpolation here
  // i use the three indexes, knowing that one will be zero

  // position of the initial RG point in the CG grid
  // points are displacement over this

  //cout << "Building msg starting from " << RGInfo.CG_x_first <<", " << RGInfo.CG_y_first << ", " <<RGInfo.CG_z_first <<" size "<< RGInfo.np_x << ", "  << RGInfo.np_y <<", " << RGInfo.np_z << endl;

  double x0= RGInfo.CG_x_first;
  double y0= RGInfo.CG_y_first;
  double z0= RGInfo.CG_z_first;

  double inv_dx= 1./dx;
  double inv_dy= 1./dy;
  double inv_dz= 1./dz;

  double xp, yp, zp;
  int count =0;

  double *RHO= new double[ns];

  //cout << "inside building msg RGInfo.np_x " << RGInfo.np_x << " RGInfo.np_y " << RGInfo.np_y << " RGInfo.np_z " << RGInfo.np_z <<endl;

  //cout <<"building msg, starts @ " << x0 << "; " << y0 << "; " << z0 << endl;

  for (int i= 0; i<RGInfo.np_x; i++){
    for (int j=0; j<RGInfo.np_y; j++){
      for (int k=0; k<RGInfo.np_z; k++){
	
	// ok, now each RG point is treated as a particle here
	// and i need the E at the point position

	xp= x0 + i*dx_Ch[ch];
	yp= y0 + j*dy_Ch[ch];
	zp= z0 + k*dz_Ch[ch]; 
	
	//cout << "building msg, point " << xp << "; " <<yp <<"; " << zp << endl;

	// this is copied from the mover
	const double ixd = floor((xp - xStart) * inv_dx);
	const double iyd = floor((yp - yStart) * inv_dy);
	const double izd = floor((zp - zStart) * inv_dz);
	int ix = 2 + int (ixd);
	int iy = 2 + int (iyd);
	int iz = 2 + int (izd);
	if (ix < 1)
	  ix = 1;
	if (iy < 1)
	  iy = 1;
	if (iz < 1)
	  iz = 1;
	if (ix > nxn - 1)
	  ix = nxn - 1;
	if (iy > nyn - 1)
	  iy = nyn - 1;
	if (iz > nzn - 1)
	  iz = nzn - 1;

	double xi  [2];
	double eta [2];
	double zeta[2];

	xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
	eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
	zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
	xi  [1] = grid->getXN(ix,iy,iz) - xp;
	eta [1] = grid->getYN(ix,iy,iz) - yp;
	zeta[1] = grid->getZN(ix,iy,iz) - zp;

	/*cout << "point: " << xp <<", " << yp <<", " << zp << endl;
	cout << "xi[0] " << xi[0] << " xi[1] " << xi[1] << " eta[0] " << eta[0] << " eta[1] " << eta[1] << " zeta[0] " << zeta[0] << " zeta[1] " << zeta[1];
	cout << "invVOL: " << invVOL << endl;*/

	double Exl = 0.0;
	double Eyl = 0.0;
	double Ezl = 0.0;

	double Bxl = 0.0;
	double Byl = 0.0;
	double Bzl = 0.0;

	double rho0l = 0.0;
	double rho1l = 0.0;
	double rho2l = 0.0;
	double rho3l = 0.0;

	const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
	const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
	const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
	const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
	const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
	const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
	const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
	const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
	
	//                                                                                                                                                       
	Exl += weight000 * (Ex[ix][iy][iz]             );
	Exl += weight001 * (Ex[ix][iy][iz - 1]         );
	Exl += weight010 * (Ex[ix][iy - 1][iz]         );
	Exl += weight011 * (Ex[ix][iy - 1][iz - 1]     );
	Exl += weight100 * (Ex[ix - 1][iy][iz]         );
	Exl += weight101 * (Ex[ix - 1][iy][iz - 1]     );
	Exl += weight110 * (Ex[ix - 1][iy - 1][iz]     );
	Exl += weight111 * (Ex[ix - 1][iy - 1][iz - 1] );

	Eyl += weight000 * (Ey[ix][iy][iz]             );
	Eyl += weight001 * (Ey[ix][iy][iz - 1]         );
	Eyl += weight010 * (Ey[ix][iy - 1][iz]         );
	Eyl += weight011 * (Ey[ix][iy - 1][iz - 1]     );
	Eyl += weight100 * (Ey[ix - 1][iy][iz]         );
	Eyl += weight101 * (Ey[ix - 1][iy][iz - 1]     );
	Eyl += weight110 * (Ey[ix - 1][iy - 1][iz]     );
	Eyl += weight111 * (Ey[ix - 1][iy - 1][iz - 1] );

	Ezl += weight000 * (Ez[ix][iy][iz]             );
	Ezl += weight001 * (Ez[ix][iy][iz - 1]         );
	Ezl += weight010 * (Ez[ix][iy - 1][iz]         );
	Ezl += weight011 * (Ez[ix][iy - 1][iz - 1]     );
	Ezl += weight100 * (Ez[ix - 1][iy][iz]         );
	Ezl += weight101 * (Ez[ix - 1][iy][iz - 1]     );
	Ezl += weight110 * (Ez[ix - 1][iy - 1][iz]     );
	Ezl += weight111 * (Ez[ix - 1][iy - 1][iz - 1] );

	// B
	Bxl += weight000 * (Bxn[ix][iy][iz]             );
	Bxl += weight001 * (Bxn[ix][iy][iz - 1]         );
	Bxl += weight010 * (Bxn[ix][iy - 1][iz]         );
	Bxl += weight011 * (Bxn[ix][iy - 1][iz - 1]     );
	Bxl += weight100 * (Bxn[ix - 1][iy][iz]         );
	Bxl += weight101 * (Bxn[ix - 1][iy][iz - 1]     );
	Bxl += weight110 * (Bxn[ix - 1][iy - 1][iz]     );
	Bxl += weight111 * (Bxn[ix - 1][iy - 1][iz - 1] );

	Byl += weight000 * (Byn[ix][iy][iz]             );
	Byl += weight001 * (Byn[ix][iy][iz - 1]         );
	Byl += weight010 * (Byn[ix][iy - 1][iz]         );
	Byl += weight011 * (Byn[ix][iy - 1][iz - 1]     );
	Byl += weight100 * (Byn[ix - 1][iy][iz]         );
	Byl += weight101 * (Byn[ix - 1][iy][iz - 1]     );
	Byl += weight110 * (Byn[ix - 1][iy - 1][iz]     );
	Byl += weight111 * (Byn[ix - 1][iy - 1][iz - 1] );

	Bzl += weight000 * (Bzn[ix][iy][iz]             );
	Bzl += weight001 * (Bzn[ix][iy][iz - 1]         );
	Bzl += weight010 * (Bzn[ix][iy - 1][iz]         );
	Bzl += weight011 * (Bzn[ix][iy - 1][iz - 1]     );
	Bzl += weight100 * (Bzn[ix - 1][iy][iz]         );
	Bzl += weight101 * (Bzn[ix - 1][iy][iz - 1]     );
	Bzl += weight110 * (Bzn[ix - 1][iy - 1][iz]     );
	Bzl += weight111 * (Bzn[ix - 1][iy - 1][iz - 1] );
	// end B

	for (int is=0; is <ns ; is ++){ RHO[is]=0;}

	for (int is=0; is <ns; is++){
	  RHO[is] += weight000 * (rhons[is][ix][iy][iz]             );
	  RHO[is] += weight001 * (rhons[is][ix][iy][iz - 1]         );
	  RHO[is] += weight010 * (rhons[is][ix][iy - 1][iz]         );
	  RHO[is] += weight011 * (rhons[is][ix][iy - 1][iz - 1]     );
	  RHO[is] += weight100 * (rhons[is][ix - 1][iy][iz]         );
	  RHO[is] += weight101 * (rhons[is][ix - 1][iy][iz - 1]     );
	  RHO[is] += weight110 * (rhons[is][ix - 1][iy - 1][iz]     );
	  RHO[is] += weight111 * (rhons[is][ix - 1][iy - 1][iz - 1] );
	}
	
	
	CGMsg_II[0*Size +count]= Exl;
	CGMsg_II[1*Size +count]= Eyl;
	CGMsg_II[2*Size +count]= Ezl;

	CGMsg_II[3*Size +count]= Bxl;
	CGMsg_II[4*Size +count]= Byl;
	CGMsg_II[5*Size +count]= Bzl;

	for (int is=0; is<ns; is ++){
	  int NN= 6+is;
	  CGMsg_II[NN*Size +count]= RHO[is];
	}

	count ++;

      }
    }
  }
  delete[]RHO;

}

void EMfields3D::ApplyInitialInterpolation(VirtualTopology3D *vct, Grid * grid){

  if (vct->getCommToParent() != MPI_COMM_NULL) {

    for (int m=0; m< RG_numBCMessages_II; m++){

      int II= RGBC_Info_II[m].np_x;
      int JJ= RGBC_Info_II[m].np_y;
      int KK= RGBC_Info_II[m].np_z;

      int ix_f= RGBC_Info_II[m].ix_first;
      int iy_f= RGBC_Info_II[m].iy_first;
      int iz_f= RGBC_Info_II[m].iz_first;

      int count=0;
      for (int i= 0; i<II; i++)
	for (int j=0; j<JJ; j++)
	  for (int k=0; k<KK; k++){

	    //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize);                   
	    Ex[ix_f + i][iy_f + j][iz_f + k]= Ex_II[m][count];
	    Ey[ix_f + i][iy_f + j][iz_f + k]= Ey_II[m][count];
	    Ez[ix_f + i][iy_f + j][iz_f + k]= Ez_II[m][count];

	    Bxn[ix_f + i][iy_f + j][iz_f + k]= Bxn_II[m][count];
	    Byn[ix_f + i][iy_f + j][iz_f + k]= Byn_II[m][count];
	    Bzn[ix_f + i][iy_f + j][iz_f + k]= Bzn_II[m][count];

	    if (ns>0)
	      rhons[0][ix_f + i][iy_f + j][iz_f + k]= rho0_II[m][count]	;
	    if (ns>1)
	      rhons[1][ix_f + i][iy_f + j][iz_f + k]= rho1_II[m][count]	;
	    if (ns>2)
	      rhons[2][ix_f + i][iy_f + j][iz_f + k]= rho2_II[m][count]	;
	    if (ns>3)
	      rhons[3][ix_f + i][iy_f + j][iz_f + k]= rho3_II[m][count]	;

	    count ++;
	    cout << "Ex[ix_f + i][iy_f + j][iz_f + k]: " << Ex[ix_f + i][iy_f + j][iz_f + k] <<endl;
	  }
    }
  
    // initialize B on centers
    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
    
    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
    
    RG_numBCMessages_II=0;
    RG_MaxMsgSize_II=0;
    delete[]RGBC_Info_II;
    delete[]RGMsg_II;
 
   }
}

void EMfields3D::DeallocateII(){
  
  if (numChildren <1) return;

  delArr2(CG_Info_II, numChildren);
  delete[]CG_numBCMessages_II;
  delete[]CGMsg_II;
}

/* index are given from left to right */
double EMfields3D::BufferFactor(char Dir, int X, int Y, int Z, int BufLenX, int BufLenY, int BufLenZ){
  // X, Y, Z are the indexes in the three directions
  // BufLen the length of the buffering area
  
  // direction of the BC Buffer - (L)eft / (R)ight - (B)ottom / (T)op - (F)ront / (b)ack  

  int WHICH;
  double ret;
  int BufLenX_I=4; int BufLenY_I=4; int BufLenZ_I=4;

  // prototype left - index 2 has to get value 1
  // distance (positive)= I-2
  // ret= 1- distance/SP
  // ret= 1.- double(I-2.)/double(SP)
  // prototype right - index nn-3 has to get value 1
  // distance (positive)= (NN-3-I)
  // ret= 1- (distance/SP)
  // ret= 1.- double(NN-3.-I)/double(SP)

  // the node with the "1"
  double start=0;

  if (Dir == 'L'){ // left
    ret= 1.-double(X-start)/double(BufLenX_I);
  } else if (Dir == 'B'){ // bottom
    ret= 1.-double(Z-start)/double(BufLenZ_I);
  } else if (Dir == 'F'){ // front
    ret= 1.-double(Y-start)/double(BufLenY_I);
  } else if (Dir == 'R'){ // right
    ret= 1.- double(nxn-1-start-X)/double(BufLenX_I);
  } else if (Dir == 'T'){ // top
    ret= 1.- double(nzn-1-start-Z)/double(BufLenZ_I);
  } else if (Dir == 'b'){ // back
    ret= 1.- double(nyn-1-start-Y)/double(BufLenY_I);
  }
  else{
    cout <<"Dramatic error in BufferFactor, I am sending unrecongnised side " <<Dir <<". Aborting ...";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (ret >1. or ret < 0.){
    cout <<"Dramatic error in BufferFactor, size "<<  Dir <<" i get " <<  ret <<". Aborting ...";
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  return ret;
    
}

void EMfields3D::averageBC_BufferSource(VirtualTopology3D * vct, double *** SX, double *** SY, double *** SZ, double *** vectX, double *** vectY, double *** vectZ, double ** Fx_BC, double ** Fy_BC, double ** Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages){


  double Fac;
  
  for (int m=0; m<RG_numBCMessages; m++){
  
    // one direction has to be 1
    int N0=0;
    if (!(RGBC_Info[m].np_x==1)) N0++;
    if (!(RGBC_Info[m].np_y==1)) N0++;
    if (!(RGBC_Info[m].np_z==1)) N0++;

    //if (N0 !=2) {cout << "WARNING: in setBC_Nodes N0!= 2 (but it may be ok)" << endl; }

    int II= RGBC_Info[m].np_x;
    int JJ= RGBC_Info[m].np_y;
    int KK= RGBC_Info[m].np_z;

    int ix_f= RGBC_Info[m].ix_first;
    int iy_f= RGBC_Info[m].iy_first;
    int iz_f= RGBC_Info[m].iz_first;

    int count=0;
    for (int i= 0; i<II; i++){
      for (int j=0; j<JJ; j++){
	for (int k=0; k<KK; k++){

	  Fac= BufferFactor(DirBuffer[m], ix_f + i, iy_f + j, iz_f + k, BufX, BufY, BufZ);

	  
	  //Ex_Ghost_BC= newArr2(double, RG_numBCMessages_Ghost, RG_MaxMsgSize); 	  
	  SX[ix_f + i][iy_f + j][iz_f + k]= Fac*Fx_BC[m][count] + (1.-Fac)* vectX[ix_f + i][iy_f + j][iz_f + k];
	  SY[ix_f + i][iy_f + j][iz_f + k]= Fac*Fy_BC[m][count] + (1.-Fac)* vectY[ix_f + i][iy_f + j][iz_f + k];
	  SZ[ix_f + i][iy_f + j][iz_f + k]= Fac*Fz_BC[m][count] + (1.-Fac)* vectZ[ix_f + i][iy_f + j][iz_f + k];
	  
	  count ++;

	} // end cycle k
      } // end cycle j
    } // end cycle i
    //cout << "set BC Nodes: count " << count << endl;
  } // end cycle on msg

  return;

}

void EMfields3D::MPI_RGBC_struct_commit(){

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

}

void EMfields3D::copyMoments(Grid * grid, VirtualTopology3D * vct, double ***P_rho, double ***P_Jx, double ***P_Jy, double ***P_Jz, double ***P_pxx, double ***P_pxy, double ***P_pxz, double ***P_pyy, double ***P_pyz, double ***P_pzz, int is){

  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++){
	P_rho[i][j][k]= rhons[is][i][j][k];
	
	P_Jx[i][j][k]=Jxs[is][i][j][k];
	P_Jy[i][j][k]=Jys[is][i][j][k];
	P_Jz[i][j][k]=Jzs[is][i][j][k];

	P_pxx[i][j][k]=pXXsn[is][i][j][k];
	P_pxy[i][j][k]=pXYsn[is][i][j][k];
	P_pxz[i][j][k]=pXZsn[is][i][j][k];
	P_pyy[i][j][k]=pYYsn[is][i][j][k];
	P_pyz[i][j][k]=pYZsn[is][i][j][k];
	P_pzz[i][j][k]=pZZsn[is][i][j][k];
      }
	
  smooth(Smooth, P_rho, 1, grid, vct);                                                                                                      
  smooth(Smooth, P_Jx, 1, grid, vct);  
  smooth(Smooth, P_Jy, 1, grid, vct);  
  smooth(Smooth, P_Jz, 1, grid, vct);                                                                
  smooth(Smooth, P_pxx, 1, grid, vct);
  smooth(Smooth, P_pxy, 1, grid, vct);  
  smooth(Smooth, P_pxz, 1, grid, vct);     
  smooth(Smooth, P_pyy, 1, grid, vct);   
  smooth(Smooth, P_pyz, 1, grid, vct);  
  smooth(Smooth, P_pzz, 1, grid, vct);

}


void EMfields3D::PostInit(){
  for (int is=0; is < ns; is++)
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
        for (int k=0; k < nzc; k++){

	  RHOINIT[is][i][j][k] = rhocs[is][i][j][k];

	}
}

void EMfields3D::smoothFaces(double value, double *** F, Grid * grid, VirtualTopology3D * vct, int type){
  // smooth only the faces, first/ last active node

  if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

  if (numGrid==0) return;
  int valueSM= 0.0;

  value= valueSM;

  int i_s, i_e, j_s, j_e, k_s, k_e;

  int nxn, nyn, nzn; // redifined
  switch (type) {
  case (0):
    nxn = grid->getNXC();
    nyn = grid->getNYC();
    nzn = grid->getNZC();
    
    break;
  case (1):
    nxn = grid->getNXN();
    nyn = grid->getNYN();
    nzn = grid->getNZN();
    break;
  }



  if (vct->getXleft_neighbor_P()== MPI_PROC_NULL) i_s= 2; else i_s=1;
  if (vct->getYleft_neighbor_P()== MPI_PROC_NULL) j_s= 2; else j_s=1;
  if (vct->getZleft_neighbor_P()== MPI_PROC_NULL) k_s= 2; else k_s=1;
  if (vct->getXright_neighbor_P()== MPI_PROC_NULL) i_e= nxn-2; else i_e= nxn-1;
  if (vct->getYright_neighbor_P()== MPI_PROC_NULL) j_e= nyn-2; else j_e= nyn-1;
  if (vct->getZright_neighbor_P()== MPI_PROC_NULL) k_e= nzn-2; else k_e= nzn-1;

  /* implementation 2: 2D smoothing */
  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;

      switch (type) {
      case (0):
	communicateCenterBoxStencil(nxn, nyn, nzn, F, vct);	
	break;
      case (1):
	communicateNodeBoxStencil(nxn, nyn, nzn, F, vct);
	break;
      }
      

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      //alpha = (1.0 - value) / 6;
      alpha = (1.0 - value) / 4;
      // left face
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j - 1][k] + F[i][j + 1][k] + F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

	      } // end left face

      // right face
      if (vct->getXright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[i][j][k] + alpha * ( F[i][j - 1][k] + F[i][j + 1][k] + F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

	      } // end right face

      // front face
      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i - 1][j][k] + F[i + 1][j][k] +  F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } // end front face

      // back face
      if (vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn - 1; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i - 1][j][k] + F[i + 1][j][k] + F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn - 1; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } // end back face

      // bottom face
      if (vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i - 1][j][k] + F[i + 1][j][k] + F[i][j - 1][k] + F[i][j + 1][k]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[i][j][k] = temp[i][j][k];

	      } // end bottom face
      
      // top face
      if (vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn - 1; k++)

	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i - 1][j][k] + F[i + 1][j][k] + F[i][j - 1][k] + F[i][j + 1][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn - 1; k++)
	      F[i][j][k] = temp[i][j][k];

      } // end top face
      
      // gli spigoli
      alpha = (1.0 - value) / 2;
      
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * ( F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = 1; i < 2; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } 

      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j - 1][k] + F[i][j + 1][k] );

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[i][j][k] = temp[i][j][k];

      }

      
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j][k - 1] + F[i][j][k + 1]);


	for (int i = 1; i < 2; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } 

      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * ( F[i][j - 1][k] + F[i][j + 1][k]);

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[i][j][k] = temp[i][j][k];

      }

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j][k - 1] + F[i][j][k + 1]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } 

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j - 1][k] + F[i][j + 1][k]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[i][j][k] = temp[i][j][k];

      }

      
      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j][k - 1] + F[i][j][k + 1]);


	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      F[i][j][k] = temp[i][j][k];

      } 

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i][j - 1][k] + F[i][j + 1][k]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[i][j][k] = temp[i][j][k];

      }

      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i-1][j][k] + F[i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = 1; k < 2; k++)
	      F[i][j][k] = temp[i][j][k];

      }
      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i-1][j][k] + F[i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[i][j][k] = temp[i][j][k];

      }

      if (vct->getYright_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i-1][j][k] + F[i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = 1; k < 2; k++)
	      F[i][j][k] = temp[i][j][k];

      }
      if (vct->getYright_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[i][j][k] + alpha * (F[i-1][j][k] + F[i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[i][j][k] = temp[i][j][k];

      }
      
      // end gli spigoli

      delArr3(temp, nxn, nyn);
    }

  }
  
  /* end implementation 2: 2D smoothing */
}




void EMfields3D::smoothFaces(double value, double **** F, Grid * grid, VirtualTopology3D * vct, int is, int type){
  // smooth only the faces, first/ last active node
  // NB: ghost nodes in moments are NOT fixed, so i cannot go and use indeces 0 and nn-1;
  // that's why 2 --> nn-2

  if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

  return;
  value=0.5;
  if (numGrid==0) return;

  int i_s, i_e, j_s, j_e, k_s, k_e;

  int nxn, nyn, nzn; // redifined
  switch (type) {
  case (0):
    nxn = grid->getNXC();
    nyn = grid->getNYC();
    nzn = grid->getNZC();
    
    break;
  case (1):
    nxn = grid->getNXN();
    nyn = grid->getNYN();
    nzn = grid->getNZN();
    break;
  }


  if (vct->getXleft_neighbor_P()== MPI_PROC_NULL) i_s= 2; else i_s=1;
  if (vct->getYleft_neighbor_P()== MPI_PROC_NULL) j_s= 2; else j_s=1;
  if (vct->getZleft_neighbor_P()== MPI_PROC_NULL) k_s= 2; else k_s=1;
  if (vct->getXright_neighbor_P()== MPI_PROC_NULL) i_e= nxn-2; else i_e= nxn-1;
  if (vct->getYright_neighbor_P()== MPI_PROC_NULL) j_e= nyn-2; else j_e= nyn-1;
  if (vct->getZright_neighbor_P()== MPI_PROC_NULL) k_e= nzn-2; else k_e= nzn-1;

  /* implementation 2: 2D smoothing */
  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;
   

      switch (type) {
      case (0):
	communicateCenterBoxStencil(nxn, nyn, nzn, F[is], vct);
	
	break;
      case (1):
	communicateNodeBoxStencil(nxn, nyn, nzn, F[is], vct);
	break;
      }
      

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      //alpha = (1.0 - value) / 6;
      alpha = (1.0 - value) / 4;
      // left face
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j - 1][k] + F[is][i][j + 1][k] + F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

	      } // end left face

      // right face
      if (vct->getXright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[is][i][j][k] + alpha * ( F[is][i][j - 1][k] + F[is][i][j + 1][k] + F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

	      } // end right face

      // front face
      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i - 1][j][k] + F[is][i + 1][j][k] +  F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } // end front face

      // back face
      if (vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn - 1; j++)
	    for (int k = k_s; k < k_e; k++)

	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i - 1][j][k] + F[is][i + 1][j][k] + F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn - 1; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } // end back face

      // bottom face
      if (vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i - 1][j][k] + F[is][i + 1][j][k] + F[is][i][j - 1][k] + F[is][i][j + 1][k]);

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[is][i][j][k] = temp[i][j][k];

	      } // end bottom face
      
      // top face
      if (vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn - 1; k++)

	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i - 1][j][k] + F[is][i + 1][j][k] + F[is][i][j - 1][k] + F[is][i][j + 1][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn - 1; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } // end top face
      
      // gli spigoli
      alpha = (1.0 - value) / 2;
      
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * ( F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = 1; i < 2; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } 

      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j - 1][k] + F[is][i][j + 1][k] );

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }

      
      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j][k - 1] + F[is][i][j][k + 1]);


	for (int i = 1; i < 2; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } 

      if (vct->getXleft_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * ( F[is][i][j - 1][k] + F[is][i][j + 1][k]);

	for (int i = 1; i < 2; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getYleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j][k - 1] + F[is][i][j][k + 1]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } 

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j - 1][k] + F[is][i][j + 1][k]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = 1; k < 2; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }

      
      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getYright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j][k - 1] + F[is][i][j][k + 1]);


	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = nyn-2; j < nyn-1; j++) 
	    for (int k = k_s; k < k_e; k++)
	      F[is][i][j][k] = temp[i][j][k];

      } 

      if (vct->getXright_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i][j - 1][k] + F[is][i][j + 1][k]);

	for (int i = nxn-2; i < nxn-1; i++)
	  for (int j = j_s; j < j_e; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }

      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i-1][j][k] + F[is][i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = 1; k < 2; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }
      if (vct->getYleft_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i-1][j][k] + F[is][i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = 1; j < 2; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }

      if (vct->getYright_neighbor_P() == MPI_PROC_NULL and  vct->getZleft_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = 1; k < 2; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i-1][j][k] + F[is][i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = 1; k < 2; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }
      if (vct->getYright_neighbor_P() == MPI_PROC_NULL and  vct->getZright_neighbor_P() == MPI_PROC_NULL){

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i-1][j][k] + F[is][i+1][j][k] );

	for (int i = i_s; i < i_e; i++)
	  for (int j = nyn-2; j < nyn-1; j++)
	    for (int k = nzn-2; k < nzn-1; k++)
	      F[is][i][j][k] = temp[i][j][k];

      }
      
      // end gli spigoli

      delArr3(temp, nxn, nyn);
    }

  }
  
  /* end implementation 2: 2D smoothing */
}

void EMfields3D::smoothAll(double value, double **** F, Grid * grid, VirtualTopology3D * vct, int is, int type){

  if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

  smoothFaces(value,  F, grid, vct, is, type);

  int i_s, i_e, j_s, j_e, k_s, k_e;
  int nxn, nyn, nzn; // redifined
  switch (type) {
  case (0):
    nxn = grid->getNXC();
    nyn = grid->getNYC();
    nzn = grid->getNZC();
    
    break;
  case (1):
    nxn = grid->getNXN();
    nyn = grid->getNYN();
    nzn = grid->getNZN();
    break;
  }


  // the faces have already been dealt with
  if (vct->getXleft_neighbor_P()== MPI_PROC_NULL) i_s= 2; else i_s=1;
  if (vct->getYleft_neighbor_P()== MPI_PROC_NULL) j_s= 2; else j_s=1;
  if (vct->getZleft_neighbor_P()== MPI_PROC_NULL) k_s= 2; else k_s=1;
  if (vct->getXright_neighbor_P()== MPI_PROC_NULL) i_e= nxn-2; else i_e= nxn-1;
  if (vct->getYright_neighbor_P()== MPI_PROC_NULL) j_e= nyn-2; else j_e= nyn-1;
  if (vct->getZright_neighbor_P()== MPI_PROC_NULL) k_e= nzn-2; else k_e= nzn-1;


  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;
      
      switch (type) {
      case (0):
	communicateCenterBoxStencil(nxn, nyn, nzn, F[is], vct);	
	break;
      case (1):
	communicateNodeBoxStencil(nxn, nyn, nzn, F[is], vct);
	break;
      }
      

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Exth
      for (int i = i_s; i < i_e; i++)
        for (int j = j_s; j < j_e; j++)
          for (int k = k_s; k < k_e; k++)
            temp[i][j][k] = value * F[is][i][j][k] + alpha * (F[is][i - 1][j][k] + F[is][i + 1][j][k] + F[is][i][j - 1][k] + F[is][i][j + 1][k] + F[is][i][j][k - 1] + F[is][i][j][k + 1]);
      for (int i = i_s; i < i_e; i++)
        for (int j = j_s; j < j_e; j++)
          for (int k = k_s; k < k_e; k++)
            F[is][i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
  
  

}

void EMfields3D::smoothAll(double value, double *** F, Grid * grid, VirtualTopology3D * vct, int type){

  if (!SmoothGrid) {
     if (vct->getCartesian_rank() == 0){
       cout << "I am not smoothing Grid " << numGrid << endl;
     }
     return;
   }

  smoothFaces(value,  F, grid, vct, type);

  int i_s, i_e, j_s, j_e, k_s, k_e;

  int nxn, nyn, nzn; // redifined
  switch (type) {
  case (0):
    nxn = grid->getNXC();
    nyn = grid->getNYC();
    nzn = grid->getNZC();
    
    break;
  case (1):
    nxn = grid->getNXN();
    nyn = grid->getNYN();
    nzn = grid->getNZN();
    break;
  }

  // the faces have already been dealt with
  if (vct->getXleft_neighbor_P()== MPI_PROC_NULL) i_s= 2; else i_s=1;
  if (vct->getYleft_neighbor_P()== MPI_PROC_NULL) j_s= 2; else j_s=1;
  if (vct->getZleft_neighbor_P()== MPI_PROC_NULL) k_s= 2; else k_s=1;
  if (vct->getXright_neighbor_P()== MPI_PROC_NULL) i_e= nxn-2; else i_e= nxn-1;
  if (vct->getYright_neighbor_P()== MPI_PROC_NULL) j_e= nyn-2; else j_e= nyn-1;
  if (vct->getZright_neighbor_P()== MPI_PROC_NULL) k_e= nzn-2; else k_e= nzn-1;

  
  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;
      
      switch (type) {
      case (0):
	communicateCenterBoxStencil(nxn, nyn, nzn, F, vct);	
	break;
      case (1):
	communicateNodeBoxStencil(nxn, nyn, nzn, F, vct);
	break;
      }
      

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Exth
      for (int i = i_s; i < i_e; i++)
        for (int j = j_s; j < j_e; j++)
          for (int k = k_s; k < k_e; k++)
            temp[i][j][k] = value * F[i][j][k] + alpha * (F[i - 1][j][k] + F[i + 1][j][k] + F[i][j - 1][k] + F[i][j + 1][k] + F[i][j][k - 1] + F[i][j][k + 1]);
      for (int i = i_s; i < i_e; i++)
        for (int j = j_s; j < j_e; j++)
          for (int k = k_s; k < k_e; k++)
            F[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
  
  

}

void EMfields3D::TESTGhost (double **** vec) {

  for (int is=0; is<ns; is++){

    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++)
	vec[is+ns][1][j][k]= vec[is][0][j][k];

    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++)
	vec[is+ns][nxn-2][j][k]= vec[is][nxn-1][j][k];

    for (int i=0; i<nxn; i++)
      for (int k=0; k<nzn; k++)
	vec[is+ns][i][1][k]= vec[is][i][0][k];

    for (int i=0; i<nxn; i++)
      for (int k=0; k<nzn; k++)
	vec[is+ns][i][nyn-2][k]= vec[is][i][nyn-1][k];

    for (int i=0; i<nxn; i++)
      for (int j=0; j<nyn; j++)
	vec[is+ns][i][j][1]= vec[is][i][j][0];

    for (int i=0; i<nxn; i++)
      for (int j=0; j<nyn; j++)
	vec[is+ns][i][j][nzn-2]= vec[is][i][j][nzn-1];
  }

}

void EMfields3D::TESTGhost (double **** dest, double **** source) {

  for (int is=0; is<ns; is++){

    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++)
	dest[is+ns][1][j][k]= source[is][0][j][k];

    for (int j=0; j<nyn; j++)
      for (int k=0; k<nzn; k++)
	dest[is+ns][nxn-2][j][k]= source[is][nxn-1][j][k];

    for (int i=0; i<nxn; i++)
      for (int k=0; k<nzn; k++)
	dest[is+ns][i][1][k]= source[is][i][0][k];

    for (int i=0; i<nxn; i++)
      for (int k=0; k<nzn; k++)
	dest[is+ns][i][nyn-2][k]= source[is][i][nyn-1][k];

    for (int i=0; i<nxn; i++)
      for (int j=0; j<nyn; j++)
	dest[is+ns][i][j][1]= source[is][i][j][0];

    for (int i=0; i<nxn; i++)
      for (int j=0; j<nyn; j++)
	dest[is+ns][i][j][nzn-2]= source[is][i][j][nzn-1];
  }

}

void EMfields3D::correctDivB(Grid * grid, VirtualTopology3D * vct){

  // correctio only for higher order grids
  if (numGrid==0) return;

  /* the idea is to do this in the RG, after passing ghost cell BC,
     so that I do not need to correct the first three cells */
  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid <<": DIVERGENCE B CLEANING ***" << endl;

  // the nodes on which to work (***included***)
  // if internal core
  iStart_BP=1; iEnd_BP= nxn-2; lenX_BP=nxn-2;
  jStart_BP=1; jEnd_BP= nyn-2; lenY_BP=nyn-2;
  kStart_BP=1; kEnd_BP= nzn-2; lenZ_BP=nzn-2;
  // if bordering MPI_PROC_NULL, consider I want to start from node=3
  // and end at node N-4 (**included**)
  if (vct->getXleft_neighbor()== MPI_PROC_NULL)
    iStart_BP= 3; lenX_BP= lenX_BP-2;
  if (vct->getXright_neighbor()== MPI_PROC_NULL)
    iEnd_BP= nxn-4; lenX_BP= lenX_BP-2;
  if (vct->getYleft_neighbor()== MPI_PROC_NULL)
    jStart_BP= 3; lenY_BP= lenY_BP-2;
  if (vct->getYright_neighbor()== MPI_PROC_NULL)
    jEnd_BP= nyn-4; lenY_BP= lenY_BP-2;
  if (vct->getZleft_neighbor()== MPI_PROC_NULL)
    kStart_BP= 3; lenZ_BP= lenZ_BP-2;
  if (vct->getZright_neighbor()== MPI_PROC_NULL)
    kEnd_BP= nzn-4; lenZ_BP= lenZ_BP-2;
  if (lenX_BP <1 or lenY_BP <1 or lenZ_BP <1){
    cout << "Fatal error: disable divB correction or increase the number of cells per core" <<endl << "ABORTING NOW..." << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid <<": lenX= "<< lenX_BP << ": lenY= "<< lenY_BP  <<": lenZ= "<< lenZ_BP<<" ***" << endl;

  double ***divB = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIX = newArr3(double, nxc, nyc, nzc);
  double ***gradPHIY = newArr3(double, nxc, nyc, nzc);
  double ***gradPHIZ = newArr3(double, nxc, nyc, nzc);
  double *xkrylovPoisson = new double[lenX_BP * lenY_BP * lenZ_BP];
  double *bkrylovPoisson = new double[lenX_BP * lenY_BP * lenZ_BP];
  double ***PHIB= newArr3(double, nxn, nyn, nzn);
  eqValue(0.0, divB, nxn, nyn, nzn);
  eqValue(0.0, gradPHIX, nxc, nyc, nzc);
  eqValue(0.0, gradPHIY, nxc, nyc, nzc);
  eqValue(0.0, gradPHIZ, nxc, nyc, nzc);
  eqValue(0.0, PHIB, nxn, nyn, nzn);
  eqValue(0.0, bkrylovPoisson, lenX_BP * lenY_BP * lenZ_BP);

  grid->divC2N(divB, Bxc, Byc, Bzc);
  phys2solver(bkrylovPoisson, divB, iStart_BP, iEnd_BP, jStart_BP, jEnd_BP, kStart_BP, kEnd_BP);
  eqValue(0.0, xkrylovPoisson, lenX_BP * lenY_BP * lenZ_BP);
  GMRES(&Field::BPoissonImage, xkrylovPoisson, lenX_BP * lenY_BP * lenZ_BP, bkrylovPoisson, 20, 200, CGtol, grid, vct, this);

  /* here, if I am bordering MPI_PROC_NULL, PHIB is OK starting from none n=3 */
  solver2phys(PHIB, xkrylovPoisson, iStart_BP, iEnd_BP, jStart_BP, jEnd_BP, kStart_BP, kEnd_BP);
  /* I do not need to pass BC */
  communicateNode(nxn, nyn, nzn, PHIB, vct);

  /* calculate the gradient */
  grid->gradN2C(gradPHIX, gradPHIY, gradPHIZ, PHIB) ;

  /* add correction, only at the required cells -- included-- 
     -1 in the end because the End's are calculated for nodes, end here we have ells */
  sub(Bxc, gradPHIX, iStart_BP, iEnd_BP-1, jStart_BP, jEnd_BP-1, kStart_BP, kEnd_BP-1);
  sub(Byc, gradPHIY, iStart_BP, iEnd_BP-1, jStart_BP, jEnd_BP-1, kStart_BP, kEnd_BP-1);
  sub(Bzc, gradPHIZ, iStart_BP, iEnd_BP-1, jStart_BP, jEnd_BP-1, kStart_BP, kEnd_BP-1);

  double sumBxc=0.0;
  double sumByc=0.0;
  double sumBzc=0.0;
  for (int i=1; i<nxc-1; i++)
    for (int j=1; j< nyc-1; j++)
      for (int k=1; k< nzc-1; k++){
	sumBxc+=fabs(Bxc[i][j][k]);
	sumByc+=fabs(Byc[i][j][k]);
	sumBzc+=fabs(Bzc[i][j][k]);
      }
  cout << "sumBxc: " << sumBxc << " sumByc: " <<sumByc << " sumBzc: " <<sumBzc << endl;
  
  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid <<": DIVERGENCE B CLEANING ENDED ***" << endl;

  delArr3(divB, nxn, nyn);
  delArr3(gradPHIX, nxc, nyc); 
  delArr3(gradPHIY, nxc, nyc);
  delArr3(gradPHIZ, nxc, nyc);
  delArr3(PHIB, nxn, nyn);
  delete[]xkrylovPoisson ;
  delete[]bkrylovPoisson ;

}

void EMfields3D::BPoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct){
  /* if RG, I enter here with active nodes (not ghost nodes) OK */
  double ***temp = newArr3(double, nxn, nyn, nzn);
  /* NB: now the image is on nodes */
  double ***im = newArr3(double, nxn, nyn, nzn);
  eqValue(0.0, temp, nxn, nyn, nzn);
  eqValue(0.0, im, nxn, nyn, nzn);
  eqValue(0.0, image, lenX_BP * lenY_BP * lenZ_BP);
  /* move form krylov to physical space */
  solver2phys(temp, vector, iStart_BP, iEnd_BP, jStart_BP, jEnd_BP, kStart_BP, kEnd_BP);
  /* calculate the laplacian */
  // the first active may be NOK
  grid->lapN2N(im, temp, vct);
  // from physical to krylov space
  phys2solver(image, im, iStart_BP, iEnd_BP, jStart_BP, jEnd_BP, kStart_BP, kEnd_BP);

  // deallocate 
  delArr3(temp, nxn, nyn);
  delArr3(im, nxn, nyn);
}

/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::calculateEB_ECSIM(Grid * grid, VirtualTopology3D * vct, Collective *col, int cycle) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** G" << numGrid  << ": E CALCULATION ECSIM ***" << endl;

  int dim = 3*(nxn - 2) * (nyn - 2) * (nzn - 2) + 3*(nxc - 2) * (nyc - 2) * (nzc - 2);

  double *xkrylov = new double[dim];
  double *bkrylov = new double[dim];

  // set to zero all the stuff       
  eqValue(0.0, xkrylov, dim);
  eqValue(0.0, bkrylov, dim);

  if (vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER ECSIM ***" << endl;
  // prepare the source                  

  MaxwellSource_ECSIM(bkrylov, grid, vct, col);
  phys2solver(xkrylov, Ex, Ey, Ez, Bxc, Byc, Bzc, nxn, nyn, nzn);
  // solver                                                       

  GMRES_ECSIM(&Field::MaxwellImage_ECSIM, xkrylov, dim, bkrylov, 20, 200, GMREStol, grid, vct, col, this);

  endEcalc_ECSIM(xkrylov, vct, col, grid);
  // deallocate temporary arrays    

  double maxintE=0.0; double maxintB=0.0;
   double intE, intB;
   for (int i=1; i< nxn-1; i++)
     for (int j=1; j< nyn-1; j++)
       for (int k=1; k< nzn-1; k++){
	 intE= Ex[i][j][k]*Ex[i][j][k]+ Ey[i][j][k]*Ey[i][j][k]+ Ez[i][j][k]*Ez[i][j][k];

	 if (intE> maxintE) maxintE=intE;
       }

   for (int i=1; i< nxc-1; i++)
     for (int j=1; j< nyc-1; j++)
       for (int k=1; k< nzc-1; k++){
	 intB= Bxc[i][j][k]*Bxc[i][j][k]+ Byc[i][j][k]*Byc[i][j][k]+ Bzc[i][j][k]*Bzc[i][j][k];

	 if (intB> maxintB) maxintB=intB;
       }

   cout << "*** G" << numGrid  << ": ECSIM algorithm, Cycle " << cycle << " maxintE: " << maxintE << ", maxintB: " << maxintB <<endl;
            
  delete[]xkrylov;
  delete[]bkrylov;
}

// source ECSIM

// image ECSIM
void EMfields3D::MaxwellSource_ECSIM(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col) {

  // brutal fix in case of no particle
  cout << "*** G" << numGrid << ": ECSIM SOURCE: brutal fix in case of no particle "<< endl;
  // first Source for Curl B eq on nodes-2
  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++)
      for (int k=0; k< nzn; k++){
	tempX[i][j][k]= Ex[i][j][k];
	tempY[i][j][k]= Ey[i][j][k];
	tempZ[i][j][k]= Ez[i][j][k];
      }
	
  // the Source for Curl E eq on centers -2
  for (int i=0; i< nxc; i++)
    for (int j=0; j< nyc; j++)
      for (int k=0; k< nzc; k++){
	tempXC[i][j][k]= Bxc[i][j][k];
	tempYC[i][j][k]= Byc[i][j][k];
	tempZC[i][j][k]= Bzc[i][j][k];
      }

  /* MLMD BC */	
  if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){
    // void EMfields3D::BoundaryConditionsECSIMSource(double *** vectNX, double ***vectNY, double *** vectNZ, double *** vectCX, double *** vectCY, double ***vectCZ, VirtualTopology3D *vct, Grid *grid)
    BoundaryConditionsECSIMSource(tempX, tempY, tempZ, tempXC, tempYC, tempZC, vct, grid);
  }

  phys2solver(bkrylov, tempX, tempY, tempZ, tempXC, tempYC, tempZC, nxn, nyn, nzn);
  return;

  /* NN #ifdef __PROFILING__
  clocks->start(11);
  #endif

  #ifdef ELECTROSTATIC
  eqValue(0, Bxc, nxc, nyc, nzc);
  eqValue(0, Byc, nxc, nyc, nzc);
  eqValue(0, Bzc, nxc, nyc, nzc);
  #else */
  scale(tempXC, Bxc, 1, nxc, nyc, nzc);
  scale(tempYC, Byc, 1, nxc, nyc, nzc);
  scale(tempZC, Bzc, 1, nxc, nyc, nzc);
  // NN #endif
  
  /*smoothJh(Smooth, vct, col);
  
  for (int i=0; i<nxn-0; i++)
    for (int j=0; j<nyn-0; j++)
      for (int k=0; k<nzn-0; k++) {
	double Jx_tot = Jxh[i][j][k] + zeroCurrent*Jx_ext[i][j][k];
	double Jy_tot = Jyh[i][j][k] + zeroCurrent*Jy_ext[i][j][k];
	double Jz_tot = Jzh[i][j][k] + zeroCurrent*Jz_ext[i][j][k];
	tempX[i][j][k] = -4*M_PI*th*dt*Jx_tot*grid->getInvVOLn(i,j,k);
	tempY[i][j][k] = -4*M_PI*th*dt*Jy_tot*grid->getInvVOLn(i,j,k);
	tempZ[i][j][k] = -4*M_PI*th*dt*Jz_tot*grid->getInvVOLn(i,j,k);
      }
  */
  addscale(1, tempX, Ex, nxn, nyn, nzn);
  addscale(1, tempY, Ey, nxn, nyn, nzn);
  addscale(1, tempZ, Ez, nxn, nyn, nzn);
  
  /*NN
  // Poisson correction                                             
  if (col->getPoissonCorrection() == "yes") {
    // calculate the gradient                                 

    if (Cylindrical) grid->gradC2N_cyl(temp2X, temp2Y, temp2Z, Phic);
    else             grid->gradC2N(temp2X, temp2Y, temp2Z, Phic);
    BC_E_Poisson(vct,  temp2X, temp2Y, temp2Z);
    addscale(-1, tempX, temp2X, nxn, nyn, nzn);
    addscale(-1, tempY, temp2Y, nxn, nyn, nzn);
    addscale(-1, tempZ, temp2Z, nxn, nyn, nzn);
  }
  // Langdon correction (simpler alternative to divergence cleanning). 
               
  // A.B. Lagndon. CPC 70 447-450 (1992)                     
  if (col->getLangdonCorrection() != 0) {
    double d = col->getLangdonCorrection();
    if (vct->getCartesian_rank() == 0) printf("--> Langdon correction\n");
    if (Cylindrical) grid->divN2C_cyl(tempC, Ex, Ey, Ez);
    else             grid->divN2C(tempC, Ex, Ey, Ez);

    grid->interpN2C(rhoc, rhon);
    communicateCenterBC(nxc, nyc, nzc, rhoc, 2, 2, 2, 2, 2, 2, vct);

    addscale(-4*M_PI, tempC, rhoc, nxc, nyc, nzc);
    scale(tempC, d, nxc, nyc, nzn);

    if (Cylindrical) grid->gradC2N_cyl(temp2X, temp2Y, temp2Z, tempC);
    else grid->gradC2N(temp2X, temp2Y, temp2Z, tempC);
    addscale(dt, tempX, temp2X, nxn, nyn, nzn);
    addscale(dt, tempY, temp2Y, nxn, nyn, nzn);
    addscale(dt, tempZ, temp2Z, nxn, nyn, nzn);
  }

  if (Cylindrical == 1)  fixBC_Cyl_Source(tempX, tempY, tempZ);
  else { 
  fixBC_Source (vct, col, tempX, tempY, tempZ);
   NNif (col->getCase() == "Dipole"){
    fixBC_SourceC(vct, col, tempXC, tempYC, tempZC);
    fixBC_PlanetSource (vct, col, grid, tempX , tempY , tempZ );
    fixBC_PlanetSourceC(vct, col, grid, tempXC, tempYC, tempZC);
    }
    NN}*/

  // physical space -> Krylov space               
  phys2solver(bkrylov, tempX, tempY, tempZ, tempXC, tempYC, tempZC, nxn, nyn, nzn);
  /*#ifdef __PROFILING__
  clocks->stop(11);
  #endif*/
}

/*! Mapping of Maxwell image to give to solver */
void EMfields3D::MaxwellImage_ECSIM(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  /* NN#ifdef __PROFILING__
  clocks->start(12);
  #endif*/
  int dim = 3*(nxn - 2) * (nyn - 2) * (nzn - 2) + 3*(nxc - 2) * (nyc - 2) * (nzc - 2);
  eqValue(0.0, im, dim);
  //eqValue(0.0, im, 6 * (nxn - 2) * (nyn - 2) * (nzn - 2));             
  // move from krylov space to physical space: E -> vect   B -> temp       
  solver2phys(tempX, tempY, tempZ, tempXC, tempYC, tempZC, vector, nxn, nyn, nzn);

  /*NN #ifdef __PROFILING__
  clocks->start(13);
  #endif*/
  communicateCenterBC(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempX, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempY, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempZ, 2, 2, 2, 2, 2, 2, vct);
  /* NN#ifdef __PROFILING__
  clocks->stop(13);
  #endif*/

  /* NN#ifdef ELECTROSTATIC
  for (int i=1; i<nxn-1; i++)
    for (int j=1; j<nyn-1; j++)
      for (int k=1; k<nzn-1; k++) {
	imageEX[i][j][k] = tempXC[i][j][k];
	imageEY[i][j][k] = tempYC[i][j][k];
	imageEZ[i][j][k] = tempZC[i][j][k];
      }
      #else*/
  // c*th*dt*curl(Eth) + Bth             
  double cthdt = c*th*dt;

  /* NNif (Cylindrical == 1) {
    grid->curlN2C_cyl(imageEX, imageEY, imageEZ, tempX, tempY, tempZ);
    } else {*/
  grid->curlN2C(imageEX, imageEY, imageEZ, tempX, tempY, tempZ);
  //NN }

  scale(imageEX, cthdt, nxc, nyc, nzc);
  scale(imageEY, cthdt, nxc, nyc, nzc);
  scale(imageEZ, cthdt, nxc, nyc, nzc);
  addscale(1, imageEX, tempXC, nxc, nyc, nzc);
  addscale(1, imageEY, tempYC, nxc, nyc, nzc);
  addscale(1, imageEZ, tempZC, nxc, nyc, nzc);
  // NN#endif

  // c*dt*th*curl(Bth) - Eth - dt*th*4*pi*Sum(M Eth)                                         
  /* NN#ifdef ELECTROSTATIC
  eqValue(0, imageBX, nxn, nyn, nzn);
  eqValue(0, imageBY, nxn, nyn, nzn);
  eqValue(0, imageBZ, nxn, nyn, nzn);
  #else*/
  /* NN if (Cylindrical == 1) {
    grid->curlC2N_cyl(imageBX, imageBY, imageBZ, tempXC, tempYC, tempZC);
    } else { */
    grid->curlC2N(imageBX, imageBY, imageBZ, tempXC, tempYC, tempZC);
    // NN}
  scale(imageBX, -cthdt, nxn, nyn, nzn);
  scale(imageBY, -cthdt, nxn, nyn, nzn);
  scale(imageBZ, -cthdt, nxn, nyn, nzn);
  addscale(1, imageBX, tempX, nxn, nyn, nzn);
  addscale(1, imageBY, tempY, nxn, nyn, nzn);
  addscale(1, imageBZ, tempZ, nxn, nyn, nzn);

  //cout << "*** G" << numGrid <<": ECSIM: Butchering stuff in the image" << endl;

  if (vct->getCommToParent() != MPI_COMM_NULL and MLMD_BC){
    /* these are the MLMD BC */
    BoundaryConditionsECSIMImage(imageBX, imageBY, imageBZ, imageEX, imageEY, imageEZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC, nxn, nyn, nzn, vct, grid);
  }

  phys2solver(im, imageBX, imageBY, imageBZ, imageEX, imageEY, imageEZ, nxn, nyn, nzn);
  return;
  
  // this not needed if not particles
  /*
  double fac = dt*th*FourPI;
  for (int i=1; i<nxn-1; i++)
    for (int j=1; j<nyn-1; j++)
      for (int k=1; k<nzn-1; k++) {

	double MEx, MEy, MEz;
        double invV = grid->getInvVOLn(i, j, k);
        productMassE(&MEx, &MEy, &MEz, tempX, tempY, tempZ, i, j, k);
        //imageBX[i][j][k] += invV*fac*MEx;                        
        //imageBY[i][j][k] += invV*fac*MEy;           
        //imageBZ[i][j][k] += invV*fac*MEz;                              
        temp2X[i][j][k] = invV*fac*MEx;
        temp2Y[i][j][k] = invV*fac*MEy;
        temp2Z[i][j][k] = invV*fac*MEz;
      }

  //smooth(Smooth, temp2X, temp2Y, temp2Z, vct, col);                       

  for (int i=1; i<nxn-1; i++)
    for (int j=1; j<nyn-1; j++)
      for (int k=1; k<nzn-1; k++) {
        imageBX[i][j][k] += temp2X[i][j][k];
        imageBY[i][j][k] += temp2Y[i][j][k];
        imageBZ[i][j][k] += temp2Z[i][j][k];
      }
  #endif

  // Lambda damping: if no PoissonCorrection                              
  if (col->getCase()=="MHDUCLA" || col->getCase() == "GEM"){
    // do nothing                                                          
  } else if (col->getPoissonCorrection() != "yes" and damping == 1) {
    sumscalprod(imageEX, 1.0, tempXC, Lambda, nxc, nyc, nzc);
    sumscalprod(imageEY, 1.0, tempYC, Lambda, nxc, nyc, nzc);
    sumscalprod(imageEZ, 1.0, tempZC, Lambda, nxc, nyc, nzc);
    sumscalprod(imageBX, 1.0, tempX, Lambda, nxn, nyn, nzn);
    sumscalprod(imageBY, 1.0, tempY, Lambda, nxn, nyn, nzn);
    sumscalprod(imageBZ, 1.0, tempZ, Lambda, nxn, nyn, nzn);
  }

  if (Cylindrical == 1)  fixBC_Cyl_Image(grid, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
    else {
    fixBC_Image (vct, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
    if (col->getCase() == "Dipole"){
      fixBC_ImageC(vct, imageEX, imageEY, imageEZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
      fixBC_PlanetImage (vct, col, grid, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
      fixBC_PlanetImageC(vct, col, grid, imageEX, imageEY, imageEZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
    }
    }

  */
  // move from physical space to krylov space                                                
  phys2solver(im, imageBX, imageBY, imageBZ, imageEX, imageEY, imageEZ, nxn, nyn, nzn);

  /* NN#ifdef __PROFILING__
  clocks->stop(12);
  #endif */

}


void EMfields3D::endEcalc_ECSIM(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid)
{
  /* NN#ifdef __PROFILING__
  clocks->start(14);
  #endif*/

  //eq(Ex_prev, Ex, nxn, nyn, nzn);          
  //eq(Ey_prev, Ey, nxn, nyn, nzn);                         
  //eq(Ez_prev, Ez, nxn, nyn, nzn);                

  // move from krylov space to physical space                 
  solver2phys(Exth, Eyth, Ezth, tempXC, tempYC, tempZC, xkrylov, nxn, nyn, nzn);

  // Lambda damping 
  /* NN if (col->getCase()=="MHDUCLA" || col->getCase()=="GEM"){
    //do nothing                                                            
  }else{
    if (col->getPoissonCorrection() == "yes" and damping == 1) {
      for (int i=0; i<nxn; i++)
	for (int j=0; j<nyn; j++)
          for (int k=0; k<nzn; k++) {
            Exth[i][j][k] +=  (Ex_poisson[i][j][k] -  Exth[i][j][k])*Lambda[i][j][k];
            Eyth[i][j][k] +=  (Ey_poisson[i][j][k] -  Eyth[i][j][k])*Lambda[i][j][k];
            Ezth[i][j][k] +=  (Ez_poisson[i][j][k] -  Ezth[i][j][k])*Lambda[i][j][k];
          }
      sumscalprod(tempXC, -1, Lambda, tempXC, nxc, nyc, nzc);
      sumscalprod(tempYC, -1, Lambda, tempYC, nxc, nyc, nzc);
      sumscalprod(tempZC, -1, Lambda, tempZC, nxc, nyc, nzc);
    }
    }*/

  addscale(1 / th, -(1.0 - th) / th, Bxc, tempXC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Byc, tempYC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Bzc, tempZC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);
  smoothE(Smooth, vct, col);

  // communicate so the interpolation can have good values                       
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct);

  /* NN#ifdef __PROFILING__
  clocks->stop(14);
  #endif */
}

void EMfields3D::BoundaryConditionsECSIMImage(double*** imageNX, double*** imageNY, double*** imageNZ, double*** imageCX, double*** imageCY, double*** imageCZ, double*** vectNX, double*** vectNY, double*** vectNZ, double***vectCX, double*** vectCY, double***vectCZ, int nxn, int nyn, int nzn, VirtualTopology3D * vct, Grid * grid){
 


  /* imageN is the image on the nodes
     imageC is the image on the centers
     vectN is the E n+th node vector being iterated
     vectC is the B n+th center vector being iterated */

  /* set image for equation on nodes */
  setBC_NodesImage(vct, imageNX, imageNY, imageNZ, vectNX, vectNY, vectNZ, RGBC_Info_Active, RG_numBCMessages_Active, nxn, nyn, nzn);
 
  /*if (!MLMD_fixBCenters){
    if (!vct->getCartesian_rank()){
      cout << "set MLMD_fixBCenters= true for ECSIM BC ---" << endl
	   << "Aborting now ..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    }*/

  /* set image for equations on centers */
  //setBC_NodesImage(vct, imageCX, imageCY, imageCZ, vectCX, vectCY, vectCZ, RGBC_Info_fix3B, RG_numBCMessages_fix3B, nxn, nyn, nzn);

  return;
}

void EMfields3D::BoundaryConditionsECSIMSource(double *** vectNX, double *** vectNY, double *** vectNZ, double *** vectCX, double *** vectCY, double *** vectCZ, VirtualTopology3D *vct, Grid *grid){

  

  /* i set here E N BC */
  setBC_Nodes(vct, vectNX, vectNY, vectNZ, Exth_Active_BC, Eyth_Active_BC, Ezth_Active_BC, RGBC_Info_Active, RG_numBCMessages_Active);

  /* i set here B C BC*/
  /*if (!MLMD_fixBCenters){
    if (!vct->getCartesian_rank()){
      cout << "set MLMD_fixBCenters= true for ECSIM BC ---" << endl
	   << "Aborting now ..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    }*/


  //setBC_Nodes(vct, vectCX, vectCY, vectCZ, Bxc_fix3B_BC, Byc_fix3B_BC, Bzc_fix3B_BC, RGBC_Info_fix3B, RG_numBCMessages_fix3B);

  return;
}

void EMfields3D::centers2nodesB(VirtualTopology3D * vct, Grid *grid, Collective * col){

  cout << "Interpolating B centers to nodes " << endl;
  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);

  /* at the moment, do not worry about ghost nodes -
     first active is calculated, not sure if B ghost cell is OK */
  communicateNode(nxn, nyn, nzn, Bxn, vct);
  communicateNode(nxn, nyn, nzn, Byn, vct);
  communicateNode(nxn, nyn, nzn, Bzn, vct);
}

void EMfields3D::divECorrection_AllFaces(double *** FX, double *** FY, double *** FZ, Grid * grid, VirtualTopology3D * vct){


  cout << "divECorrection_AllFaces " << endl;

  if (Case!= "LightWave"){
    cout << "You have asked for divE correction, but i still have to add the particle part ..." << endl
	 << "Exiting now ... " << endl;
    return;
  }
  if (vct->getXLEN()!= 1 or vct->getYLEN()!=1 or vct->getZLEN()!=1){
    cout << "Not sure divECorrection_Face will work" << endl;
  }

  // FX, FY, FZ are the vector for which i want to apply the divE correction 

  // declare only once, instantiate differently all the times 
  
  int nxn_RED, nyn_RED, nzn_RED;
  double *** FX_RED; double *** FY_RED; double *** FZ_RED;
  // DIR: 0 -> X; 1 -> Y; 2 -> Z
  int DIR;
  // SIDE: 0: LEFT; 1: RIGHT
  int SIDE;

  // FACE X LEFT

  if (vct->getCommField_XLeft()!= MPI_COMM_NULL){
    nxn_RED=4; nyn_RED= nyn; nzn_RED=nzn;
    DIR=0; SIDE=0;

    /* test */
    double ** divE_Pre= newArr2(double, nyc, nzc);
    grid->divN2C_XSide(divE_Pre, FX, FY, FZ, 1);
    double maxdivE_Pre=0.0;
    int YCoordMAX_Pre=1000;
    int ZCoordMAX_Pre=1000;
    for (int j=1; j< nyc-1; j++)
      for (int k=1; k< nzc-1; k++){
	if (divE_Pre[j][k]> maxdivE_Pre){
	  maxdivE_Pre= divE_Pre[j][k];
	  YCoordMAX_Pre=j;
	  ZCoordMAX_Pre=k;
	}
      }
    /* end test */

    cout << "nxn: " << nxn << ", nyn: " << nyn << ", nzn: " << nzn << endl;
    cout << "nxn_RED: " << nxn_RED << ", nyn_RED: " << nyn_RED << ", nzn_RED: " << nzn_RED << endl;
    FX_RED= newArr3(double, nxn_RED, nyn_RED, nzn_RED);
    FY_RED= newArr3(double, nxn_RED, nyn_RED, nzn_RED);
    FZ_RED= newArr3(double, nxn_RED, nyn_RED, nzn_RED);
    
    int I=1;

    for (int i=0; i<nxn_RED; i++)
      for (int j=0; j<nyn_RED; j++)
	for (int k=0; k<nzn_RED; k++){
	  FX_RED[i][j][k]= FX[I][j][k];
	  FY_RED[i][j][k]= FY[I][j][k];
	  FZ_RED[i][j][k]= FZ[I][j][k];
      }
    
    divECorrection_OneFace(FX_RED, FY_RED, FZ_RED, grid, vct, nxn_RED, nyn_RED, nzn_RED, DIR, SIDE);
    

    // i remove boundaries to avoid problems; this will work only in 1D
    for (int j=2; j<nyn-2; j++)
      for (int k=2; k<nzn-2; k++){
	FX[I][j][k]=FX_RED[I][j][k]; // I??
	FY[I][j][k]=FY_RED[I][j][k];
	FZ[I][j][k]=FZ_RED[I][j][k];
      }

    /* test */
    double ** divE_Post= newArr2(double, nyc, nzc);
    grid->divN2C_XSide(divE_Post, FX, FY, FZ, 1);
    double maxdivE_Post=0.0;
    int YCoordMAX_Post=1000;
    int ZCoordMAX_Post=1000;
    for (int j=1; j< nyc-1; j++)
      for (int k=1; k< nzc-1; k++){
	if (divE_Post[j][k]> maxdivE_Post){
	  maxdivE_Post= divE_Post[j][k];
	  YCoordMAX_Post=j;
	  ZCoordMAX_Post=k;
	}
      }
    /* end test */

    cout << "maxdivE_Pre: " << maxdivE_Pre <<" @ [ "<< YCoordMAX_Pre << ", " << ZCoordMAX_Pre << "]; " << endl  << "maxdivE_Post: " << maxdivE_Post <<" @ [ "<< YCoordMAX_Post << ", " << ZCoordMAX_Post << "]; " << endl;
      

    delArr3(FX_RED, nxn_RED, nyn_RED);
    delArr3(FY_RED, nxn_RED, nyn_RED);
    delArr3(FZ_RED, nxn_RED, nyn_RED);

    //delArr2(divE_Pre, nyc);
    delArr2(divE_Post, nyc);
  } // X LEFT
}

void EMfields3D::divECorrection_OneFace(double *** FX, double *** FY, double *** FZ, Grid * grid, VirtualTopology3D * vct, int nxn, int nyn, int nzn, int DIR, int SIDE){

  cout << "INSIDE divECorrection_OneFace: nxn: " << nxn << ", nyn: " << nyn << ", nzn: " << nzn << endl;

  // redefining nxc, nyc, nzc
  int nxc= nxn-1; int nyc= nyn-1; int nzc= nzn-1;
  
  /* these I have redifined here */
  double ***PHI = newArr3(double, nxc, nyc, nzc);
  double ***tempC = newArr3(double, nxc, nyc, nzc);
  /* end these I have redifined here */

  double ***divE = newArr3(double, nxc, nyc, nzc);
  double ***gradPHIX = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIY = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIZ = newArr3(double, nxn, nyn, nzn);

  double *xkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  double *bkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  // set to zero all the stuff 
  eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, divE, nxc, nyc, nzc);
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, gradPHIX, nxn, nyn, nzn);
  eqValue(0.0, gradPHIY, nxn, nyn, nzn);
  eqValue(0.0, gradPHIZ, nxn, nyn, nzn);
  // Adjust E calculating laplacian(PHI) = div(E) -4*PI*rho DIVERGENCE CLEANING
  if (PoissonCorrection_RGFace) {
    if (vct->getCartesian_rank() == 0)
      cout << "*** G" << numGrid <<": DIVERGENCE CLEANING FOR RG FACE***" << endl;
    //grid->divN2C(divE, Ex, Ey, Ez);
    //!!!
    // only central cell in the reduced direction is OK
    grid->divN2C(divE, FX, FY, FZ, nxc, nyc, nzc);
    /* AT THE MOMENT, I REMOVE PARTICLE FEEDBACK
    scale(tempC, rhoc, -FourPI, nxc, nyc, nzc);
    sum(divE, tempC, nxc, nyc, nzc); */
    // move to krylov space 
    /* BC on the source: first active node, df =0 */
    // this only for the X left
    for (int i=0; i< nxc; i++)
      for (int j=0; j< nyc; j++)
	for (int k=0; k< nzc; k++){
	  if (j==1 or k==1 or j==nyc-2 or k==nzc-2){
	    divE[i][j][k]=0.0;
	  }
	}

    phys2solver(bkrylovPoisson, divE, nxc, nyc, nzc);

    MPI_Comm COMM;
    
    if (DIR==0 and SIDE==0) COMM= vct->getCommField_XLeft();
    else{
      cout << "Screwed up the communicator, aborting " << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    GMRES_FACE(&Field::PoissonImage_2D, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 20, 200, GMREStol, grid, vct, this, COMM, nxc, nyc, nzc);

    solver2phys(PHI, xkrylovPoisson, nxc, nyc, nzc);
    // to put back communicateCenterBC(nxc, nyc, nzc, PHI, 2, 2, 2, 2, 2, 2, vct);
    // calculate the gradient
    grid->gradC2N(gradPHIX, gradPHIY, gradPHIZ, PHI, nxn, nyn, nzn);
    // sub
    sub(FX, gradPHIX, nxn, nyn, nzn);
    sub(FY, gradPHIY, nxn, nyn, nzn);
    sub(FZ, gradPHIZ, nxn, nyn, nzn);

  }                             // end of divergence cleaning 

  delete[]xkrylovPoisson;
  delete[]bkrylovPoisson;
  delArr3(divE, nxc, nyc);
  delArr3(gradPHIX, nxn, nyn);
  delArr3(gradPHIY, nxn, nyn);
  delArr3(gradPHIZ, nxn, nyn);

  /* delete the redifined */
  delArr3(PHI, nxc, nyc);
  delArr3(tempC, nxc, nyc);
  /* end delete the redifined */
  
}
/*void EMfields3D::divECorrection_Face(double *** FX, double *** FY, double *** FZ, Grid * grid, VirtualTopology3D * vct){

  // at the moment i am writing this as a test for LW, so i won't put rho -
    // exit if I am not doing LW 
  if (Case!= "LightWave"){
    cout << "You have asked for divE correction, but i still have to add the particle part ..." << endl
	 << "Exiting now ... " << endl;
  }

  // may need a communicateNodes here

  // FX, FY, FZ are the vector for which i want to apply the divE correction 

  // create vectors for E n+theta on the faces
  //   I create them for all, used them only if you are on a face 
  //   - contain only the active nodes 

 
  // declare only once, instantiate differently all the times 
  double ** divE_F;
  double ** PHI_F;

  cout << "for now, this works only with one core" << endl;
  // the GMRES per face
  if (vct->getCommField_XLeft() != MPI_COMM_NULL){
    cout << "Grid " << numGrid << " local rank " << vct->getCartesianRank() <<": i am on field boundary and I am about to do divE correction" << endl;

    double * xkP_Xleft= new double[(nyc - 2) * (nzc - 2)];
    double * bkP_Xleft= new double[(nyc - 2) * (nzc - 2)];

    divE_F=  newArr2(double, nyn, nzn);
    PHI_F= newArr2(double, nyc, nzc);
      
    eqValue(0.0, xkP_Xleft, (nyc - 2) * (nzc - 2));
    
    grid->divN2C_XSide(divE_F, FX, FY, FZ, 1);
    // here i should add the rhoc part 
    
    phys2solver(bkP_Xleft, divE_F, nyc, nzc);

    GMRES(&Field::PoissonImage_2D, xkP_Xleft, (nyc - 2) * (nzc - 2), bkP_Xleft, 20, 200, GMREStol, grid, vct, this);

    solver2phys(PHI_F, xkP_Xleft, nyc, nzc);
    // i won't touch the ghost - arrange communicate somehow 
    // communicateCenter(nxc, nyc, nzc, PHI, vct);

    grid->gradC2N_XSide(gradPHIX_F, gradPHIY_F, gradPHIZ_F, PHI_F);

    int I=1;
    //  still to take care of many cores 
    for (int j=2; j< nyn-1; j++)
      for (int k=2; k< nzn-1; k++){
	FX[I][j][k]-= gradPHIX_F[j][k];
	FY[I][j][k]-= gradPHIY_F[j][k];
	FZ[I][j][k]-= gradPHIZ_F[j][k];
      }
	
    

    delete[]xkP_Xleft; delete[]bkP_Xleft;
    delArr2(divE_F, nyn);
    delArr2(gradPHIX_F, nyn); delArr2(gradPHIY_F, nyn); delArr2(gradPHIZ_F, nyn);
  }
  

}

void EMfields3D::PoissonImage_2D_X(double *image, double *vector, Grid * grid, VirtualTopology3D * vct, int I) {
  // allocate 2 three dimensional service vectors
  double **temp = newArr2(double,  nyc, nzc);
  double **im = newArr2(double,  nyc, nzc);
  eqValue(0.0, image, (nyc - 2) * (nzc - 2));
  eqValue(0.0, temp,  nyc, nzc);
  eqValue(0.0, im,  nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, nxc, nyc, nzc);
  // calculate the laplacian
  grid->lapC2Cpoisson_XSide(im, temp, vct, I);
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);
  // deallocate temporary array and objects
  delArr2(temp, nyc);
  delArr2(im, nyc);
  }*/

/* GEM challenge setting, with the X point at Lx/4 */
void EMfields3D::initGEM_Shifted(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;


  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDx_mlmd(0) ;
  const double coarsedy= grid->getDy_mlmd(0) ;
  const double coarsedz= grid->getDz_mlmd(0) ;

  // this local
  double Lx= col->getLx_mlmd(0);
  double Ly= col->getLy_mlmd(0);
  // end local

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;


  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge (MLMD-ready) with Pertubation located at Lx/4, Ly/2 (global)" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
	cout << "rho species " << i << " = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
	for (int k = 0; k < nzn; k++) {

	  globalx= grid->getXN(i, j, k) + grid->getOx_SW();
	  globaly= grid->getYN(i, j, k) + grid->getOy_SW();
	  globalz= grid->getZN(i, j, k) + grid->getOz_SW();
         
	  double xpert;
	  double ypert;

	  // initialize the density for species
	  for (int is = 0; is < ns; is++) {
	    if (DriftSpecies[is])
	      rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((globaly - Ly / 2) / delta) * cosh((globaly - Ly / 2) / delta)))) / FourPI;
	    else
	      rhons[is][i][j][k] = rhoINIT[is] / FourPI;
	  }
	  // electric field
	  Ex[i][j][k] = 0.0;
	  Ey[i][j][k] = 0.0;
	  Ez[i][j][k] = 0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x * tanh((globaly - Ly / 2) / delta);
	  // add the initial GEM perturbation
	  Byn[i][j][k] = B0y;  
	  // add the initial X perturbation
	  xpert = globalx - Lx / 4; /* THIS IS THE ONLY DIFFERECE W.R.T. INITGEM */
	  ypert = globaly - Ly / 2;
	  exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
	  Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  // guide field
	  Bzn[i][j][k] = B0z;
	}
    // initialize B on centers
   
    communicateNode(nxn, nyn, nzn, Bxn, vct);
    communicateNode(nxn, nyn, nzn, Byn, vct);
    communicateNode(nxn, nyn, nzn, Bzn, vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C_GC(Bxc, Bxn);
    grid->interpN2C_GC(Byc, Byn);
    grid->interpN2C_GC(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenter(nxc, nyc, nzc, Bxc, vct);
    communicateCenter(nxc, nyc, nzc, Byc, vct);
    communicateCenter(nxc, nyc, nzc, Bzc, vct);

    for (int is = 0; is < ns; is++)
      grid->interpN2C_GC(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}
