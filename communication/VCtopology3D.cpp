
#include "VCtopology3D.h"

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
VCtopology3D::VCtopology3D(Collective *col) {
  // *******************************************
  // *******************************************

  /** mlmd: these set's have to be done AFTER you now your grid number **/
  // here you have to set the topology for the fields
  /*PERIODICX = col->getPERIODICX();
  PERIODICY = col->getPERIODICY();
  PERIODICZ = col->getPERIODICZ();*/
  // here you have to set the topology for the Particles
  /*PERIODICX_P = col->getPERIODICX();
  PERIODICY_P = col->getPERIODICY();
  PERIODICZ_P = col->getPERIODICZ();*/
  /** end mlmd: these set's have to be done AFTER you now your grid number **/

  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;

  reorder = 1;

  /** mlmd: these set's have to be done AFTER you now your grid number **/
  /*periods[0] = PERIODICX;
  periods[1] = PERIODICY;
  periods[2] = PERIODICZ;

  periods_P[0] = PERIODICX_P;
  periods_P[1] = PERIODICY_P;
  periods_P[2] = PERIODICZ_P;*/
  /** mlmd: these set's have to be done AFTER you now your grid number **/

  cVERBOSE = false;             // communication verbose ?

  /*! MLMD specific part */
  Ngrids = col->getNgrids();

  /*! mlmd: inter- and intra- communicators */
  CommToChildren= new MPI_Comm[Ngrids];
  rank_CommToChildren= new int[Ngrids];

  /*! everything initialised to MPI_COMM_NULL */
  for (int ng=0; ng< Ngrids; ng++){
    CommToChildren[ng]= MPI_COMM_NULL;
    rank_CommToChildren[ng]= -1;
  }

  CommToParent= MPI_COMM_NULL;
  rank_CommToParent= -1;
  numChildren= 0;
  /* end mlmd: inter- and intra- communicators */
  
  verboseMLMD= true;
  /*! end MLMD specific part */

  /* TAGS: definitive values assigned in setup_vctopology */
  TagsForInit_Parent= -1;
  TagsForInit_Children= new int[Ngrids];
  for (int t=0; t< Ngrids; t++){
    TagsForInit_Children[t]= -2;
  }

  /** mlmd: XLEN, YLEN, ZLEN different per grids **/
  XLEN_mlmd= new int[Ngrids];
  YLEN_mlmd= new int[Ngrids];
  ZLEN_mlmd= new int[Ngrids];
  
  MaxGridCoreN= 0;
  for (int g=0; g< Ngrids; g++){
    XLEN_mlmd[g]= col->getXLEN_mlmd(g);
    YLEN_mlmd[g]= col->getYLEN_mlmd(g);
    ZLEN_mlmd[g]= col->getZLEN_mlmd(g);

    if (XLEN_mlmd[g]*YLEN_mlmd[g]*ZLEN_mlmd[g] > MaxGridCoreN) 
      MaxGridCoreN= XLEN_mlmd[g]*YLEN_mlmd[g]*ZLEN_mlmd[g];
  }

  int rr;
  MPI_Comm_rank(MPI_COMM_WORLD, &rr);
  if (rr==0){
    cout <<"Inside VCtopology3D, MaxGridCoreN= " << MaxGridCoreN << endl;
  }

}

/** Destructor */
VCtopology3D::~VCtopology3D(){
  
  /* communicator stuff*/
  delete[]CommToChildren;
  delete[]rank_CommToChildren;

  /* TAGS */
  delete[]TagsForInit_Children;
  
  /* topology lists */
  delete[]XLEN_mlmd;
  delete[]YLEN_mlmd;
  delete[]ZLEN_mlmd;

  /* lookup maps for coordinates */
  delete[]Xcoord_CommToParent;
  delete[]Ycoord_CommToParent;
  delete[]Zcoord_CommToParent;
  delete[]numGrid_CommToParent;

  delArr2(Xcoord_CommToChildren, numChildren);
  delArr2(Ycoord_CommToChildren, numChildren);
  delArr2(Zcoord_CommToChildren, numChildren);
  delArr2(numGrid_CommToChildren, numChildren);
  /* end lookup maps for coordinates */
}

/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  
  pre-mlmd:
  inline void VCtopology3D::setup_vctopology(MPI_Comm old_comm) { */

inline void VCtopology3D::setup_vctopology(MPI_Comm old_comm, Collective *col) {
  int size_CW; //on MPI_COMM_WORLD
  
  /* build cartesian communicators for the grids */
  MPI_Comm_size(MPI_COMM_WORLD, &size_CW);
  MPI_Comm_rank(MPI_COMM_WORLD, &systemWide_rank);

  int TotalSize=0;
  for (int g=0; g<Ngrids; g++){
    TotalSize+= XLEN_mlmd[g]*YLEN_mlmd[g]*ZLEN_mlmd[g];  }


  if (size_CW != TotalSize) {
    if (systemWide_rank==0){
      cout << "The number of MPI processes must be " << TotalSize << ", aborting ..." << flush;
      abort();
    }
  }

  numGrid= col->getnumGrid_clt() ;

  /*for (int i=0; i< Ngrids; i++){
    cout << "I am rank " << systemWide_rank << " in VCtopology, I belong to grid " << numGrid << endl;
    }*/

  /*MPI_Barrier(MPI_COMM_WORLD);                                                             
  cout << "exiting now..."<< endl;                                                        
  MPI_Finalize; exit(EXIT_SUCCESS);*/


  // now i knwo my XLEN, YLEN, XLEN and can update the divisions
  XLEN= XLEN_mlmd[numGrid]; 
  YLEN= YLEN_mlmd[numGrid];
  ZLEN= ZLEN_mlmd[numGrid];
  nprocs = XLEN*YLEN*ZLEN;

  divisions[0] = XLEN;
  divisions[1] = YLEN;
  divisions[2] = ZLEN;

  /* periodicity moved here, AFTER i know my numGrid */
  /* this is the periodicity of the local grid */
  PERIODICX = col->getPERIODICX(numGrid);                                                                               
  PERIODICY = col->getPERIODICY(numGrid);
  PERIODICZ = col->getPERIODICZ(numGrid);

  PERIODICX_P = col->getPERIODICX(numGrid);   
  PERIODICY_P = col->getPERIODICY(numGrid); 
  PERIODICZ_P = col->getPERIODICZ(numGrid);

  periods[0] = PERIODICX;                                                                                                      
  periods[1] = PERIODICY;   
  periods[2] = PERIODICZ; 

  periods_P[0] = PERIODICX_P; 
  periods_P[1] = PERIODICY_P; 
  periods_P[2] = PERIODICZ_P;
  /* end periodicity moved here, AFTER i know my numGrid */
  /*! MPI_COMM_GRID is the non cartesian communicator per grid */
  MPI_Comm_split(MPI_COMM_WORLD, numGrid, systemWide_rank, &MPI_COMM_GRID); 
			
  /*! this entire chunk is lifted from the non-mlmd version, with MPI_COMM_GRID instead of old_comm */
  // create a matrix with ranks, and neighbours for fields
  MPI_Cart_create(MPI_COMM_GRID, 3, divisions, periods, reorder, &CART_COMM);
  // create a matrix with ranks, and neighbours for Particles
  MPI_Cart_create(MPI_COMM_GRID, 3, divisions, periods_P, reorder, &CART_COMM_P);
  // field Communicator
  if (CART_COMM != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
    MPI_Comm_size(CART_COMM, &nproc);
    MPI_Cart_coords(CART_COMM, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM, XDIR, RIGHT, &xleft_neighbor, &xright_neighbor);
    MPI_Cart_shift(CART_COMM, YDIR, RIGHT, &yleft_neighbor, &yright_neighbor);
    MPI_Cart_shift(CART_COMM, ZDIR, RIGHT, &zleft_neighbor, &zright_neighbor);
  }
  else {
    // EXCEPTION
    cout << "A process is trown away from the new topology for fields. VCtopology3D.h" << endl;
  }
  // Particles Communicator
  if (CART_COMM_P != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM_P, &cartesian_rank);
    MPI_Comm_size(CART_COMM, &nproc);
    MPI_Cart_coords(CART_COMM_P, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM_P, XDIR, RIGHT, &xleft_neighbor_P, &xright_neighbor_P);
    MPI_Cart_shift(CART_COMM_P, YDIR, RIGHT, &yleft_neighbor_P, &yright_neighbor_P);
    MPI_Cart_shift(CART_COMM_P, ZDIR, RIGHT, &zleft_neighbor_P, &zright_neighbor_P);
  }
  else {
    // EXCEPTION
    cout << "A process is trown away from the new topology for Particles. VCtopology3D.h" << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (systemWide_rank==0) cout << "Before creating children" << endl;

  /*! end: this entire chunk is lifted from the non-mlmd version, with MPI_COMM_GRID instead of old_comm */
  /* end build cartesian communicators for the grids */
  
  /*! build the INTER- and INTRA-communicators between parent and children */
  
  /* tmp variables */

  /*! list of the parent grids; size: Ngrids; parentList[0]=-1 (from Collective::getParentGrid(0)) */
  int *parentList=new int[Ngrids];
  /*! matrix of the children grids; size: [Ngrids][Ngrids]; unused slots = -1 (from collective gets)*/
  int **childrenList= newArr2(int, Ngrids, Ngrids);
  /*! list with the number of children of each grid */
  int * childrenNum=new int[Ngrids];
  /*! list of parent-child INTER-communicators; their are ordered with respect to the numGrid of the child
    it contains the value MPI_COMM_NULL if the specific core is not involved in the communicator 
    if the entry n != MPI_COMM_NULL, if numGrid=n that's the communicator to the parent
    if numGrid!=n that's the communicator to a child */
  MPI_Comm *ChildParentInterComm= new MPI_Comm[Ngrids];
  /* same but for communicators */
  MPI_Comm *ChildParentComm= new MPI_Comm[Ngrids];
  /* end tmp variables for mlmd init */

  // refer to collective for unused slots in childreList and parentList[0]
  for (int ng=0; ng <Ngrids; ng++){
    parentList[ng]= col->getParentGrid(ng);
    childrenNum[ng]= col->getChildrenNum(ng);

    for (int c=0; c< Ngrids; c++)
      childrenList[ng][c]= col->getChildrenGrids(ng, c);
  }
  /* end parent- children list, tmp */
  
                                                                                           
  /*cout << "I am rank " << systemWide_rank << " in VCtopology, I belong to grid " << numGrid << endl;           
  cout << "My parent is " << parentList[numGrid] <<endl;
  cout << "My children are " ;
  for (int j=0;j< childrenNum[numGrid]; j++){ cout << childrenList[numGrid][j] <<" " <<endl; }*/
      
  // i set up my parent and my children
  parentGrid= parentList[numGrid];
  childrenGrid= new int[Ngrids];
  for (int i=0; i<Ngrids; i++) childrenGrid[i]= -1;// to provoke segm fault when i make mistakes
  for (int i=0; i<childrenNum[numGrid]; i++) {childrenGrid[i]= childrenList[numGrid][i];}

  /*MPI_Barrier(MPI_COMM_WORLD);
    if (systemWide_rank==0) cout << "Before creating children comm" << endl;*/

  /*! everything initialised to MPI_COMM_NULL */
  for (int ng=0; ng< Ngrids; ng++){
    ChildParentInterComm[ng]= MPI_COMM_NULL;
    ChildParentComm[ng]= MPI_COMM_NULL;
  }


  // as a child
  if (numGrid>0){ //gridNum=0 is not a child
    int LowestRankParent= col->getLowestRankOfGrid(parentGrid);
    if (cartesian_rank==0){
      cout << "Grid " << numGrid << ", as a child, LRP: " <<LowestRankParent << endl; }
    // as remote leader, I have to put the lowest rank in the common communicator, MPI_COMM_WORLD
    MPI_Intercomm_create(CART_COMM, 0, MPI_COMM_WORLD, LowestRankParent, numGrid, ChildParentInterComm+numGrid); 
    // due to true, the children cores are ranked after the parent cores
    MPI_Intercomm_merge(ChildParentInterComm[numGrid], true, ChildParentComm + numGrid);
  }

  // as a parent
  for (int nc=0; nc< childrenNum[numGrid]; nc++){
    int child= childrenGrid[nc];
    int LowestRankChild= col->getLowestRankOfGrid(child);
    if (cartesian_rank==0){
      cout << "Grid " <<numGrid <<", as a parent, tag " << child << ", LowestRankParent: " << LowestRankChild<<endl;
    }
    MPI_Intercomm_create(CART_COMM, 0, MPI_COMM_WORLD, LowestRankChild, child, ChildParentInterComm + child);
    // due to false, the parent cores are ranked before the children cores
    MPI_Intercomm_merge(ChildParentInterComm[child], false, ChildParentComm + child);
  }  
  
  /*! end tmp communicators array */
  

  MPI_Barrier(MPI_COMM_WORLD);
  if (systemWide_rank==0) cout << "after creating children comm" << endl;

  /* assign value to permanent variables */
  for (int ng=0; ng< Ngrids; ng++){
    // number of children of the current grid
    if (ng==numGrid){
      numChildren= childrenNum[ng];
    }
    
    if (ChildParentInterComm[ng] != MPI_COMM_NULL ){ // ng is either the parent or the child
      if (ng== numGrid){ // ng=numGrid is the child
	CommToParent= ChildParentComm[ng];
	MPI_Comm_rank(CommToParent, &rank_CommToParent);
      } else {
	// ng is the child, numGrid is the parent
	for (int c=0; c<childrenNum[numGrid]; c++ ){ // keep the same children order as in childrenNum
	  if (childrenList[numGrid][c]== ng) {
	    CommToChildren[c]= ChildParentComm[ng];
	    MPI_Comm_rank(CommToChildren[c], rank_CommToChildren + c);
	  }
	}
      } // end of the else
    } // end if != MPI_COMM_NULL
  } // end cycle ng

  /* end assign values to permanent variables */
  
  /* set coordinates of cores in the parent-child communicator,
     to set up the inter-grid operations */

  int my_Xcoord= coordinates[0];
  int my_Ycoord= coordinates[1];
  int my_Zcoord= coordinates[2];


  /*int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int recvcount, MPI_Datatype recvtype,
		    MPI_Comm comm)*/


  if (CommToParent != MPI_COMM_NULL){ // I have a parent 
    Xcoord_CommToParent= new int[2*MaxGridCoreN];
    Ycoord_CommToParent= new int[2*MaxGridCoreN];
    Zcoord_CommToParent= new int[2*MaxGridCoreN];
    numGrid_CommToParent= new int[2*MaxGridCoreN];

    for (int i=0; i< 2*MaxGridCoreN; i++){
      Xcoord_CommToParent[i]=-1;
      Ycoord_CommToParent[i]=-1;
      Zcoord_CommToParent[i]=-1;
      numGrid_CommToParent[i]=-1;
    }

    int size; MPI_Comm_size(CommToParent, &size) ;
    int * vectX= new int[size];
    int * vectY= new int[size];
    int * vectZ= new int[size];
    int * vectG= new int[size];

    MPI_Allgather(&my_Xcoord, 1, MPI_INT, vectX, 1, MPI_INT, CommToParent);
    MPI_Allgather(&my_Ycoord, 1, MPI_INT, vectY, 1, MPI_INT, CommToParent);
    MPI_Allgather(&my_Zcoord, 1, MPI_INT, vectZ, 1, MPI_INT, CommToParent);
    MPI_Allgather(&numGrid, 1, MPI_INT, vectG, 1, MPI_INT, CommToParent);
    //int Start= XLEN_mlmd[parentGrid]*YLEN_mlmd[parentGrid]*ZLEN_mlmd[parentGrid];
    //for (int i=0; i< XLEN*YLEN*ZLEN; i++){ // copy only the child part, considering how big the parent grid is
    for (int i=0; i< size; i++){
      Xcoord_CommToParent[i]= vectX[i];
      Ycoord_CommToParent[i]= vectY[i];
      Zcoord_CommToParent[i]= vectZ[i];
      numGrid_CommToParent[i]= vectG[i];
    }
    delete []vectX;
    delete []vectY;
    delete []vectZ;
    delete []vectG;

  } // end if (CommToParent != MPI_COMM_NULL){ // I have a parent  

  
  if (numChildren > 0){ // I have children
    Xcoord_CommToChildren= newArr2(int, numChildren, 2*MaxGridCoreN);
    Ycoord_CommToChildren= newArr2(int, numChildren, 2*MaxGridCoreN);
    Zcoord_CommToChildren= newArr2(int, numChildren, 2*MaxGridCoreN);
    numGrid_CommToChildren= newArr2(int, numChildren, 2*MaxGridCoreN);
    
    for (int i=0; i< numChildren; i++){
      for (int j=0; j< 2*MaxGridCoreN; j++){
	Xcoord_CommToChildren[i][j]=-1;
	Ycoord_CommToChildren[i][j]=-1;
	Zcoord_CommToChildren[i][j]=-1;
	numGrid_CommToChildren[i][j]=-1;
      }
    }

    for (int nc=0; nc< numChildren; nc++){
      int size; MPI_Comm_size(CommToChildren[nc], &size) ;
      int * vectX= new int[size];
      int * vectY= new int[size];
      int * vectZ= new int[size];
      int * vectG= new int[size];

      MPI_Allgather(&my_Xcoord, 1, MPI_INT, vectX, 1, MPI_INT, CommToChildren[nc]);
      MPI_Allgather(&my_Ycoord, 1, MPI_INT, vectY, 1, MPI_INT, CommToChildren[nc]);
      MPI_Allgather(&my_Zcoord, 1, MPI_INT, vectZ, 1, MPI_INT, CommToChildren[nc]);
      MPI_Allgather(&numGrid, 1, MPI_INT, vectG, 1, MPI_INT, CommToChildren[nc]);

      //for (int i=0; i< XLEN*YLEN*ZLEN; i++){ // copy only the parent part
      for (int i=0; i< size; i++){ // copy only the parent part
	Xcoord_CommToChildren[nc][i]= vectX[i];
	Ycoord_CommToChildren[nc][i]= vectY[i];
	Zcoord_CommToChildren[nc][i]= vectZ[i];
	numGrid_CommToChildren[nc][i]= vectG[i];
      }
    
      delete []vectX;
      delete []vectY;
      delete []vectZ;
      delete []vectG;
    } // end of for (int nc=0; nc< numChildren; nc++){

  } // end if (numChildren > 0){


  /* quick test to check that the same info is shared on the parent and child side,
     check that the outputs are the same 
     NB: for this to work, grid 1 must be a child of grid 0
  */
  /*if (systemWide_rank==XLEN_mlmd[0]*YLEN_mlmd[0]*ZLEN_mlmd[0]){
    int size;
    MPI_Comm_size(CommToParent, &size) ;
    for (int i=0; i< size+2 ;i++){
      cout <<"R" << systemWide_rank <<", i: " << i <<"->["<< Xcoord_CommToParent[i] <<", " <<Ycoord_CommToParent[i] <<", "<<Zcoord_CommToParent[i]<<"]"  <<", numGrid_CommToParent[i]: " << numGrid_CommToParent[i]<<endl;
    }
  }
  
  if (systemWide_rank==0){
    int size; 
    MPI_Comm_size(CommToChildren[0], &size) ;
    for (int i=0; i<size+2; i++){
      cout <<"R" << systemWide_rank <<", i: " << i <<"->["<< Xcoord_CommToChildren[0][i] <<", " <<Ycoord_CommToChildren[0][i] <<", "<<Zcoord_CommToChildren[0][i]<<"]"  <<", numGrid_CommToChildren[0][i]: " << numGrid_CommToChildren[0][i]<<endl;
    }
    }*/
  /* end quick test to check that the same info is shared on the parent and child side */
  
  /*   if (verboseMLMD){
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int i=0; i< size_CW; i++){
      if (systemWide_rank==i){ // this only not to mess up the outputs
	cout <<"Rank in MPI_COMM_WORLD: " << systemWide_rank << endl;
	cout <<"grid number: " << numGrid << endl;
	cout <<"Rank in local communicator: " <<  cartesian_rank << endl;
	cout <<"XLEN=" << XLEN << ", YLEN= " << YLEN << ", ZLEN= " << ZLEN << endl;
	cout <<"Intercomms: " << endl; 
	for (int ng=0; ng< Ngrids; ng++){
	  cout << "# " << ng <<" : " << ((ChildParentInterComm[ng] == MPI_COMM_NULL)? "No" : "Yes") << ": " << ChildParentInterComm[ng] << endl;	
	}
      	
	// now, communicators
	cout <<"comms: " << endl; 
	for (int ng=0; ng< Ngrids; ng++){
	  cout << "# " << ng <<" : " << ((ChildParentComm[ng] == MPI_COMM_NULL)? "No" : "Yes") << ": " << ChildParentComm[ng] << endl;	
	}      
	cout << "comm to parent: " <<((CommToParent == MPI_COMM_NULL)? "No" : "Yes") <<" : " << CommToParent << ", rank: " << rank_CommToParent << endl;
	cout << "comms to children: ";
	for (int c=0; c< numChildren; c++) {
	  cout<< ((CommToChildren[c] == MPI_COMM_NULL)? "No" : "Yes") <<" :  " << CommToChildren[c] <<", rank: " << rank_CommToChildren[c] << "  ";
	}
	cout << endl;
	cout <<"-----------------------------------------------" << endl;
      } // end print process by process
      MPI_Barrier(MPI_COMM_WORLD);
    }
    } // end verboseMLMD*/



  // uncomment for the test: LookUpMap4Coords
  int size;
  MPI_Barrier(MPI_COMM_WORLD);
  cout <<"R" << systemWide_rank <<": numGrid: " << numGrid << " ,rank in local comm: " << cartesian_rank <<endl;
  cout <<"R" << systemWide_rank <<": " << "coords got locally: [" <<coordinates[0] <<"; " <<coordinates[1] <<"; " <<coordinates[2] <<"]" <<endl;
  if (CommToParent!= MPI_COMM_NULL){
    MPI_Comm_size(CommToParent, &size); 
    cout <<"R" << systemWide_rank <<": " << "Rank as a child: " << rank_CommToParent <<"/ " <<  size << endl;
    cout <<"R" << systemWide_rank <<": " << "coords from *coord_CommToParent[rank_CommToParent]: [" << Xcoord_CommToParent[rank_CommToParent] <<"; " << Ycoord_CommToParent[rank_CommToParent] <<"; " << Zcoord_CommToParent[rank_CommToParent] <<"], numGrid_CommToParent[rank_CommToParent]: " << numGrid_CommToParent[rank_CommToParent] <<endl;
    }
  for (int nc=0; nc<numChildren; nc++){
    MPI_Comm_size(CommToChildren[nc], &size);
    cout <<"R" << systemWide_rank <<": " << "Rank as a parent of child " << nc <<": " << rank_CommToChildren[nc] <<"/ " <<size << endl;
    cout <<"R" << systemWide_rank <<": " << "coords from *coord_CommToChildren[nc][rank_CommToChildren[nc]]: [" << Xcoord_CommToChildren[nc][rank_CommToChildren[nc]] <<"; " << Ycoord_CommToChildren[nc][rank_CommToChildren[nc]] <<"; " << Zcoord_CommToChildren[nc][rank_CommToChildren[nc]] <<"], numGrid_CommToChildren[nc][rank_CommToChildren[nc]]: " << numGrid_CommToChildren[nc][rank_CommToChildren[nc]] <<endl;
    }

  // end uncomment for the test: LookUpMap4Coords 
   /* SETUP OF TAGS VALUES */
   // as a child, TagsForInit_Parent is to communicate with your parent: your gridNum
   TagsForInit_Parent= numGrid;
   // as a parent, TagsForInit_Children[c] is to communicate with chukd c: its gridNum
   for(int c=0; c<numChildren; c++){
     TagsForInit_Children[c]= childrenList[numGrid][c];
   }
   /* END SETUP OF TAGS VALUES */

   /* cleaning up tmp arrays*/
   delete[]parentList;
   delete[]childrenList;
   delete[]childrenNum;
   delete[]ChildParentInterComm;
   delete[]ChildParentComm;

}



/** print topology info */
inline void VCtopology3D::Print() {
  cout << endl;
  cout << "Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Processors grid: " << XLEN << "x" << YLEN << "x" << ZLEN << endl;
  cout << "Periodicity Field X: " << periods[0] << endl;
  cout << "Periodicity Field Y: " << periods[1] << endl;
  cout << "Periodicity Field z: " << periods[2] << endl;
  cout << "Periodicity Particles X: " << periods_P[0] << endl;
  cout << "Periodicity Particles Y: " << periods_P[1] << endl;
  cout << "Periodicity Particles z: " << periods_P[2] << endl;
  cout << endl;
}
/** print cartesian rank of neighbors and coordinate of process */
inline void VCtopology3D::PrintMapping() {
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << "; Z = " << coordinates[2] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = " << yright_neighbor << "; zLeft = " << zleft_neighbor << "; zRight = " << zright_neighbor << endl;
  cout << endl;
}


/*! MLMD specific functions */

/* mlmd test functions */
void VCtopology3D::testMlmdCommunicators(){

  //  cout << "killing the test... " << endl;
  //return;

  int *FromParent= new int[2];
  MPI_Status statusS, statusR;
  MPI_Request requestS, requestR;
  int tag;

  MPI_Barrier(MPI_COMM_WORLD);

  // rank 0 of parents send
  for (int c=0; c<numChildren; c++){
    if (rank_CommToChildren[c]==0){ 
      FromParent[0]= numGrid;
      FromParent[1]= c;
      tag= TagsForInit_Children[c];
      MPI_Isend(FromParent, 2, MPI_INT, XLEN*YLEN*ZLEN, tag, CommToChildren[c], &requestS );
      MPI_Wait(&requestS, &status);
      //cout << "I am grid "<< numGrid << " in testMlmdCommunicator "<< endl;
    }
  }

  // rank XLEN*YLEN*ZLEN of children receives
  if (rank_CommToParent == XLEN_mlmd[parentGrid]*YLEN_mlmd[parentGrid]*ZLEN_mlmd[parentGrid]){
    tag= TagsForInit_Parent;
    MPI_Recv(FromParent, 2, MPI_INT, 0, tag, CommToParent, &statusR);
    cout << "testMlmdCommunicators: Grid " << numGrid << ", child " << FromParent[1] << " of ParentGrid " << FromParent[0] << "has received a message on the ParentChild communicator, tag " << tag << endl;
  }

  MPI_Barrier(CART_COMM);
  MPI_Barrier(MPI_COMM_WORLD);
  if (systemWide_rank==0){
    cout << "Preliminary test on communicators passed " << endl;
  }

  
  delete[]FromParent;

  
}

/** get XLEN */
int VCtopology3D::getXLEN() {return(XLEN);};
/** get YLEN */
int VCtopology3D::getYLEN() {return(YLEN);};
/** get ZLEN */
int VCtopology3D::getZLEN() {return(ZLEN);};
/** end values local to the grid **/
/* mlmd: access XLEN, YLEN, ZLEN @ grid level */
int VCtopology3D::getXLEN(int N) {return XLEN_mlmd[N]; }
int VCtopology3D::getYLEN(int N) {return YLEN_mlmd[N]; }
int VCtopology3D::getZLEN(int N) {return ZLEN_mlmd[N]; }


/* end mlmd test functions */
/*! end MLMD specific functions */
