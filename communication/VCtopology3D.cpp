
#include "VCtopology3D.h"

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
VCtopology3D::VCtopology3D(Collective *col) {
  // *******************************************
  // *******************************************
  // change these values to change the topology
  XLEN = col->getXLEN();
  YLEN = col->getYLEN();
  ZLEN = col->getZLEN();
  nprocs = XLEN * YLEN * ZLEN;  // mlmd: @grid level
  // here you have to set the topology for the fields
  PERIODICX = col->getPERIODICX();
  PERIODICY = col->getPERIODICY();
  PERIODICZ = col->getPERIODICZ();
  // here you have to set the topology for the Particles
  PERIODICX_P = col->getPERIODICX();
  PERIODICY_P = col->getPERIODICY();
  PERIODICZ_P = col->getPERIODICZ();
  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;

  reorder = 1;

  divisions[0] = XLEN;
  divisions[1] = YLEN;
  divisions[2] = ZLEN;

  periods[0] = PERIODICX;
  periods[1] = PERIODICY;
  periods[2] = PERIODICZ;

  periods_P[0] = PERIODICX_P;
  periods_P[1] = PERIODICY_P;
  periods_P[2] = PERIODICZ_P;

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
  
}

/** Destructor */
VCtopology3D::~VCtopology3D(){
  
  /* communicator stuff*/
  delete[]CommToChildren;
  delete[]rank_CommToChildren;

  /* TAGS */
  delete[]TagsForInit_Children;
  
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

  if (size_CW != XLEN*YLEN*ZLEN*Ngrids) {
    if (systemWide_rank==0){
      cout << "The number of MPI processes must be XLEN*YLEN*ZLEN*Ngrids= " << XLEN*YLEN*ZLEN*Ngrids << ", aborting ..." << flush;
      abort();
    }
  }

  numGrid = systemWide_rank % Ngrids;  /*! this is the ID of the current grid  */
  
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
  /*! end: this entire chunk is lifted from the non-mlmd version, with MPI_COMM_GRID instead of old_comm */
  /* end build cartesian communicators for the grids */
  
  /*! build the INTER- and INTRA-communicators between parent and children */
  
  /* tmp variables */

  /*! list of the parent grids; size: Ngrids; parentList[0]=-1 (from Collective::getParentGrid(0)) */
  int *parentList;
  /*! matrix of the children grids; size: [Ngrids][Ngrids]; unused slots = -1 (from collective gets)*/
  int **childrenList;
  /*! list with the number of children of each grid */
  int * childrenNum;
  /*! list of parent-child INTER-communicators; their are ordered with respect to the numGrid of the child
    it contains the value MPI_COMM_NULL if the specific core is not involved in the communicator 
    if the entry n != MPI_COMM_NULL, if numGrid=n that's the communicator to the parent
    if numGrid!=n that's the communicator to a child */
  MPI_Comm *ChildParentInterComm;
  /* same but for communicators */
  MPI_Comm *ChildParentComm; 
  /* end tmp variables for mlmd init */


  /* parent- children list, tmp*/
  parentList= new int[Ngrids];
  childrenNum= new int[Ngrids];
  childrenList = newArr2(int, Ngrids, Ngrids);

  // refer to collective for unused slots in childreList and parentList[0]
  for (int ng=0; ng <Ngrids; ng++){
    parentList[ng]= col->getParentGrid(ng);
    childrenNum[ng]= col->getChildrenNum(ng);

    for (int c=0; c< Ngrids; c++)
      childrenList[ng][c]= col->getChildrenGrids(ng, c);
  }
  /* end parent- children list, tmp */
  
  /*! tmp communicators arrays */
  ChildParentInterComm= new MPI_Comm[Ngrids];
  ChildParentComm= new MPI_Comm[Ngrids];

  /*! everything initialised to MPI_COMM_NULL */
  for (int ng=0; ng< Ngrids; ng++){
    ChildParentInterComm[ng]= MPI_COMM_NULL;
    ChildParentComm[ng]= MPI_COMM_NULL;
  }


  // as a child
  if (numGrid>0){ //gridNum=0 is not a child
    /*if (cartesian_rank==0){
      cout << "Grid " << numGrid << ", as a child" << endl; }*/
    // this value for the remote leader depends on the MPI_Comm_Split
    MPI_Intercomm_create(CART_COMM, 0, MPI_COMM_WORLD, parentList[numGrid], numGrid, ChildParentInterComm+numGrid);
    MPI_Intercomm_merge(ChildParentInterComm[numGrid], true, ChildParentComm + numGrid);
  }

  // as a parent
  for (int nc=0; nc< childrenNum[numGrid]; nc++){
    int child= childrenList[numGrid][nc];
    /*if (cartesian_rank==0){
      cout << "Grid " <<numGrid <<", as a parent, tag " << child <<endl;
      }*/
    // this value for the remote leader depends on the MPI_Comm_Split
    MPI_Intercomm_create(CART_COMM, 0, MPI_COMM_WORLD, child, child, ChildParentInterComm + child);
    MPI_Intercomm_merge(ChildParentInterComm[child], false, ChildParentComm + child);
  }  
  
  /*! end tmp communicators array */
  
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
  

   if (verboseMLMD){
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int i=0; i< size_CW; i++){
      if (systemWide_rank==i){ // this only not to mess up the outputs
	cout <<"Rank in MPI_COMM_WORLD: " << systemWide_rank << endl;
	cout <<"grid number: " << numGrid << endl;
	cout <<"Rank in local communicator: " <<  cartesian_rank << endl;
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
   } // end verboseMLMD

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
      cout << "I am grid "<< numGrid << " in testMlmdCommunicator "<< endl;
    }
  }

  // rank XLEN*YLEN*ZLEN of children receives
  if (rank_CommToParent == XLEN*YLEN*ZLEN){
    tag= TagsForInit_Parent;
    MPI_Recv(FromParent, 2, MPI_INT, 0, tag, CommToParent, &statusR);
    cout << "Grid " << numGrid << ", child " << FromParent[1] << " of ParentGrid " << FromParent[0] << "has received a message on the ParentChild communicator, tag " << tag << endl;
  }

  MPI_Barrier(CART_COMM);
  MPI_Barrier(MPI_COMM_WORLD);
  if (systemWide_rank==0){
    cout << "Preliminary test on communicators passed " << endl;
  }

  
  delete[]FromParent;

  
}
/* end mlmd test functions */
/*! end MLMD specific functions */
