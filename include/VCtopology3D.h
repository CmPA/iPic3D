/***************************************************************************
  VCtopology3D.h  -  a 3D Virtual cartesian topology
  A virtual topology is a mechanism for naming the processes
  in a communicator in a way that fits the communication
  pattern better. Since our processes will communicate mainly
  with the nearest neighbours after the fashion of a two-dimensional
  grid, we create a virtual topology to reflect this fact
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL Luveun
developers           : Stefano Markidis, Giovanni Lapenta
***************************************************************************/

/*! the mlmd communicator system:
  MPI_COMM_WORLD: non cartesian communicator for the entire system 
  MPI_COMM_GRID: non cartesian communicator for the grid
  MPI_COMM: cartesian communicator for the grid */


#ifndef VCtopology3D_H
#define VCtopology3D_H

#include "MPIdata.h"
#include "VirtualTopology3D.h"
#include "Collective.h"
#include "Alloc.h"
#include <iostream>

using std::cout;
using std::endl;

/**
 *  
 * Virtual cartesian topology
 * A virtual topology is a mechanism for naming the processes
 * in a communicator in a way that fits the communication
 * pattern better. Since our processes will communicate mainly
 * with the nearest neighbours after the fashion of a two-dimensional
 * grid, we create a virtual topology to reflect this fact
 * @version 2.0
 */


class VCtopology3D:public VirtualTopology3D {
public:
  /** constructor: Define topology parameters: dimension, domain decomposition,... */
  VCtopology3D(Collective *col);
  /** destructor */
  ~VCtopology3D();
  /** pre-mlmd: Find the neighbors in the new communicator 
      void setup_vctopology(MPI_Comm comm_old);
      mlmd: build all communicators here */
  void setup_vctopology(MPI_Comm comm_old, Collective *col);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();
  /** values local to the grid
  /** get XLEN */
  int getXLEN();
  /** get YLEN */
  int getYLEN();
  /** get ZLEN */
  int getZLEN();
  /** end values local to the grid **/
  /* mlmd: access XLEN, YLEN, ZLEN @ grid level */
  int getXLEN(int N);
  int getYLEN(int N);
  int getZLEN(int N);
  /** get nprocs 
      mlmd: @ grid level */
  int getNprocs() {return(nprocs);};
  /** get periodicity on boundaries - DIRECTION X*/
  bool getPERIODICX() {return(PERIODICX);};
  /** get periodicity on boundaries - DIRECTION Y*/
  bool getPERIODICY() {return(PERIODICY);};
  /** get periodicity on boundaries - DIRECTION Z*/
  bool getPERIODICZ() {return(PERIODICZ);};
  /** get periodicity on boundaries - DIRECTION X*/
  bool getPERIODICX_P() {return(PERIODICX_P);};
  /** get periodicity on boundaries - DIRECTION Y*/
  bool getPERIODICY_P() {return(PERIODICY_P);};
  /** get periodicity on boundaries - DIRECTION Z*/
  bool getPERIODICZ_P() {return(PERIODICZ_P);};
  /** get the cartesian rank of the process */
  int getCartesian_rank() {return(cartesian_rank);};
  /** get the total number of process */
  int getNproc() {return(nproc);};
  /** get the cartesian rank of XLEFT neighbor */
  int getXleft_neighbor() {return(xleft_neighbor);};
  /** get the cartesian rank of XRIGHT neighbor */
  int getXright_neighbor() {return(xright_neighbor);};
  /** get the cartesian rank of YLEFT neighbor */
  int getYleft_neighbor() {return(yleft_neighbor);};
  /** get the cartesian rank of YRIGHT neighbor */
  int getYright_neighbor() {return(yright_neighbor);};
  /** get the cartesian rank of ZLEFT neighbor */
  int getZleft_neighbor() {return(zleft_neighbor);};
  /** get the cartesian rank of ZRIGHT neighbor */
  int getZright_neighbor() {return(zright_neighbor);};
  /** get the cartesian rank of XLEFT neighbor */
  int getXleft_neighbor_P() {return(xleft_neighbor_P);};
  /** get the cartesian rank of XRIGHT neighbor */
  int getXright_neighbor_P() {return(xright_neighbor_P);};
  /** get the cartesian rank of YLEFT neighbor */
  int getYleft_neighbor_P() {return(yleft_neighbor_P);};
  /** get the cartesian rank of YRIGHT neighbor */
  int getYright_neighbor_P() {return(yright_neighbor_P);};
  /** get the cartesian rank of ZLEFT neighbor */
  int getZleft_neighbor_P() {return(zleft_neighbor_P);};
  /** get the cartesian rank of ZRIGHT neighbor */
  int getZright_neighbor_P() {return(zright_neighbor_P);};
  /** get the coordinates in dir direction of process*/
  int getCoordinates(int dir) {return(coordinates[dir]);};
  /** get the coordinates of process*/
  int *getCoordinates() {return(coordinates);};
  /** get the number of procs in each direction*/
  int *getDivisions() {return(divisions);};
  /** get Periodicity condition in dir direction */
  int getPeriods(int dir) {return(periods[dir]);};
  /** if cVERBOSE == true, print to the screen all the comunication */
  bool getcVERBOSE() {return(cVERBOSE);};
  /** get the MPI communicator */
  MPI_Comm getComm() {return(CART_COMM);};    /*! pre-mlmd: cartesian communicator for the entire system
						mlmd: cartesian communicator per grid level */
  /*! specific MLMD functions */
  /*! mlmd gets */
  /*! returns the non cartesian communicator per grid;
    to have the cartesian comm per grid, use getComm()*/
  MPI_Comm getCommGrid() {return(MPI_COMM_GRID);}; 
  /*! return the number of the current grid in the mlmd hierarchy 
   NB: the number is not necessarely the same as the level */
  int getNumGrid() {return(numGrid); };
  /*! return the rank of the process in the system-wide communicator, MPI_COMM_WORLD */
  int getSystemWide_rank() {return(systemWide_rank); };
  /** mlmd: get rank on the communicator to parent **/
  int getRank_CommToParent() {return rank_CommToParent;}
  /** mlmd: get rank on the communicator to parent FOR PARTICLES **/
  int getRank_CommToParent_P(int is) {return rank_CommToParent_P[is];}
  /*! return the communicator to the parent; it's MPI_COMM_NULL for the coarse grid */
  MPI_Comm getCommToParent() {return CommToParent; }
  MPI_Comm getCommToParent_BCGhost() {return CommToParent_BCGhost; }
  MPI_Comm getCommToParent_BCBuffer() {return CommToParent_BCBuffer; }
  MPI_Comm getCommToParent_Proj() {return CommToParent_Proj; }
  /*! return the communicator to the parent FOR PARTICLES; it's MPI_COMM_NULL for the coarse grid or if you chose other BC's*/ 
  MPI_Comm getCommToParent_P(int is) {return CommToParent_P[is]; }
  /* returns the number of children in the mlmd hierarchy */
  int getNumChildren() {return (numChildren);}
  /* returns the n-th communicator to child form CommToChildren */
  MPI_Comm getCommToChild(int n) {return CommToChildren[n];}
  MPI_Comm getCommToChild_BCGhost(int n) {return CommToChildren_BCGhost[n];}
  MPI_Comm getCommToChild_BCBuffer(int n) {return CommToChildren_BCBuffer[n];}
  MPI_Comm getCommToChild_Proj(int n) {return CommToChildren_Proj[n];}
  /* return the rank as a parent of the n-th child */
  int getRank_CommToChildren(int nc){return rank_CommToChildren[nc];}
  /* returns the n-th communicator to child form CommToChildren FOR PARTICLES*/
  MPI_Comm getCommToChild_P(int n, int is) {return CommToChildren_P[n][is];}
  /* return the rank as a parent of the n-th child FOR PARTICLES*/
  int getRank_CommToChildren_P(int nc, int is){return rank_CommToChildren_P[nc][is];}

  /*! mlmd test functions */
  /*! tries some basic communication on parent-child inter-communicators and communicators */
  void testMlmdCommunicators();
  /*! end mlmd test functions */

  /* return the number - not the level - of the parent grid */
  int getParentGridNum(){return parentGrid;}

  /* returns the number - not the level or the order in the children vector - of the child grid n */
  int getChildGridNum(int n){return childrenGrid[n];}

  /* return the max number of cores used of a single grid */
  int getMaxGridCoreN() {return MaxGridCoreN;}

  /* return the values of the cartesian coordinate lookup table 
     NB: N is the rank number in the inter-communicator (CommToParent or CommToChildren)*/
  // for fields
  int getXcoord_CommToParent(int N){return Xcoord_CommToParent[N];}
  int getYcoord_CommToParent(int N){return Ycoord_CommToParent[N];}
  int getZcoord_CommToParent(int N){return Zcoord_CommToParent[N];}

  // for particles
  int getXcoord_CommToParent_P(int N, int is){return Xcoord_CommToParent_P[N][is];}
  int getYcoord_CommToParent_P(int N, int is){return Ycoord_CommToParent_P[N][is];}
  int getZcoord_CommToParent_P(int N, int is){return Zcoord_CommToParent_P[N][is];}
  /*! end specific MLMD functions */

private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;    /*! pre-mlmd: field cartesian communicator for the entire system;
                           mlmd: field cartesian communicator per grid */
  /** New communicator with virtual cartesian topology for Particles*/
  MPI_Comm CART_COMM_P;  /*! pre-mlmd: particle cartesian communicator for the entire system;
			   mlmd: particle cartesian communicator per grid */
  /** MPI status during sending and receiving communication */
  MPI_Status status;
  /** Direction X for shift MPI_Cart_Shift*/
  int XDIR;
  /** Direction Y for shift MPI_Cart_Shift*/
  int YDIR;
  /** Direction Z for shift MPI_Cart_Shift*/
  int ZDIR;
  /** RIGHT = +    upwards   shift */
  int RIGHT;
  /** LEFT  = -    downwards shift */
  int LEFT;
  /** dimension of virtual topology */
  int PROCDIM;
  /** local to this grid, structure for all the levels is also accessible **/ 
  /** number of subdomains - Direction X */
  int XLEN;
  /** number of subdomains - Direction Y */
  int YLEN;
  /** number of subdomains - Direction Z */
  int ZLEN;
  /** nprocs = number of processors
   mlmd: on the local grid*/
  int nprocs;
  /** periodicity on boundaries - DIRECTION X*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICX;
  /** periodicity on boundaries - DIRECTION Y*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICY;
  /** periodicity on boundaries - DIRECTION Z*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICZ;
  /** periodicity on boundaries - DIRECTION X*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICX_P;
  /** periodicity on boundaries - DIRECTION Y*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICY_P;
  /** periodicity on boundaries - DIRECTION Z*/
  /** NB: when mlmd, this is the periodicity of the LOCAL grod **/
  bool PERIODICZ_P;
  /** rank may be reordered     */
  int reorder;
  /** arrays for Create_Cart_create  */
  int divisions[3];
  /** periodicity */
  int periods[3];
  int periods_P[3];
  /** coordinates on processors grid */
  int coordinates[3];
  /** cartesian rank 
   mlmd: cartesian rank on the grid cartesian communicator*/
  int cartesian_rank;
  /** Number of processors requested at running time */
  int nproc;
  /** cartesian rank of XLEFT neighbor */
  int xleft_neighbor;
  /** cartesian rank of XRIGHT neighbor */
  int xright_neighbor;
  /** cartesian rank of YLEFT neighbor */
  int yleft_neighbor;
  /** cartesian rank of YRIGHT neighbor */
  int yright_neighbor;
  /** cartesian rank of ZRIGHT neighbor */
  int zleft_neighbor;
  /** cartesian rank of ZLEFT neighbor */
  int zright_neighbor;
  /** cartesian rank of XLEFT neighbor */
  int xleft_neighbor_P;
  /** cartesian rank of XRIGHT neighbor */
  int xright_neighbor_P;
  /** cartesian rank of YLEFT neighbor */
  int yleft_neighbor_P;
  /** cartesian rank of YRIGHT neighbor */
  int yright_neighbor_P;
  /** cartesian rank of ZRIGHT neighbor */
  int zleft_neighbor_P;
  /** cartesian rank of ZLEFT neighbor */
  int zright_neighbor_P;

  /** if cVERBOSE == true, print to the screen all the comunication */
  bool cVERBOSE;
  bool verboseMLMD;

  /*! mlmd specific variables */
  /*! non cartesian communicator at grid level */
  MPI_Comm MPI_COMM_GRID;
  /*! number of mlmd grids */
  int Ngrids;
  /*! number of the current grid
    possible values: 0 -> Ngrids -1 */
  int numGrid;
  /* mlmd specific ranks 
   - see also cartesian_rank*/
  /*! rank in the system-wide communicator - MPI_COMM_WORLD */
  int systemWide_rank;
  /* rank in CommToParent communicator */
  int rank_CommToParent;
  /* rank in CommToChildren communicator */
  int *rank_CommToChildren;
  /* rank in CommToParent communicator FOR PARTICLES*/
  int * rank_CommToParent_P;
  /* rank in CommToChildren communicator FOR PARTICLES*/
  int **rank_CommToChildren_P;

  /* end ranks */
  
  /* topology for the grid system */
  int *XLEN_mlmd;
  int *YLEN_mlmd;
  int *ZLEN_mlmd;
  /* end topology for the grid system */

  /* here, the max number of cores used of a single grid; 
   needed to size communication buffers in EMfields*/
  int MaxGridCoreN;

  /* number of children of the current grid */
  int numChildren;

  /* communicator to parent */
  MPI_Comm CommToParent;
  MPI_Comm CommToParent_BCGhost;
  MPI_Comm CommToParent_BCBuffer;
  MPI_Comm CommToParent_Proj;

  /* communicator to children */
  MPI_Comm *CommToChildren;
  MPI_Comm *CommToChildren_BCGhost;
  MPI_Comm *CommToChildren_BCBuffer;
  MPI_Comm *CommToChildren_Proj;
  /* communicator to parent -- for Particles 
     != from MPI_COMM_NULL if the grid is a child which wants to receive PBC
     but it's != MPI_COMM_NULL also if this particular core does not need to receive msg-
     check RG_numPBCMessages also */
  
  // number of particle species
  int ns;

  // divided by species
  MPI_Comm *CommToParent_P;

  /* communicator to children -- for Particles, divided by species*/
  /* [numChildren][ns]*/
  MPI_Comm **CommToChildren_P;

  /** in case topologies do weird things, the cartesian coordinates (x of Xlen, y of Ylen, z of Zlen)
      of each core in the parent - child communicator is set here **/

  

  /*! TAGS 
    as an experiment, let's start a system of tags to mark the different kind of mlmd communication
    TagsFor*_Parent will be used for communication with the parent
    TagsFor*_Children will be used for communication with the children */
    int TagsForInit_Parent;
    int *TagsForInit_Children;
      
  /*TagsForFieldBC_, TagsForParticleBC_, TagsForProjection ... */
  
    // only mine, not all
    // # of parent grid
    int parentGrid;
    // # of children grid
    int *childrenGrid;

    /* use these variables as a look-up table to find the core you need to communicate with on the inter-grid communicators 
       this is done in case the mapping coordinates <-> rank changes in different implementation
       This mapping is intended to be used in the Init* functions only
       Stop to think if you need to use it somewhere else 
       In each vector/ matrix, one finds the map for the entire communicator
       the index is the core rank in the communicator

       - for particles, also indexed on species
    */

    // this if the grid is a child
    int * Xcoord_CommToParent;
    int * Ycoord_CommToParent;
    int * Zcoord_CommToParent;
    int * numGrid_CommToParent; // this is redundant, just to check      

    // for particles
    int ** Xcoord_CommToParent_P;
    int ** Ycoord_CommToParent_P;
    int ** Zcoord_CommToParent_P;
    int ** numGrid_CommToParent_P; // this is redundant, just to check      

    // this if the grid is a parent
    // rows are the number of children
    int ** Xcoord_CommToChildren;
    int ** Ycoord_CommToChildren;
    int ** Zcoord_CommToChildren;
    int ** numGrid_CommToChildren; // this is redundant, just to check   

    int *** Xcoord_CommToChildren_P;
    int *** Ycoord_CommToChildren_P;
    int *** Zcoord_CommToChildren_P;
    int *** numGrid_CommToChildren_P; // this is redundant, just to check   
    /* end use these variables as a look-up table to find the core you need to communicate to on the inter-grid communicators */ 
};

#endif
