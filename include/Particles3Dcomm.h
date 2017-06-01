/*******************************************************************************************
  Particles3Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef Part3DCOMM_H
#define Part3DCOMM_H

#include "Particles.h"

/* mlmd-related struct: RG P BC */
struct RGPBC_struct {  // when changing this, change MPI_RGPBC_struct_commit also           

  // number of Refined Grid point in the x, y, z direction for which repopulation is needed
  int np_x;
  int np_y;
  int np_z;

  // CG coordinates corresponding to the first point for this GC core           
  double CG_x_first;
  double CG_y_first;
  double CG_z_first;

  /* CG core which sends this set of BCs                                                                                                      
     important: one core per message;     
     the rank is the one on the PARENT-CHILD communicator   
  */
  int CG_core;
  /* RG core involved in the communication;      
     i need it because i plan to have one RG core collecting all the info and   
     sending it to one CG core, to minimise CG-RG communication;   
     the rank is on the PARENT-CHILD communicator*/
  int RG_core;
  // so RG grid knows what she is dealing with  
  // when sending BC, send it back as tag  
  // NB: the MsgID is the order in which that particle RG core builds the msg in init Phase1;  
  int MsgID;
}; // end structure   

// structure to send repopulated particles from the coarse to the refined grids
struct RepP_struct{
  
  // position
  double x;
  double y;
  double z;

  // velocity 
  double u;
  double v;
  double w;

  // charge
  double q;

  // ID
  unsigned long ID;

}; // end structure

// structure to exchange repopulated particles (after split) within RG
struct CRP_struct{
  
  // position
  double x;
  double y;
  double z;

  // velocity 
  double u;
  double v;
  double w;

  // charge
  double q;

  // ID
  unsigned long ID;

  // rank of destination- on the local RG communicator
  int Destination;
}; // end structure


/** Class for particle distribution calculation in 3D */
class c_vDist {
public:
  c_vDist() {;};
  ~c_vDist() {;};
  /** Initialization of particle distribution arrays in 3D */
  void init(int ispec, double vX, double vY, double vZ, int bi, int bj, int bk, double vR, double vFact, Collective * col, Grid * grid);
  void add (double x, double y, double z, double u, double v, double w);
  bool get_doVdist()     { return dovDist3D;};
  double get_dim_i()     { return (vBinEnd_i-vBinBeg_i);};
  double get_dim_j()     { return (vBinEnd_j-vBinBeg_j);};
  double get_dim_k()     { return (vBinEnd_k-vBinBeg_k);};
  double get_vBinBeg_i() { return vBinBeg_i;};
  double get_vBinBeg_j() { return vBinBeg_j;};
  double get_vBinBeg_k() { return vBinBeg_k;};
  double get_vBinEnd_i() { return vBinEnd_i;};
  double get_vBinEnd_j() { return vBinEnd_j;};
  double get_vBinEnd_k() { return vBinEnd_k;};
  double get_dvi()       { return dv_i;};
  double get_dvj()       { return dv_j;};
  double get_dvk()       { return dv_k;};
  int    get_ntotBins()  { return (nBins_i*nBins_j*nBins_k);};
  int    get_nBinsi()    { return nBins_i;};
  int    get_nBinsj()    { return nBins_j;};
  int    get_nBinsk()    { return nBins_k;};
  long   get(int i, int j, int k) { return vDist3D[i][j][k];};
  double get_vDistLoc_x()         { return vDistLoc_x;};
  double get_vDistLoc_y()         { return vDistLoc_y;};
  double get_vDistLoc_z()         { return vDistLoc_z;};
private:
  double           vDistRad;
  double           vDistLoc_x;
  double           vDistLoc_y;
  double           vDistLoc_z;
  double           dv_i;
  double           dv_j;
  double           dv_k;
  bool             dovDist3D;
  int              nBins_i;
  int              nBins_j;
  int              nBins_k;
  unsigned long*** vDist3D;
  double           vBinBeg_i;
  double           vBinEnd_i;
  double           vBinBeg_j;
  double           vBinEnd_j;
  double           vBinBeg_k;
  double           vBinEnd_k;
};

/**
 * 
 * Abstract class for particles of the same species, in a 2D space and 3component velocity with communications methods
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3Dcomm:public Particles {
public:
  /** constructor */
  Particles3Dcomm();
  /** destructor */
  ~Particles3Dcomm();
  /** allocate particles */
  void allocate(int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid);

  /** calculate the weights given the position of particles */
  void calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid);
  /** interpolation method GRID->PARTICLE order 1: CIC */
  void interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes */
  int communicate(VirtualTopology3D * ptVCT);
  /** what it needs to do, is to remove particles from the PRA **/
  int communicateAfterMover(VirtualTopology3D * ptVCT);
  /** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXright(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYright(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZright(double *b_, long long np, VirtualTopology3D * vct);
  /** Delete the a particle from a list(array) and pack the list(array) */
  void del_pack(long long np, long long *nplast);

  /** method to debuild the buffer received */
  /*! mlmd: i also need the communicator */
  //int unbuffer(double *b_);
  int unbuffer(double *b_, MPI_Comm Comm);

  /** resize the receiving buffer */
  void resize_buffers(int new_buffer_size);
  /** a method to compute how many particles are not in the right domain */
  int isMessagingDone(VirtualTopology3D * ptVCT);
  /** calculate the maximum number exiting from this domain */
  int maxNpExiting();
  /** calculate the weights given the position of particles */
  // void calculateWeights(double*** weight, double xp, double yp, double zp,int ix, int iy, int iz, Grid* grid);
  /** get X-position array for all the particles */
  double *getXall() const;
  /** get Y-position array for all the particles */
  double *getYall() const;
  /** get Z-position array for all the particles */
  double *getZall() const;
  /** get u (X-velocity) array for all the particles */
  double *getUall() const;
  /** get v (Y-velocity) array for all the particles */
  double *getVall() const;
  /** get w (Z-velocity) array for all the particles */
  double *getWall() const;
  /** get X-position array for all the particles by reference */
  double *& getXref();
  /** get Y-position array for all the particles by reference */
  double *& getYref();
  /** get Z-position array for all the particles by reference */
  double *& getZref();
  /** get u (X-velocity) array for all the particles by reference */
  double *& getUref();
  /** get v (Y-velocity) array for all the particles by reference */
  double *& getVref();
  /** get w (Z-velocity) array for all the particles by reference */
  double *& getWref();
  /** get q array for all the particles by reference */
  double *& getQref();
  /** get the ID array   */
  unsigned long *getParticleIDall() const;
  /** get X-position of particle with label indexPart */
  double getX(long long indexPart) const;
  /** get Y-position of particle with label indexPart */
  double getY(long long indexPart) const;
  /** get Z-position of particle with label indexPart */
  double getZ(long long indexPart) const;
  /** get u (X-velocity) of particle with label indexPart */
  double getU(long long indexPart) const;
  /** get v (Y-velocity) of particle with label indexPart */
  double getV(long long indexPart) const;
  /** get w (Z-velocity) of particle with label indexPart */
  double getW(long long indexPart) const;
  /** get ID of particle with label indexPart */
  unsigned long getParticleID(long long indexPart) const;
  /**get charge of particle with label indexPart */
  double getQ(long long indexPart) const;
  /** get charge of array for ID particles */
  double *getQall() const;
  /** get the number of particles of this subdomain */
  long long getNOP() const;
  /** return the Kinetic energy */
  /*! mlmd: now i need the communicator also */
  //double getKe();
  double getKe(MPI_Comm Comm);
  /** return the maximum kinetic energy */
  /*! mlmd: now i need the communicator also */
  //double getMaxVelocity();
  double getMaxVelocity(MPI_Comm Comm); 
  /** return energy distribution */
  /*! mlmd: now i need the communicator also */
  //unsigned long *getVelocityDistribution(int nBins, double maxVel);
  unsigned long *getVelocityDistribution(int nBins, double maxVel, MPI_Comm Comm);
  /** return the momentum */
  /*! mlmd: now I need the communicator also */
  //double getP();
  double getP(MPI_Comm Comm);
  /** Print particles info: positions, velocities */
  void Print(VirtualTopology3D * ptVCT) const;
  /** Print the number of particles of this subdomain */
  void PrintNp(VirtualTopology3D * ptVCT) const;
  /** Add distributions in this iteration to the total */
  void Add_vDist3D();
  void Write_vDist3D(string SaveDirName);
  /* initialise the structures for particle BC exchnage */
  void initWeightPBC(Grid * grid, VirtualTopology3D * ptVCT);
  /* some checks after initWeightPBC - comment for production runs */
  void CheckAfterInitWeightPBC(VirtualTopology3D * ptVCT);
  /* Phase 1: RG cores build their side of the map for PBC */
  void initWeightPBC_Phase1(Grid *grid, VirtualTopology3D * vct, RGPBC_struct* RGPBC_Info, int *RG_numPBCMessages);   
  /* commit the structure created for initial CG/RG handshake as MPI_Datatype*/
  void MPI_RGPBC_struct_commit();
  /* commit the structure for the particle CG/RG exchange as MPI_Datatype */
  void MPI_RepP_struct_commit();
  /* commit the structure for repopulated particle exchange within the RG */
  void MPI_CRP_struct_commit();
  // each RG core buid its own PBC communication map
  void Explore3DAndCommit(Grid *grid, int i_s, int i_e, int j_s, int j_e, int k_s, int k_e, int *numMsg, int *MaxSizeMsg, VirtualTopology3D * vct);
  /* add one handshake msg to the list */
  void Assign_RGBC_struct_Values(RGPBC_struct *s, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp, int MsgID_tmp);
  /* build and send particle BC msg -- CG to RG */
  void SendPBC(Grid* grid, VirtualTopology3D * vct);
  /* RG receives PBC msg and acts accordingly */
  void ReceivePBC(Grid* grid, VirtualTopology3D * vct);
  /* prepares P BC msg to the refined grid */
  void buildPBCMsg(Grid* grid, VirtualTopology3D * vct, int ch);
  /* add a particle to the PBC msg */
  void addP(int ch, int n, double x, double y, double z, double u, double v, double w, double q, unsigned long ID, VirtualTopology3D* vct);
  /* resize the buffers responsible for sending repopulated particle from the CG to the RG -- CG side*/
  void resize_CG_MLMD_buffers(VirtualTopology3D * vct);
  /* resize the buffers responsible for sending repopulated particles from teh CG to the RG -- RG side 
     the difference: I do not need to preserve info, only the resize buffers */
  void resize_RG_MLMD_buffers(int NEW_sizePBCMsg);
  /* resize the buffers responsible for communicating repopulated particles -- called by addP_CRP*/
  void resize_CRP_buffers(VirtualTopology3D * vct);
  /* resize the buffers responsible for communicating repopulated particles -- when size is avaialble*/

  void resize_CRP_buffers(VirtualTopology3D * vct, int NewSize);

  /* add a particle for the exchange of the repopulated particles within the R*/
  void addP_CRP(CRP_struct * Vec, int *num, double x, double y, double z, double u, double v, double w, double q, unsigned long ID, int DestinationRank, VirtualTopology3D* vct);			       
  /* add a particle packed into the CRP_struct into the core particle vectors 
     -- lots of checks inside, which can be disabled */
  void unpack_CRP(CRP_struct p, VirtualTopology3D * vct);
  /* ApplyPBC: RG, after receiving PBC, applies them */
  void ApplyPBC(VirtualTopology3D* vct, Grid* grid);
  /* split a particle received form the CG into RG particles */
  void SplitPBC(VirtualTopology3D * vct, Grid* grid, RepP_struct p);
  /* for repopulation method -2 */
  bool RepopulatedParticleHasEnteredRG(Grid* grid, double xTmp, double yTmp, double zTmp, double u, double v, double w);
  /* a barrier on both parent and child side of the particle communicator, to prevent messages for different particle speciesfrom crossing */
  void MPI_Barrier_ParentChild(VirtualTopology3D* vct);
  /* check number of particle sent/ received from each grid */
  void CheckSentReceivedParticles(VirtualTopology3D* vct);
  /* distribute split particles to the right RG core*/
  void communicateRepopulatedParticles(Grid* grid, VirtualTopology3D * vct);

  // moved here from Particles3D.h
  void get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]);
  void get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext);
  void get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez);

protected:
  /** number of species */
  /*! comment: the number of THIS species, not the total number of particle species */
  int ns;
  /** maximum number of particles of this species on this domain. used for memory allocation */
  long long npmax;
  /** number of particles of this species on this domain */
  long long nop;
  /** total number of particles */
  long long np_tot;
  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;
  /** charge to mass ratio */
  double qom;
  /** recon thick */
  double delta;
  /** thermal velocity  - Direction X*/
  double uth;
  /** thermal velocity  - Direction Y*/
  double vth;
  /** thermal velocity  - Direction Z*/
  double wth;
  /** u0 Drift velocity - Direction X */
  double u0;
  /** v0 Drift velocity - Direction Y */
  double v0;
  /** w0 Drift velocity - Direction Z */
  double w0;
  /** Positions arra - X component */
  double *x;
  /** Positions array - Y component */
  double *y;
  /** Positions array - Z component */
  double *z;
  /** Velocities array - X component */
  double *u;
  /** Velocities array - Y component */
  double *v;
  /** Velocities array - Z component */
  double *w;
  /** TrackParticleID */
  bool TrackParticleID;
  /** ParticleID */
  unsigned long *ParticleID;
  /** rank of processor in which particle is created (for ID) */
  int BirthRank[2];
  /** number of variables to be stored in buffer for communication for each particle  */
  int nVar;
  /** Charge array */
  double *q;

  /** Initial charge density */
  double rhoINIT;
  /** Injection charge density */
  double rhoINJECT;

  /** Simulation domain lengths */
  double xstart, xend, ystart, yend, zstart, zend, invVOL;
  /** simulation domain lengths, including ghost cells if bordering MPI_PROC_NULL **/
  double xStart_GC, yStart_GC, zStart_GC, xEnd_GC, yEnd_GC, zEnd_GC;
  /** time step */
  double dt;
  /** Lx = simulation box length - x direction   */
  double Lx;
  /** Ly = simulation box length - y direction   */
  double Ly;
  /** Lz = simulation box length - z direction   */
  double Lz;
  /** grid spacings */
  double dx, dy, dz;
  /** parent grid spacing **/
  double DxP, DyP, DzP;
  /* RF, used in split */
  double RFx, RFy, RFz;
  /** number of grid 
          nodes */
  int nxn, nyn, nzn;
  /** buffers for communication */
  /** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
  int buffer_size;
  /** smaller buffer size */
  int buffer_size_small;
  /** buffer with particles going to the right processor - Direction X */
  double *b_X_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_X_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction X */
  double *b_X_LEFT;
  /** pointer to the buffer for resizing */
  double *b_X_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Y */
  double *b_Y_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Y_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Y */
  double *b_Y_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Y_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Z */
  double *b_Z_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Z_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Z */
  double *b_Z_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Z_LEFT_ptr;

  /** number of particles exiting per cycle*/
  int npExitXright;
  /** number of particles exiting to X-LEFT per cycle*/
  int npExitXleft;
  /** number of particles exiting to Y-RIGHT per cycle*/
  int npExitYright;
  /** number of particles exiting to Y-LEFT per cycle*/
  int npExitYleft;
  /** number of particles exiting to Z-RIGHT per cycle*/
  int npExitZright;
  /** number of particles exiting to Z-LEFT per cycle*/
  int npExitZleft;
  /** total number of particles exiting per cycle */
  int npExit;
  /** number of particles not in the right domain   */
  int rightDomain;


  /** bool for communication verbose */
  bool cVERBOSE;
  /** Boundary condition on particles:
          <ul>
          <li>0 = exit</li>
          <li>1 = perfect mirror</li>
          <li>2 = riemission</li>
          <li>3 = periodic condition </li>
          </ul>
          */
  /** Boundary Condition Particles: FaceXright */
  int bcPfaceXright;
  /** Boundary Condition Particles: FaceXleft */
  int bcPfaceXleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceYright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceYleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceZright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceZleft;
  /** speed of light in vacuum */
  double c;
  /** restart variable for loading particles from restart file */
  int restart;
  /** Number of iteration of the mover*/
  int NiterMover;
  /** velocity of the injection of the particles */
  double Vinj;
  /** removed charge from species */
  double Q_removed;
  /** density of the injection of the particles */
  double Ninj;

  int nvDistLoc;
  c_vDist* vDist;

  // local to this grid
  int XLEN;
  int YLEN;
  int ZLEN;

  /*! mlmd specific variables */
  /*! if true, mlmd related output */
  bool MLMDVerbose;
  /*! rank0 in MPI_COMM_WORLD, for system-wide output */
  int SpokePerson;
  /*! number of the current grid in the mlmd hierarchy */
  int numGrid;
  /* CG SIDE */
  /*! number of children */
  int numChildren;
  /* msg receiving structure 
     [numChildren] [MAX_RG_numPBCMessages] */
  RGPBC_struct ** CG_Info;
  int *CG_numPBCMessages;
  /* struct hosting info on particles to send to the RG 
     [number of children][MaxNumMsg][sizeCG_PBCMsg] */
  RepP_struct *** PCGMsg;
  RepP_struct *** PCGMsg_ptr; // alias to make the resize
  int MaxNumMsg;
  /* number of particles to send to each child core
     [number of children][MaxNumMsg] */
  int ** nopPCGMsg;
  /* END CG SIDE */
  /* BOTH SIDES */
  /* size of the vector -- can be resized up to MAXsizePCMsg -- value set in initWeightPBC 
     broken for CG and RG side, in case one core is both */
  int sizeCG_PBCMsg; // CG side
  int sizeRG_PBCMsg; // RG side
  /* resizing up to here -- value set in initWeightPBC
     same for RG and CG side */
  int MAXsizePBCMsg;
  /* wether to allow the resize of the buffers containing the particles to send CG->RG for repopulation
     false may be heavy on memory consumption                 
     true may be heavy on performance */
  bool AllowPMsgResize;
  /* END BOTH SIDES */
  /* RG SIDE */
  /* number of PBC Msgs to receive, as a child */
  int RG_numPBCMessages;
  /* intermediate, for handshake ops */
  int RG_numPBCMessages_LevelWide;
  /* struct for RG PBC handshake, child side */
  RGPBC_struct * RGPBC_Info;
  /* intermediate, for handshake ops */
  RGPBC_struct * RGPBC_Info_LevelWide;
  /* max size struct */
  int MAX_RG_numPBCMessages;
  /* intermediate, for handshake ops */
  int MAX_RG_numPBCMessages_LevelWide;
  /* struct hosting particles received from the RG (here, they are not split yet)
     [RG_numPBCMessages][sizeRG_PCMsg]*/
  RepP_struct ** PRGMsg;
  /* number of particles received as part of each msg
     [RG_numPBCMessages] */
  int * nopPRGMsg;
  /* debug structure to verify if all PBC msg arrived 
     here, put a 'true' in the line corresponding to a particular msg if it arrived */
  bool * PRGMsgArrived;
  /* general buffer to receive particles BC; before putting them in the proper slot in PRGMsg
     [sizeRG_PCMsg]*/
  RepP_struct* PRGMsg_General;
  /* particles added to the PRA area from MLMD BC */
  int PRA_PAdded;
  /* number of particles in the core at the end of communicate */
  int nop_EndCommunicate;
  /* number of particles in the core AFTER ReceivePBC */
  int nop_AfterReceivePBC;
  /* number of particles to transfer to core HighestRank in communicateRepopulatedParticles*/
  int np_ToCoreH_CRP; 
  /* vectors to exchange repopulated particles within the RG */
  CRP_struct * CRP_ToCoreH;
  CRP_struct * CRP_ToCoreH_ptr; // for resize
  /* Core H counterpart */
  /* [XLEN*YLEN*ZLEN][size_CRP] */
  /* in each row, the particles to send to that core */
  CRP_struct ** H_CRP_Msg;
  CRP_struct ** H_CRP_Msg_ptr; // for resize
  /* [XLEN*YLEN*ZLEN] 
     the number of particles to send to that core */
  int *H_num_CRP_nop;
  CRP_struct * H_CRP_General;
  CRP_struct * H_CRP_General_ptr; // for resize
  /* size of ToCoreH */
  int size_CRP;
  /** cannot resize CRP above this **/
  int MAX_np_CRP;
  /* only in HighestRank: tells you how many msgs to expect while communicating repopulated particles and which cores should participate */
  /* if here 1, this core participates in particle exchange in communicateRepopulatedParticles */
  int * H_CRP_cores;
  /* the number of these msgs that HighestRank should expect */
  int num_H_CRP_cores;
  /* number of CRP particles receive by rank < HighestRank */
  int NH_num_CRP_nop;
  /* particles deleted at the end of communicate*/
  int PRA_deleted;
  /* END RG SIDE */

  /* Index (NOT coordinates, not number of cells) at which PRA starts / ends 
     set in allocate */
  int PRA_XLeft_Start;
  int PRA_XLeft_End;
  int PRA_XRight_Start;
  int PRA_XRight_End;
  int PRA_YLeft_Start;
  int PRA_YLeft_End;
  int PRA_YRight_Start;
  int PRA_YRight_End;
  int PRA_ZLeft_Start;
  int PRA_ZLeft_End;
  int PRA_ZRight_Start;
  int PRA_ZRight_End;

  /* corresponding coordinates --
     calculated with respect to the entire grid, not this particular core */
  double Coord_XLeft_Start;
  double Coord_XLeft_End;
  double Coord_XRight_Start;
  double Coord_XRight_End;
  double Coord_YLeft_Start;
  double Coord_YLeft_End;
  double Coord_YRight_Start;
  double Coord_YRight_End;
  double Coord_ZLeft_Start;
  double Coord_ZLeft_End;
  double Coord_ZRight_Start;
  double Coord_ZRight_End;

  /* MPI Datatype associated to RGPBC_struct; init in MPI_RGBC_struct_commit  */
  MPI_Datatype MPI_RGPBC_struct;
  MPI_Datatype MPI_RepP_struct;
  MPI_Datatype MPI_CRP_struct;

  // communicators to parent/ children divided for particle species 
  MPI_Comm CommToParent_P;
  MPI_Comm *CommToChild_P;

  int *Rank_CommToChildren_P;
  int Rank_CommToParent_P;

  // communicators involving only the cores in the grid involved in particle BC communication
  MPI_Comm COMM_RG_PBCSubset_P;
  /* the CG core (rank in the PG communicator) exchanging PBC info that will communicate with RG cores if communication with PBC function has to be initiated e.g. for the resize of the particle buffers */
  /* as saved by CG */
  /* [numChildren] */
  int *CGSide_CGLeader_PBCSubset;
  /* as saved by RG */
  int RGSide_CGLeader_PBCSubset;
  /* [numChildren] */
  MPI_Comm *COMM_CG_PBCSubset_P;

  /* number of RG cores which communicate with this CG*/
  /* [numChildren]; instantiate only if you- as a CG core- are the spokeperson to a RG 
     only the entry relative to that grid has a valid number */
  int *numRcv;
  /* [numChildren][local variable > the biggest of numRcv] 
     instantiate only if you- as a CG core- are the spokeperson to a RG 
     only the entry relative to that grid has a valid number*/
  int **RcvList;
  /*! end mlmd specific variables */

  // i keep here a copy of Ex, Ey, Ez, Bxn, Byn, Bzn, Bx_ext, By_ext, Bz_ext for the mover but also for moving repopulated particles
  double ***Ex;
  double ***Ey;
  double ***Ez;
  double ***Bx;
  double ***By;
  double ***Bz;
                                                                                                                        
  double ***Bx_ext;
  double ***By_ext;
  double ***Bz_ext;

  double Fext;
  

};

#endif
