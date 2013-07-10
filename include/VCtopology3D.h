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

#ifndef VCtopology3D_H
#define VCtopology3D_H

#include "mpi.h"
#include "VirtualTopology3D.h"
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
  VCtopology3D();
  /** destructor */
  ~VCtopology3D();
  /** Find the neighbors in the new communicator  */
  void setup_vctopology(MPI_Comm comm_old);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();
  /** get XLEN */
  int getXLEN();
  /** get YLEN */
  int getYLEN();
  /** get ZLEN */
  int getZLEN();
  /** get nprocs */
  int getNprocs();
  /** get periodicity on boundaries - DIRECTION X*/
  bool getPERIODICX();
  /** get periodicity on boundaries - DIRECTION Y*/
  bool getPERIODICY();
  /** get periodicity on boundaries - DIRECTION Z*/
  bool getPERIODICZ();
  /** get the cartesian rank of the process */
  int getCartesian_rank();
  /** get the cartesian rank of XLEFT neighbor */
  int getXleft_neighbor();
  /** get the cartesian rank of XRIGHT neighbor */
  int getXright_neighbor();
  /** get the cartesian rank of YLEFT neighbor */
  int getYleft_neighbor();
  /** get the cartesian rank of YRIGHT neighbor */
  int getYright_neighbor();
  /** get the cartesian rank of ZLEFT neighbor */
  int getZleft_neighbor();
  /** get the cartesian rank of ZRIGHT neighbor */
  int getZright_neighbor();
  /** get the cartesian rank of XLEFT neighbor */
  int getXleft_neighbor_P();
  /** get the cartesian rank of XRIGHT neighbor */
  int getXright_neighbor_P();
  /** get the cartesian rank of YLEFT neighbor */
  int getYleft_neighbor_P();
  /** get the cartesian rank of YRIGHT neighbor */
  int getYright_neighbor_P();
  /** get the cartesian rank of ZLEFT neighbor */
  int getZleft_neighbor_P();
  /** get the cartesian rank of ZRIGHT neighbor */
  int getZright_neighbor_P();
  /** get the coordinates in dir direction of process*/
  int getCoordinates(int dir);
  /** get the coordinates of process*/
  int *getCoordinates();
  /** get Periodicity condition in dir direction */
  int getPeriods(int dir);
  /** if cVERBOSE == true, print to the screen all the comunication */
  bool getcVERBOSE();
  /** get the MPI communicator */
  MPI_Comm getComm();

private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;
  /** New communicator with virtual cartesian topology for Particles*/
  MPI_Comm CART_COMM_P;
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
  /** number of subdomains - Direction X */
  int XLEN;
  /** number of subdomains - Direction Y */
  int YLEN;
  /** number of subdomains - Direction Z */
  int ZLEN;
  /** nprocs = number of processors */
  int nprocs;
  /** periodicity on boundaries - DIRECTION X*/
  bool PERIODICX;
  /** periodicity on boundaries - DIRECTION Y*/
  bool PERIODICY;
  /** periodicity on boundaries - DIRECTION Z*/
  bool PERIODICZ;
  /** periodicity on boundaries - DIRECTION X*/
  bool PERIODICX_P;
  /** periodicity on boundaries - DIRECTION Y*/
  bool PERIODICY_P;
  /** periodicity on boundaries - DIRECTION Z*/
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
  /** cartesian rank */
  int cartesian_rank;
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
};

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
inline VCtopology3D::VCtopology3D() {
  // *******************************************
  // *******************************************
  // change these values to change the topology
  XLEN = 2;
  YLEN = 2;
  ZLEN = 1;
  nprocs = XLEN * YLEN * ZLEN;
  // here you have to set the topology for the fields
  PERIODICX = true;
  PERIODICY = false;
  PERIODICZ = true;
  // here you have to set the topology for the Particles
  PERIODICX_P = true;
  PERIODICY_P = false;
  PERIODICZ_P = true;
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

}





/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  */
inline void VCtopology3D::setup_vctopology(MPI_Comm old_comm) {
  // create a matrix with ranks, and neighbours for fields
  MPI_Cart_create(old_comm, 3, divisions, periods, reorder, &CART_COMM);
  // create a matrix with ranks, and neighbours for Particles
  MPI_Cart_create(old_comm, 3, divisions, periods_P, reorder, &CART_COMM_P);
  // field Communicator
  if (CART_COMM != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
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
    MPI_Cart_coords(CART_COMM_P, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM_P, XDIR, RIGHT, &xleft_neighbor_P, &xright_neighbor_P);
    MPI_Cart_shift(CART_COMM_P, YDIR, RIGHT, &yleft_neighbor_P, &yright_neighbor_P);
    MPI_Cart_shift(CART_COMM_P, ZDIR, RIGHT, &zleft_neighbor_P, &zright_neighbor_P);
  }
  else {
    // EXCEPTION
    cout << "A process is trown away from the new topology for Particles. VCtopology3D.h" << endl;
  }

}
/** destructor */
inline VCtopology3D::~VCtopology3D() {

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
/** get XLEN */
inline int VCtopology3D::getXLEN() {
  return (XLEN);
}
/** get YLEN */
inline int VCtopology3D::getYLEN() {
  return (YLEN);
}
/** get ZLEN */
inline int VCtopology3D::getZLEN() {
  return (ZLEN);
}
/** get nprocs */
inline int VCtopology3D::getNprocs() {
  return (nprocs);
}
/** get periodicity on boundaries - DIRECTION X*/
inline bool VCtopology3D::getPERIODICX() {
  return (PERIODICX);
}
/** get periodicity on boundaries - DIRECTION Y*/
inline bool VCtopology3D::getPERIODICY() {
  return (PERIODICY);
}
/** get periodicity on boundaries - DIRECTION Z*/
inline bool VCtopology3D::getPERIODICZ() {
  return (PERIODICZ);
}
/** get the cartesian rank of the process */
inline int VCtopology3D::getCartesian_rank() {
  return (cartesian_rank);
}
/** get the cartesian rank of XLEFT neighbor */
inline int VCtopology3D::getXleft_neighbor() {
  return (xleft_neighbor);
}
/** get the cartesian rank of XRIGHT neighbor */
inline int VCtopology3D::getXright_neighbor() {
  return (xright_neighbor);
}
/** get the cartesian rank of YLEFT neighbor */
inline int VCtopology3D::getYleft_neighbor() {
  return (yleft_neighbor);
}
/** get the cartesian rank of YRIGHT neighbor */
inline int VCtopology3D::getYright_neighbor() {
  return (yright_neighbor);
}
/** get the cartesian rank of ZLEFT neighbor */
inline int VCtopology3D::getZleft_neighbor() {
  return (zleft_neighbor);
}
/** get the cartesian rank of ZRIGHT neighbor */
inline int VCtopology3D::getZright_neighbor() {
  return (zright_neighbor);
}
/** get the cartesian rank of XLEFT neighbor */
inline int VCtopology3D::getXleft_neighbor_P() {
  return (xleft_neighbor_P);
}
/** get the cartesian rank of XRIGHT neighbor */
inline int VCtopology3D::getXright_neighbor_P() {
  return (xright_neighbor_P);
}
/** get the cartesian rank of YLEFT neighbor */
inline int VCtopology3D::getYleft_neighbor_P() {
  return (yleft_neighbor_P);
}
/** get the cartesian rank of YRIGHT neighbor */
inline int VCtopology3D::getYright_neighbor_P() {
  return (yright_neighbor_P);
}
/** get the cartesian rank of ZLEFT neighbor */
inline int VCtopology3D::getZleft_neighbor_P() {
  return (zleft_neighbor_P);
}
/** get the cartesian rank of ZRIGHT neighbor */
inline int VCtopology3D::getZright_neighbor_P() {
  return (zright_neighbor_P);
}
/** if cVERBOSE == true, print to the screen all the comunication */
inline bool VCtopology3D::getcVERBOSE() {
  return (cVERBOSE);
}
/** get the coordinates in dir direction of process*/
inline int VCtopology3D::getCoordinates(int dir) {
  return (coordinates[dir]);
}
/** get the coordinates in dir direction of process*/
inline int *VCtopology3D::getCoordinates() {
  return (coordinates);
}
/** get Periodicity condition in dir direction */
inline int VCtopology3D::getPeriods(int dir) {
  return (periods[dir]);
}
/** Get the MPI communicator */
inline MPI_Comm VCtopology3D::getComm(){
  return (CART_COMM);
}

#endif
