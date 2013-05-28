/***************************************************************************
  VCtopology3Dparticles.h  -  a 3D Virtual cartesian topology
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

#ifndef VCtopology3Dparticles_H
#define VCtopology3Dparticles_H

#include "mpi.h"
#include "VirtualTopology3D.h"
#include "../utility/Alloc.h"
#include <iostream>



using std::cout;
using std::endl;

/**
 *  
 * Virtual cartesian topology for Particles
 * A virtual topology is a mechanism for naming the processes
 * in a communicator in a way that fits the communication
 * pattern better. Since our processes will communicate mainly
 * with the nearest neighbours after the fashion of a two-dimensional
 * grid, we create a virtual topology to reflect this fact
 * @version 2.0
 */


class VCtopology3Dparticles:public VirtualTopology3D {
public:
  /** constructor: Define topology parameters: dimension, domain decomposition,... */
  VCtopology3Dparticles();
  /** destructor */
  ~VCtopology3Dparticles();
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
  /** get the coordinates in dir direction of process*/
  int getCoordinates(int dir);
  /** get Periodicity condition in dir direction */
  int getPeriods(int dir);
  /** if cVERBOSE == true, print to the screen all the comunication */
  bool getcVERBOSE();

private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;
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
  /** rank may be reordered     */
  int reorder;
  /** arrays for Create_Cart_create  */
  int *divisions;
  /** periodicity */
  int *periods;
  /** coordinates on processors grid */
  int *coordinates;
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

  /** if cVERBOSE == true, print to the screen all the comunication */
  bool cVERBOSE;
};

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
inline VCtopology3Dparticles::VCtopology3Dparticles() {
  // here you have to set the topology
  PERIODICX = true;
  PERIODICY = true;
  PERIODICZ = true;





  // *******************************************
  // *******************************************
  // change these values to change the topology
  XLEN = 4;
  YLEN = 2;
  ZLEN = 4;
  nprocs = XLEN * YLEN * ZLEN;

  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;
  int reorder = 1;
  divisions = new int[3];
  periods = new int[3];
  divisions[0] = XLEN;
  divisions[1] = YLEN;
  divisions[2] = ZLEN;
  periods[0] = PERIODICX;
  periods[1] = PERIODICY;
  periods[2] = PERIODICZ;
  coordinates = new int[3];
  cVERBOSE = false;             // communication verbose ?

}





/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  */
inline void VCtopology3Dparticles::setup_vctopology(MPI_Comm old_comm) {
  // create a matrix with ranks, and neighbours
  MPI_Cart_create(old_comm, 3, divisions, periods, reorder, &CART_COMM);
  if (CART_COMM != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
    MPI_Cart_coords(CART_COMM, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM, XDIR, RIGHT, &xleft_neighbor, &xright_neighbor);
    MPI_Cart_shift(CART_COMM, YDIR, RIGHT, &yleft_neighbor, &yright_neighbor);
    MPI_Cart_shift(CART_COMM, ZDIR, RIGHT, &zleft_neighbor, &zright_neighbor);
  }
  else {
    // EXCEPTION
    cout << "A process is trown away from the new topology. VCtopology3Dparticles.h" << endl;
  }

}
/** destructor */
inline VCtopology3Dparticles::~VCtopology3Dparticles() {
  delete[]periods;
  delete[]divisions;
  delete[]coordinates;
}
/** print topology info */
inline void VCtopology3Dparticles::Print() {
  cout << endl;
  cout << "Particle Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Processors grid: " << XLEN << "x" << YLEN << "x" << ZLEN << endl;
  cout << "Periodicity X: " << periods[0] << endl;
  cout << "Periodicity Y: " << periods[1] << endl;
  cout << "Periodicity z: " << periods[2] << endl;
  cout << endl;
}
/** print cartesian rank of neighbors and coordinate of process */
inline void VCtopology3Dparticles::PrintMapping() {
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << "; Z = " << coordinates[2] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = " << yright_neighbor << "; zLeft = " << zleft_neighbor << "; zRight = " << zright_neighbor << endl;
  cout << endl;
}
/** get XLEN */
inline int VCtopology3Dparticles::getXLEN() {
  return (XLEN);
}
/** get YLEN */
inline int VCtopology3Dparticles::getYLEN() {
  return (YLEN);
}
/** get ZLEN */
inline int VCtopology3Dparticles::getZLEN() {
  return (ZLEN);
}
/** get nprocs */
inline int VCtopology3Dparticles::getNprocs() {
  return (nprocs);
}
/** get periodicity on boundaries - DIRECTION X*/
inline bool VCtopology3Dparticles::getPERIODICX() {
  return (PERIODICX);
}
/** get periodicity on boundaries - DIRECTION Y*/
inline bool VCtopology3Dparticles::getPERIODICY() {
  return (PERIODICY);
}
/** get periodicity on boundaries - DIRECTION Z*/
inline bool VCtopology3Dparticles::getPERIODICZ() {
  return (PERIODICZ);
}
/** get the cartesian rank of the process */
inline int VCtopology3Dparticles::getCartesian_rank() {
  return (cartesian_rank);
}
/** get the cartesian rank of XLEFT neighbor */
inline int VCtopology3Dparticles::getXleft_neighbor() {
  return (xleft_neighbor);
}
/** get the cartesian rank of XRIGHT neighbor */
inline int VCtopology3Dparticles::getXright_neighbor() {
  return (xright_neighbor);
}
/** get the cartesian rank of YLEFT neighbor */
inline int VCtopology3Dparticles::getYleft_neighbor() {
  return (yleft_neighbor);
}
/** get the cartesian rank of YRIGHT neighbor */
inline int VCtopology3Dparticles::getYright_neighbor() {
  return (yright_neighbor);
}
/** get the cartesian rank of ZLEFT neighbor */
inline int VCtopology3Dparticles::getZleft_neighbor() {
  return (zleft_neighbor);
}
/** get the cartesian rank of ZRIGHT neighbor */
inline int VCtopology3Dparticles::getZright_neighbor() {
  return (zright_neighbor);
}
/** if cVERBOSE == true, print to the screen all the comunication */
inline bool VCtopology3Dparticles::getcVERBOSE() {
  return (cVERBOSE);
}
/** get the coordinates in dir direction of process*/
inline int VCtopology3Dparticles::getCoordinates(int dir) {
  return (coordinates[dir]);
}
/** get Periodicity condition in dir direction */
inline int VCtopology3Dparticles::getPeriods(int dir) {
  return (periods[dir]);
}

#endif
