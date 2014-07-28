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
  /** Find the neighbors in the new communicator  */
  void setup_vctopology(MPI_Comm comm_old);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();
  /** get XLEN */
  int getXLEN() {return(XLEN);};
  /** get YLEN */
  int getYLEN() {return(YLEN);};
  /** get ZLEN */
  int getZLEN() {return(ZLEN);};
  /** get nprocs */
  int getNprocs() {return(nprocs);};
  /** get periodicity on boundaries - DIRECTION X*/
  bool getPERIODICX() {return(PERIODICX);};
  /** get periodicity on boundaries - DIRECTION Y*/
  bool getPERIODICY() {return(PERIODICY);};
  /** get periodicity on boundaries - DIRECTION Z*/
  bool getPERIODICZ() {return(PERIODICZ);};
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
  MPI_Comm getComm() {return(CART_COMM);};

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
};

#endif
