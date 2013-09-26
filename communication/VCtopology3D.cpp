#include "mpi.h"
#include "Alloc.h"
#include "VCtopology3D.h"
#include <iostream>

using std::cout;
using std::endl;

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
VCtopology3D::VCtopology3D() {
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
void VCtopology3D::setup_vctopology(MPI_Comm old_comm) {
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
VCtopology3D::~VCtopology3D() {

}
/** print topology info */
void VCtopology3D::Print() {
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
void VCtopology3D::PrintMapping() {
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << "; Z = " << coordinates[2] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = " << yright_neighbor << "; zLeft = " << zleft_neighbor << "; zRight = " << zright_neighbor << endl;
  cout << endl;
}
