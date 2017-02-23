/***************************************************************************
  MPIdata.h  -  MPI data and methods wrapper
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef MPIDATA_H
#define MPIDATA_H

#include <mpi.h>
#include <iostream>

using std::cout;
using std::endl;

/**
 * MPI Data Structure. This class contains:
 *
 * - rank of process
 * - number of processors
 * - size of communication buffer
 *
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
class MPIdata {
public:
  /** constructor: setup MPI environment */
  MPIdata(int *, char ***);
  /** destructor */
   ~MPIdata();
  /** initialize MPIdata */
  void init(int *, char ***);
  /** close MPI environment */
  void finalize_mpi();
  /** print MPI data structure */
  void Print(void);
  /** MPI status during the communication */
  MPI_Status status;
  /** rank of the process */
  int rank;    /*! in MPI_COMM_WORLD */
  /** number of processes */
  int nprocs;  /*! in MPI_COMM_WORLD */ 

  char *buffer;
  int buffer_size;
};
inline MPIdata::MPIdata(int *argc, char ***argv) {
  /* Initialize the MPI API */
  MPI_Init(argc, argv);

  /* Set rank, in MPI_COMM_WORLD */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Set total number of processors, in MPI_COMM_WORLD */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /*! mlmd: this is used in PSKOutput: decide wether MPI_COMM_WORLD or MPI_COMM_GRID
    is better here; at the moment, MPI_COMM_WORLD */
}

inline MPIdata::~MPIdata() {
}

inline void MPIdata::finalize_mpi() {
  MPI_Finalize();
}

inline void MPIdata::Print(void) {
  cout << endl;
  //cout << "Number of processes = " << nprocs << endl;
  cout << "-------------------------" << endl;
  cout << endl;
}

// extern MPIdata *mpi; // instantiated in iPIC3D.cpp

#endif
