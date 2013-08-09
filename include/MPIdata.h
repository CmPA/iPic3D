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
 *
 * I made this class a singleton.  It should only be created once,
 * since MPI_Init should be called only once. -Alec
 */
class MPIdata {
public:
  static MPIdata& instance();
private:
  // disable constructor and destructor of this singleton
  // by making them private.
  ~MPIdata(){}
  MPIdata(){}
public:
  /** initialize MPI environment */
  static void init(int *, char ***);
  /** close MPI environment */
  void finalize_mpi();
  /** print MPI data structure */
  void Print(void);
  /** MPI status during the communication */
  MPI_Status status;
public:
  static int get_rank(){return instance().rank;}
  static int get_nprocs(){return instance().nprocs;}
private:
  /** rank of the process */
  static int rank;
  /** number of processes */
  static int nprocs;

  // evidently unused...
  char *buffer;
  int buffer_size;
};
#endif
