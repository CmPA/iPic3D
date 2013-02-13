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

#include "mpi.h"
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
        MPIdata(int*, char***);
        /** destructor */
        ~MPIdata();
        /** initialize MPIdata */
        void init (int*, char***);
        /** close MPI environment */
        void finalize_mpi();
        /** print MPI data structure */
        void  Print(void);
        /** MPI status during the communication */
        MPI_Status status;
        /** rank of the process */
        int rank;
        /** number of processes */
        int nprocs;

        char* buffer;
        int buffer_size;
}  ;
inline MPIdata::MPIdata(int *argc, char*** argv){
    /* Initialize the MPI API
       In general, all nodes have a copy of all the variables defined.
       MPI_Init prepares the program run to communicate between all the
       nodes.  It is necessary to have this function call in all MPI
       code */
    MPI_Init(argc,argv);
    /* I ask how many other processors are out there = nprocs
       MPI_Comm_size initializes the numprocs variable to be the number
       of processors alloocated to the program run.  MPI_COMM_WORLD is
       a macro that MPI uses to address the correct MPI netowrk.  Since
       multiple MPI jobs could be running at once, we want to have a way
       of addressing only the processors in our program run (or our
       "world"). */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* Request my ID number = myrank
       MPI_Comm_rank initializes each nodes myrank variable to be it's
       processor number. */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

}
inline MPIdata::~MPIdata(){

}
inline void MPIdata::finalize_mpi(){
    MPI_Finalize();
}
inline void MPIdata::Print(void){
    cout << endl;
    cout << "Number of processes = " << nprocs << endl;
    cout << "-------------------------" << endl;
    cout << endl;
}

#endif

