// Developed by Stefano Markidis, and Giovanni Lapenta
#ifndef _RESTART3D_H_
#define _RESTART3D_H_

#include "MPIdata.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid.h"
#include "Field.h"
#include "Particles3Dcomm.h"
#include "PSKhdf5adaptor.h"

#include <string>

using std::string;
using std::stringstream;


/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part);

/** this restart function writes the last restart with the last cycle */
void writeRESTART(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part, bool fool);

/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART_ES(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part);

/** this restart function writes the last restart with the last cycle */
void writeRESTART_ES(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part, bool fool);

#endif // _RESTART_H_
