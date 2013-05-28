#!/bin/bash
#MSUB -r iPic3D_sya # Request name
#MSUB -n 4 # Number of tasks to use
#MSUB -T 3600 # Elapsed time limit in seconds
#MSUB -o mpi.out # Standard output. %I is the job id
#MSUB -e mpi.err # Error output. %I is the job id
#MSUB -q xlarge # The queue choice

BRIDGE_MSUB_PWD=$CCCWORKDIR'/runs/gem1'
export BRIDGE_MSUB_PWD

set -x
cd ${BRIDGE_MSUB_PWD}
ccc_mprun ./iPIC3D inputfile
