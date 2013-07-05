
#ifndef __WOPARALLEL_H__
#define __WOPARALLEL_H__

#include <iomanip>
#include <string>

#include "phdf5.h"
#include "Alloc.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "CollectiveIO.h"
#include "VCtopology3D.h"

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle);

#endif
