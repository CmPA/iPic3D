
#ifndef __WOH5HUT_H__
#define __WOH5HUT_H__

#include <iomanip>
#include <string>

#ifdef USEH5HUT
#include "H5hut-io.h"
#endif
#include "Alloc.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "CollectiveIO.h"
#include "VCtopology3D.h"
#include "Particles3Dcomm.h"

void WriteOutputH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);

#endif
