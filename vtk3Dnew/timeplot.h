/*
 *  timeplot.h
 *
 *
 *  Created by Giovanni Lapenta on 7/29/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

int nxn, nyn, nzn;
int nxc, nyc, nzc;
int XLEN, YLEN, ZLEN;
// x position in processor topology. it will read the first cell in x in that processor
int jproc;
// z position in processor topology. it will read the first cell in x in that processor
int kproc;

int nproc;
int ns;
double* qom;

double Lx, Ly, Lz;

int MaxLevel;
int InitLevel;
int DeltaT;
int nlevels;

double *temp_storageX;
double *temp_storageY;
double *temp_storageZ;

// hdf stuff
hid_t    file_id;
hid_t    dataset_id;
herr_t   status;
