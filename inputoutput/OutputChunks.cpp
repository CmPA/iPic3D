#include "OutputChunks.h"
#include <hdf5.h>


// PrintTest_2D
// #define DATASETNAME_2D "FloatData"
//#define RANK_2D 2   
//#define NX     7                      /* dataset dimensions */
//#define NY     3 


OutputChunks::OutputChunks(Collective * col, Grid * grid){

  //PrintTest_2D();
  
  // hard-coded, for now

  //PrintingChunks= true; // this should be read from inputfile

  PrintingChunks= false;  // this should be read from inputfile
  
  NSampling= 5;  

  NXSampling=4;
  NYSampling=3;
  NZSampling=1;
  // end hard-coded, for now


  if (PrintingChunks== false)
    return;

  X_index= new int[NSampling];
  Y_index= new int[NSampling];
  Z_index= new int[NSampling];
  
  X_Coord_first= new double[NSampling];
  Y_Coord_first= new double[NSampling];
  Z_Coord_first= new double[NSampling];

  
  // for now I am only sampling from 2D sims
  /* read grid info */
  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();

  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();

  ns= col->getNs();
      
  // where X_index, Y_index, Z_index are assigned
  AssignSamplingPoints();

  // check all the chunks to print are inside the core
  // also calculates the coordinate of the first point per direction
  CheckChunks_calcCoords(grid);

  /* now allocate field Chunks */
  /*Chunk_Ex= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);
  Chunk_Ey= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);
  Chunk_Ez= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);

  Chunk_Bxn= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);
  Chunk_Byn= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);
  Chunk_Bzn= newArr4(double, NSampling, NXSampling, NYSampling, NZSampling);

  Chunk_rhons= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);

  Chunk_Jxs= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_Jys= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_Jzs= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);

  Chunk_EFxs= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_EFys= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_EFzs= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);

  Chunk_pXXsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_pXYsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_pXZsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_pYYsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_pYYsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);
  Chunk_pYZsn= newArr5(double, NSampling, ns, NXSampling, NYSampling, NZSampling);*/
  /* end allocate */


}

// TO DO: read sampling points from somewhere; now, hard-coded
void OutputChunks::AssignSamplingPoints(){
  //here, the first indexes of the sampling points per core (X_index, Y_index, Z_index) are assigned

  if (PrintingChunks== false) return;
  
  X_index[0]= 3;
  Y_index[0]= 3;
  Z_index[0]= 2;

  X_index[1]= 3;
  Y_index[1]= nyn- NYSampling -3;
  Z_index[1]= 2;

  X_index[2]= nxn- NXSampling -3;
  Y_index[2]= nyn- NYSampling -3;
  Z_index[2]= 2;

  X_index[3]= nxn- NXSampling -3;;
  Y_index[3]= 3;
  Z_index[3]= 2;

  X_index[4]= floor(nxn/2);
  Y_index[4]= floor(nyn/2);
  Z_index[4]= 2;
  
  return;
  
}

// check that the chunks are inside the core, for simplicity
// if any of them is out, I set PrintingChunks= false
// which will prevent operations
// also calculates the inital coordinate of the points
void OutputChunks::CheckChunks_calcCoords(Grid * grid){

  if (PrintingChunks== false) return;

  for (int s=0; s<NSampling; s++){
    if (X_index[s] <1 or X_index[s]+NXSampling - 1 > nxn-2 )
      PrintingChunks= false;
    else
      X_Coord_first[s]= grid->getXN(X_index[s], 1, 1);

    if (Y_index[s]  <1 or Y_index[s] +NYSampling - 1 > nyn-2 )
      PrintingChunks= false;
    else
      Y_Coord_first[s]= grid->getYN(1, Y_index[s], 1);

    if (Z_index[s]  <1 or Z_index[s] +NZSampling - 1 > nzn-2 )
      PrintingChunks= false;
    else
      Z_Coord_first[s]= grid->getZN(1, 1, Z_index[s]);
  }
  
  return;
}
/*
// extract chunks to print
void OutputChunks::extractFieldChunk(double *** field, double **** chunk) {
  // input: whatever field
  // output

  for (int s=0; s< NSampling; s++){
    // improve this copy

    for (int i= 0; i< NXSampling; i ++)
      for (int j= 0; j< NYSampling; j ++)
	for (int k= 0; k< NZSampling; k ++)
	  chunk[s][i][j][k]= field[X_index[s]+i][Y_index[s]+j][Z_index[s]+k];
  }   
}

void OutputChunks::extractMomentChunk(double **** field, double ***** chunk) {
  // input: whatever field
  // outputs

  for (int s=0; s< NSampling; s++){
    // improve this copy
    for (int is=0; is< ns; is ++)
      for (int i= 0; i< NXSampling; i ++)
	for (int j= 0; j< NYSampling; j ++)
	  for (int k= 0; k< NZSampling; k ++)
	    chunk[s][is][i][j][k]= field[is][X_index[s]+i][Y_index[s]+j][Z_index[s]+k];
  }   
  }*/

// here, using EMf, I pack all the field
/*void OutputChunks::PackChunks(Field * EMf){

  if (PrintingChunks== false) return;

  cout << "inside OutputChunks::PackChunks" << endl;
  
  extractFieldChunk(EMf->Ex, Chunk_Ex); 
  extractFieldChunk(EMf->Ey, Chunk_Ey);
  extractFieldChunk(EMf->Ez, Chunk_Ez);

  extractFieldChunk(EMf->Bxn, Chunk_Bxn); 
  extractFieldChunk(EMf->Byn, Chunk_Byn);
  extractFieldChunk(EMf->Bzn, Chunk_Bzn);

  extractMomentChunk(EMf->rhons, Chunk_rhons);

  extractMomentChunk(EMf->Jxs, Chunk_Jxs);
  extractMomentChunk(EMf->Jys, Chunk_Jys);
  extractMomentChunk(EMf->Jzs, Chunk_Jzs);

  extractMomentChunk(EMf->EFxs, Chunk_EFxs);
  extractMomentChunk(EMf->EFys, Chunk_EFys);
  extractMomentChunk(EMf->EFzs, Chunk_EFzs);

  extractMomentChunk(EMf->pXXsn, Chunk_pXXsn);
  extractMomentChunk(EMf->pXYsn, Chunk_pXYsn);
  extractMomentChunk(EMf->pXZsn, Chunk_pXZsn);
  extractMomentChunk(EMf->pYYsn, Chunk_pYYsn);
  extractMomentChunk(EMf->pYZsn, Chunk_pYZsn);
  extractMomentChunk(EMf->pZZsn, Chunk_pZZsn);
  
  return;
  
  }*/

// the gets
int OutputChunks::get_NSampling(){
  return NSampling;
}

int OutputChunks::get_X_index(int s){
  return X_index[s];
}

int OutputChunks::get_Y_index(int s){
  return Y_index[s];
}

int OutputChunks::get_Z_index(int s){
  return Z_index[s];
}

int OutputChunks::get_NXSampling(){
  return NXSampling;
}
int OutputChunks::get_NYSampling(){
  return NYSampling;
}
int OutputChunks::get_NZSampling(){
  return NZSampling;
}

double OutputChunks::get_X_Coord_first(int s){
  return X_Coord_first[s];
}
double OutputChunks::get_Y_Coord_first(int s){
  return Y_Coord_first[s];
}
double OutputChunks::get_Z_Coord_first(int s){
  return Z_Coord_first[s];
}
double OutputChunks::get_dx(){
  return dx;
}
double OutputChunks::get_dy(){
  return dy;
}
double OutputChunks::get_dz(){
  return dz;
}
/*
double **** OutputChunks::get_ExChunk(){
  return Chunk_Ex;
}
double **** OutputChunks::get_EyChunk(){
  return Chunk_Ey;
}
double **** OutputChunks::get_EzChunk(){
  return Chunk_Ez;
}

double **** OutputChunks::get_BxnChunk(){
  return Chunk_Bxn;
}
double **** OutputChunks::get_BynChunk(){
  return Chunk_Byn;
}
double **** OutputChunks::get_BznChunk(){
  return Chunk_Bzn;
}
*/

// end the gets


