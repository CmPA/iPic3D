/*************************************************************************** 
 OutputChucks.h  -  output chunks of the simulation, at all cycles        
 -------------------                   
begin                : March 4, 2021   
copyright            : (C) ???
developer            : Maria Elena Innocenti
email                : ***
***************************************************************************/

#ifndef OutputChunks_H
#define OutputChunks_H

#include "mpi.h"
#include "Collective.h"
#include <hdf5.h>
#include "Grid.h"
#include "Field.h"


class OutputChunks                                            
{
  
public:
  /** constructor */
  OutputChunks(Collective * col, Grid * grid);
  /** destructor */
  ~OutputChunks();

  // extract all field chunks                                                                                                                                                   
  //void PackChunks(Field * EMf);

  // checked a number of time, proceeds only if True
  // if there are issues with the points, false
  bool PrintingChunks;

  // the gets
  int get_NSampling();
  
  int get_X_index(int);
  int get_Y_index(int);
  int get_Z_index(int);

  int get_NXSampling();
  int get_NYSampling();
  int get_NZSampling();

  double get_X_Coord_first(int);
  double get_Y_Coord_first(int);
  double get_Z_Coord_first(int);

  double get_dx();
  double get_dy();
  double get_dz();

  /*double ****get_ExChunk();
  double ****get_EyChunk();
  double ****get_EzChunk();

  double ****get_BxnChunk();
  double ****get_BynChunk();
  double ****get_BznChunk();

  double ***&get_rhonsChunk(int s, int is);

  double ***&get_JxsChunk(int s, int is);
  double ***&get_JysChunk(int s, int is);
  double ***&get_JzsChunk(int s, int is);

  double ***&get_EFxsChunk(int s, int is);
  double ***&get_EFysChunk(int s, int is);
  double ***&get_EFzsChunk(int s, int is);

  double ***&get_pXXsnChunk(int s, int is);
  double ***&get_pXYsnChunk(int s, int is);
  double ***&get_pXZsnChunk(int s, int is);
  double ***&get_pYYsnChunk(int s, int is);
  double ***&get_pYZsnChunk(int s, int is);
  double ***&get_pZZsnChunk(int s, int is);*/

  
private:
  
  /* how many sampling points per core */
  int NSampling;

  /* how many points to save in each direction */
  int NXSampling;
  int NYSampling;
  int NZSampling;

  /* index of the first point of the chunk, per direction  */
  int *X_index;
  int *Y_index;
  int *Z_index;

  /* coordinate of the first point of the chunk ** IN THE ENTIRE GRID **, 
     per direction (this is just to assist in the post-processing) */
  double *X_Coord_first;
  double *Y_Coord_first;
  double *Z_Coord_first;

  int nxn, nyn, nzn;

  double dx, dy, dz;

  // number of particle species
  int ns;
  /* where the fields are stored*/

  /*double ****Chunk_Ex;
  double ****Chunk_Ey;
  double ****Chunk_Ez;

  double ****Chunk_Bxn;
  double ****Chunk_Byn;
  double ****Chunk_Bzn;

  double *****Chunk_rhons;
  
  double *****Chunk_Jxs;
  double *****Chunk_Jys;
  double *****Chunk_Jzs;

  double *****Chunk_EFxs;
  double *****Chunk_EFys;
  double *****Chunk_EFzs;

  double *****Chunk_pXXsn;
  double *****Chunk_pXYsn;
  double *****Chunk_pXZsn;
  double *****Chunk_pYYsn;
  double *****Chunk_pYZsn;
  double *****Chunk_pZZsn;*/
    
  /* end where the fields are stored */


  // here, the first indexes of the sampling points per core (X_index, Y_index, Z_index) are assigned
  void AssignSamplingPoints();

  // check that the chunks are inside the core, for simplicity
  // also calculates the coordinate of the first point per direction                
  void CheckChunks_calcCoords(Grid * grid);

  /*
  // extract field chunk
  void extractFieldChunk(double *** field, double **** chunk) ;
  // extract moment chunk
  void extractMomentChunk(double **** field, double ***** chunk);*/
  
};
#endif
