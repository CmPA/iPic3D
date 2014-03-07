
#include "mpi.h"
#include "ComNodes3D.h"
#include "TimeTasks.h"
#include "ipicdefs.h"
#include "Alloc.h"
#include "debug.h"
#include "parallel.h"

/** communicate ghost cells (FOR NODES) */
void communicateNode(int nx, int ny, int nz, arr3_double _vector, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);


  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
/** communicate ghost cells (FOR NODES) */
void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector = _vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);

  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
/** communicate ghost cells (FOR NODES) with particles BC*/
void communicateNodeBC_P(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);


  // prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  // prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  // prepare ghost cell Edge X for communication: these are communicated in Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);

  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  // apply boundary condition to 8 Ghost Corners and communicate if necessary to 8 processors
  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);

  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}

/** SPECIES: communicate ghost cells */
void communicateNode(int nx, int ny, int nz, arr4_double _vector, int ns, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ****vector = _vector.fetch_arr4();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);



/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);


  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}

// PARTICLES
/** SPECIES: communicate ghost cells */
void communicateNode_P(int nx, int ny, int nz, arr4_double _vector, int ns, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ****vector = _vector.fetch_arr4();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);



/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);


  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}

// 
/** communicate ghost cells (FOR CENTERS) */
void communicateCenter(int nx, int ny, int nz, arr3_double _vector, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector = _vector.fetch_arr3();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);


/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);


  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
void communicateCenterBoxStencilBC(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);


  // deallocate
  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
}
// particles
/** communicate ghost cells (FOR CENTERS) with BOX stencil*/
void communicateCenterBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);


  // deallocate
  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
}

// 


void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);
  // deallocate
  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
}

void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();
  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeNodeFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);
  // deallocate
  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
}



/** SPECIES: communicate ghost cells */
void communicateCenter(int nx, int ny, int nz, arr4_double _vector, int ns, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ****vector=_vector.fetch_arr4();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ns, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ns, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ns, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ns, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, ns, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);


  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
// /////////// communication + BC ////////////////////////////
void communicateCenterBC(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor(), vct->getZleft_neighbor(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);

  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
// /////////// communication + BC ////////////////////////////
void communicateCenterBC_P(int nx, int ny, int nz, arr3_double _vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct) {
  timeTasks_set_communicating();
//  static int counter=0; if(is_output_thread()) { counter++; dprint(counter); }
  double ***vector=_vector.fetch_arr3();

  // allocate 6 ghost cell Faces
  double *ghostXrightFace = new double[(ny - 2) * (nz - 2)];
  double *ghostXleftFace = new double[(ny - 2) * (nz - 2)];
  double *ghostYrightFace = new double[(nx - 2) * (nz - 2)];
  double *ghostYleftFace = new double[(nx - 2) * (nz - 2)];
  double *ghostZrightFace = new double[(nx - 2) * (ny - 2)];
  double *ghostZleftFace = new double[(nx - 2) * (ny - 2)];
  // allocate 12 ghost cell Edges
  // X EDGE
  double *ghostXsameYleftZleftEdge = new double[nx - 2];
  double *ghostXsameYrightZleftEdge = new double[nx - 2];
  double *ghostXsameYleftZrightEdge = new double[nx - 2];
  double *ghostXsameYrightZrightEdge = new double[nx - 2];
  // Y EDGE
  double *ghostXrightYsameZleftEdge = new double[ny - 2];
  double *ghostXleftYsameZleftEdge = new double[ny - 2];
  double *ghostXrightYsameZrightEdge = new double[ny - 2];
  double *ghostXleftYsameZrightEdge = new double[ny - 2];
  // Z EDGE
  double *ghostXrightYleftZsameEdge = new double[nz - 2];
  double *ghostXrightYrightZsameEdge = new double[nz - 2];
  double *ghostXleftYleftZsameEdge = new double[nz - 2];
  double *ghostXleftYrightZsameEdge = new double[nz - 2];
  // allocate 8 ghost cell corner
  double ghostXrightYrightZrightCorner, ghostXleftYrightZrightCorner, ghostXrightYleftZrightCorner, ghostXleftYleftZrightCorner;
  double ghostXrightYrightZleftCorner, ghostXleftYrightZleftCorner, ghostXrightYleftZleftCorner, ghostXleftYleftZleftCorner;

  // apply boundary condition to 6 Ghost Faces and communicate if necessary to 6 processors: along 3 DIRECTIONS
  makeCenterFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);
  communicateGhostFace((ny - 2) * (nz - 2), vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightFace, ghostXleftFace);
  communicateGhostFace((nx - 2) * (nz - 2), vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostYrightFace, ghostYleftFace);
  communicateGhostFace((nx - 2) * (ny - 2), vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostZrightFace, ghostZleftFace);
  parseFace(nx, ny, nz, vector, ghostXrightFace, ghostXleftFace, ghostYrightFace, ghostYleftFace, ghostZrightFace, ghostZleftFace);

/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
  makeNodeEdgeY(nx, ny, nz, ghostZleftFace, ghostZrightFace, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
  makeNodeEdgeZ(nx, ny, nz, ghostXleftFace, ghostXrightFace, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
  makeNodeEdgeX(nx, ny, nz, ghostYleftFace, ghostYrightFace, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);

  // communicate twice each direction
  // X-DIRECTION: Z -> X
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZleftEdge, ghostXleftYsameZleftEdge);
  communicateGhostFace(ny - 2, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYsameZrightEdge, ghostXleftYsameZrightEdge);
  // Y-DIRECTION: X -> Y
  former_MPI_Barrier(MPI_COMM_WORLD);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXleftYrightZsameEdge, ghostXleftYleftZsameEdge);
  communicateGhostFace(nz - 2, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXrightYrightZsameEdge, ghostXrightYleftZsameEdge);
  // Z-DIRECTION: Y -> Z
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYleftZrightEdge, ghostXsameYleftZleftEdge);
  communicateGhostFace(nx - 2, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), ghostXsameYrightZrightEdge, ghostXsameYrightZleftEdge);
  // parse
  former_MPI_Barrier(MPI_COMM_WORLD);
  parseEdgeZ(nx, ny, nz, vector, ghostXrightYrightZsameEdge, ghostXleftYleftZsameEdge, ghostXrightYleftZsameEdge, ghostXleftYrightZsameEdge);
  parseEdgeY(nx, ny, nz, vector, ghostXrightYsameZrightEdge, ghostXleftYsameZleftEdge, ghostXleftYsameZrightEdge, ghostXrightYsameZleftEdge);
  parseEdgeX(nx, ny, nz, vector, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge);



  makeNodeCorner(nx, ny, nz, ghostXsameYrightZrightEdge, ghostXsameYleftZleftEdge, ghostXsameYleftZrightEdge, ghostXsameYrightZleftEdge, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  // communicate only in the X-DIRECTION
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);
  communicateGhostFace(1, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner);
  // parse
  parseCorner(nx, ny, nz, vector, &ghostXrightYrightZrightCorner, &ghostXleftYrightZrightCorner, &ghostXrightYleftZrightCorner, &ghostXleftYleftZrightCorner, &ghostXrightYrightZleftCorner, &ghostXleftYrightZleftCorner, &ghostXrightYleftZleftCorner, &ghostXleftYleftZleftCorner);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXright, bcFaceXleft, bcFaceYright, bcFaceYleft, bcFaceZright, bcFaceZleft, vct);

  delete[]ghostXrightFace;
  delete[]ghostXleftFace;
  delete[]ghostYrightFace;
  delete[]ghostYleftFace;
  delete[]ghostZrightFace;
  delete[]ghostZleftFace;
  // X EDGE
  delete[]ghostXsameYleftZleftEdge;
  delete[]ghostXsameYrightZleftEdge;
  delete[]ghostXsameYleftZrightEdge;
  delete[]ghostXsameYrightZrightEdge;
  // Y EDGE
  delete[]ghostXrightYsameZleftEdge;
  delete[]ghostXleftYsameZleftEdge;
  delete[]ghostXrightYsameZrightEdge;
  delete[]ghostXleftYsameZrightEdge;
  // Z EDGE
  delete[]ghostXrightYleftZsameEdge;
  delete[]ghostXrightYrightZsameEdge;
  delete[]ghostXleftYleftZsameEdge;
  delete[]ghostXleftYrightZsameEdge;
}
