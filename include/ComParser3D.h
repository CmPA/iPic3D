/***************************************************************************
  ComParser3D.h  -  
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef ComParser_H
#define ComParser_H

#include <math.h>

#include "Alloc.h"
#include "VirtualTopology3D.h"

/** swap the buffer */
void swapBuffer(int buffer_size, double *b_left, double *b_right);

/** swap the buffer */
void swapBuffer(double *b_left, double *b_right);

/** swap ghost cells */
void swapGhostFace(int n1, int n2, double **ghostFaceLeft, double **ghostFaceRight);

/** prepare ghost cells on 6 faces for communication of Nodes when ther is periodicity */
void makeNodeFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

/** prepare ghost cells on 6 faces for communication */
void makeNodeFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

/** prepare ghost cells on 6 faces for communication */
void makeCenterFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

// / SPECIES for interpolation
/** prepare ghost cells on 6 faces for communication */
void makeCenterFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

// ////////////////////
// ///////////////////
// EDGES
// ///////////////////
// ///////////////////

/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
void makeNodeEdgeZ(int nx, int ny, int nz, double *faceXleft, double *faceXright, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge);

// /
// 
// Y EDGE
/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
void makeNodeEdgeY(int nx, int ny, int nz, double *faceZleft, double *faceZright, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge);

// 
// 
// X EDGE
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
void makeNodeEdgeX(int nx, int ny, int nz, double *faceYleft, double *faceYright, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge);

// /////////////////////////////
// ////////////////////////////
// CORNER
// ///////////////////////////
// ///////////////////////////
/** prepare ghost cell Edge X for communication */
void makeNodeCorner(int nx, int ny, int nz, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner);

// /////////////////////////////
// ////////////////////////////
// PARSE
// ////////////////////////////
// ////////////////////////////
/** insert the ghost cells in the 3D phisical vector */
void parseFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

/** insert the ghost cells in the 3D phisical vector */
void parseFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace);

/** add the values of ghost cells faces to the 3D phisical vector */
void addFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace, VirtualTopology3D * vct);

/** add the values of ghost cells faces to the 3D phisical vector */
void addFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace, VirtualTopology3D * vct);

/** insert the ghost cells Edge Z in the 3D phisical vector */
void parseEdgeZ(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge);

/** insert the ghost cells Edge Z in the 3D phisical vector */
void parseEdgeZ(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge);

/** insert the ghost cells Edge Z in the 3D phisical vector */
void addEdgeZ(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge, VirtualTopology3D * vct);

/** insert the ghost cells Edge Z in the 3D phisical vector */
void addEdgeZ(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge, VirtualTopology3D * vct);

/** prepare ghost cell Edge Y for communication */
void parseEdgeY(int nx, int ny, int nz, double ***vector, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge);

/** prepare ghost cell Edge Y for communication */
void parseEdgeY(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge);

/** add the ghost cell values Edge Y to the 3D phisical vector */
void addEdgeY(int nx, int ny, int nz, double ***vector, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge, VirtualTopology3D * vct);

/** SPECIES: add the ghost cell values Edge Y to the 3D phisical vector */
void addEdgeY(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge, VirtualTopology3D * vct);

/** insert the ghost cells Edge X in the 3D physical vector */
void parseEdgeX(int nx, int ny, int nz, double ***vector, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge);

/** insert the ghost cells Edge X in the 3D phisical vector */
void parseEdgeX(int nx, int ny, int nz, double ****vector, int ns, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge);

/** add the ghost values Edge X to the 3D phisical vector */
void addEdgeX(int nx, int ny, int nz, double ***vector, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, VirtualTopology3D * vct);

/** add the ghost values Edge X to the 3D phisical vector */
void addEdgeX(int nx, int ny, int nz, double ****vector, int ns, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, VirtualTopology3D * vct);

// /////////////////
// ////////////////
// PARSE CORNERS
// ///////////////
// //////////////

/** insert the ghost cells Edge X in the 3D phisical vector */
void parseCorner(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner);

/** insert the ghost cells Edge X in the 3D physical vector */
void parseCorner(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner);

/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner, VirtualTopology3D * vct);

/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner, VirtualTopology3D * vct);

#endif
