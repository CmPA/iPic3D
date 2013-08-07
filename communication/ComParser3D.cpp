
#include <mpi.h>
#include "ComParser3D.h"

/** swap the buffer */
void swapBuffer(int buffer_size, double *b_left, double *b_right) {
  double *temp = new double[buffer_size];
  for (register int i = 0; i < buffer_size; i++) {
    temp[i] = b_left[i];
    b_left[i] = b_right[i];
    b_right[i] = temp[i];
  }
  delete[]temp;
}
/** swap the buffer */
void swapBuffer(double *b_left, double *b_right) {
  double temp = *b_left;
  *b_left = *b_right;
  *b_right = temp;

}
/** swap ghost cells */
void swapGhostFace(int n1, int n2, double **ghostFaceLeft, double **ghostFaceRight) {
  double **temp = newArr2(double, n1, n2);
  for (register int i = 0; i < n1; i++) {
    for (register int j = 0; j < n2; j++) {
      temp[i][j] = ghostFaceLeft[i][j];
      ghostFaceLeft[i][j] = ghostFaceRight[i][j];
      ghostFaceRight[i][j] = temp[i][j];
    }
  }
  delArr2(temp, n1);
}
/** prepare ghost cells on 6 faces for communication of Nodes when ther is periodicity */
void makeNodeFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostXleftFace[counter] = vector[2][j][k];
      ghostXrightFace[counter] = vector[nx - 3][j][k];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostYleftFace[counter] = vector[i][2][k];
      ghostYrightFace[counter] = vector[i][ny - 3][k];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      ghostZleftFace[counter] = vector[i][j][2];
      ghostZrightFace[counter] = vector[i][j][nz - 3];
      counter++;
    }

}

/** prepare ghost cells on 6 faces for communication */
void makeNodeFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostXleftFace[counter] = vector[ns][2][j][k];
      ghostXrightFace[counter] = vector[ns][nx - 3][j][k];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostYleftFace[counter] = vector[ns][i][2][k];
      ghostYrightFace[counter] = vector[ns][i][ny - 3][k];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      ghostZleftFace[counter] = vector[ns][i][j][2];
      ghostZrightFace[counter] = vector[ns][i][j][nz - 3];
      counter++;
    }
}

/** prepare ghost cells on 6 faces for communication */
void makeCenterFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostXleftFace[counter] = vector[1][j][k];
      ghostXrightFace[counter] = vector[nx - 2][j][k];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostYleftFace[counter] = vector[i][1][k];
      ghostYrightFace[counter] = vector[i][ny - 2][k];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      ghostZleftFace[counter] = vector[i][j][1];
      ghostZrightFace[counter] = vector[i][j][nz - 2];
      counter++;
    }
}

// / SPECIES for interpolation
/** prepare ghost cells on 6 faces for communication */
void makeCenterFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostXleftFace[counter] = vector[ns][1][j][k];
      ghostXrightFace[counter] = vector[ns][nx - 2][j][k];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      ghostYleftFace[counter] = vector[ns][i][1][k];
      ghostYrightFace[counter] = vector[ns][i][ny - 2][k];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      ghostZleftFace[counter] = vector[ns][i][j][1];
      ghostZrightFace[counter] = vector[ns][i][j][nz - 2];
      counter++;
    }
}

// ////////////////////
// ///////////////////
// EDGES
// ///////////////////
// ///////////////////

/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
void makeNodeEdgeZ(int nx, int ny, int nz, double *faceXleft, double *faceXright, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge) {
  int counter = 0;
  int counterLeft = 0;
  int counterRight = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      if (j == 1) {             // YLEFT
        ghostXleftYleftZsameEdge[counterLeft] = faceXleft[counter];
        ghostXrightYleftZsameEdge[counterLeft] = faceXright[counter];
        counterLeft++;
      }
      if (j == ny - 2) {        // YRIGHT 
        ghostXleftYrightZsameEdge[counterRight] = faceXleft[counter];
        ghostXrightYrightZsameEdge[counterRight] = faceXright[counter];
        counterRight++;
      }
      counter++;
    }
}



// /
// 
// Y EDGE
/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
void makeNodeEdgeY(int nx, int ny, int nz, double *faceZleft, double *faceZright, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge) {
  int counter = 0;
  int counterLeft = 0;
  int counterRight = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      if (i == 1) {             // XLEFT
        ghostXleftYsameZleftEdge[counterLeft] = faceZleft[counter];
        ghostXleftYsameZrightEdge[counterLeft] = faceZright[counter];
        counterLeft++;
      }
      if (i == nx - 2) {        // XRIGHT 
        ghostXrightYsameZleftEdge[counterRight] = faceZleft[counter];
        ghostXrightYsameZrightEdge[counterRight] = faceZright[counter];
        counterRight++;
      }
      counter++;
    }
}



// 
// 
// X EDGE
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
void makeNodeEdgeX(int nx, int ny, int nz, double *faceYleft, double *faceYright, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge) {
  int counter = 0;
  int counterLeft = 0;
  int counterRight = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      if (k == 1) {             // ZLEFT
        ghostXsameYleftZleftEdge[counterLeft] = faceYleft[counter];
        ghostXsameYrightZleftEdge[counterLeft] = faceYright[counter];
        counterLeft++;
      }
      if (k == nz - 2) {        // ZRIGHT 
        ghostXsameYleftZrightEdge[counterRight] = faceYleft[counter];
        ghostXsameYrightZrightEdge[counterRight] = faceYright[counter];
        counterRight++;
      }
      counter++;
    }
}


// /////////////////////////////
// ////////////////////////////
// CORNER
// ///////////////////////////
// ///////////////////////////
/** prepare ghost cell Edge X for communication */
void makeNodeCorner(int nx, int ny, int nz, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner) {
  *ghostXleftYrightZrightCorner = ghostXsameYrightZrightEdge[0];
  *ghostXrightYrightZrightCorner = ghostXsameYrightZrightEdge[nx - 3];

  *ghostXleftYleftZrightCorner = ghostXsameYleftZrightEdge[0];
  *ghostXrightYleftZrightCorner = ghostXsameYleftZrightEdge[nx - 3];

  *ghostXleftYrightZleftCorner = ghostXsameYrightZleftEdge[0];
  *ghostXrightYrightZleftCorner = ghostXsameYrightZleftEdge[nx - 3];

  *ghostXleftYleftZleftCorner = ghostXsameYleftZleftEdge[0];
  *ghostXrightYleftZleftCorner = ghostXsameYleftZleftEdge[nx - 3];
}
// /////////////////////////////
// ////////////////////////////
// PARSE
// ////////////////////////////
// ////////////////////////////
/** insert the ghost cells in the 3D phisical vector */
void parseFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      vector[0][j][k] = ghostXleftFace[counter];
      vector[nx - 1][j][k] = ghostXrightFace[counter];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      vector[i][0][k] = ghostYleftFace[counter];
      vector[i][ny - 1][k] = ghostYrightFace[counter];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      vector[i][j][0] = ghostZleftFace[counter];
      vector[i][j][nz - 1] = ghostZrightFace[counter];
      counter++;
    }
}

/** insert the ghost cells in the 3D phisical vector */
void parseFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace) {
  // XFACES 
  int counter = 0;
  for (register int j = 1; j < ny - 1; j++)
    for (register int k = 1; k < nz - 1; k++) {
      vector[ns][0][j][k] = ghostXleftFace[counter];
      vector[ns][nx - 1][j][k] = ghostXrightFace[counter];
      counter++;
    }
  // YFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int k = 1; k < nz - 1; k++) {
      vector[ns][i][0][k] = ghostYleftFace[counter];
      vector[ns][i][ny - 1][k] = ghostYrightFace[counter];
      counter++;
    }
  // ZFACES
  counter = 0;
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++) {
      vector[ns][i][j][0] = ghostZleftFace[counter];
      vector[ns][i][j][nz - 1] = ghostZrightFace[counter];
      counter++;
    }

}



/** add the values of ghost cells faces to the 3D phisical vector */
void addFace(int nx, int ny, int nz, double ***vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace, VirtualTopology3D * vct) {
  int counter;


  if (vct->getXright_neighbor_P() != MPI_PROC_NULL) { // XRIGHT
    counter = 0;
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[nx - 2][j][k] += ghostXrightFace[counter];
        counter++;
      }
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL) {  // XLEFT
    counter = 0;
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[1][j][k] += ghostXleftFace[counter];
        counter++;
      }
  }

  if (vct->getYright_neighbor_P() != MPI_PROC_NULL) { // YRIGHT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[i][ny - 2][k] += ghostYrightFace[counter];
        counter++;
      }
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL) {  // Y LEFT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[i][1][k] += ghostYleftFace[counter];
        counter++;
      }
  }
  if (vct->getZright_neighbor_P() != MPI_PROC_NULL) { // ZRIGHT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int j = 1; j < ny - 1; j++) {
        vector[i][j][nz - 2] += ghostZrightFace[counter];
        counter++;
      }
  }
  if (vct->getZleft_neighbor_P() != MPI_PROC_NULL) {  // ZLEFT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int j = 1; j < ny - 1; j++) {
        vector[i][j][1] += ghostZleftFace[counter];
        counter++;
      }

  }

}
/** add the values of ghost cells faces to the 3D phisical vector */
void addFace(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace, double *ghostZrightFace, double *ghostZleftFace, VirtualTopology3D * vct) {
  int counter;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL) { // XRIGHT
    counter = 0;
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[ns][nx - 2][j][k] += ghostXrightFace[counter];
        counter++;
      }
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL) {  // XLEFT
    counter = 0;
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[ns][1][j][k] += ghostXleftFace[counter];
        counter++;
      }
  }

  if (vct->getYright_neighbor_P() != MPI_PROC_NULL) { // YRIGHT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[ns][i][ny - 2][k] += ghostYrightFace[counter];
        counter++;
      }
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL) {  // Y LEFT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int k = 1; k < nz - 1; k++) {
        vector[ns][i][1][k] += ghostYleftFace[counter];
        counter++;
      }
  }
  if (vct->getZright_neighbor_P() != MPI_PROC_NULL) { // ZRIGHT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int j = 1; j < ny - 1; j++) {
        vector[ns][i][j][nz - 2] += ghostZrightFace[counter];
        counter++;
      }
  }
  if (vct->getZleft_neighbor_P() != MPI_PROC_NULL) {  // ZLEFT
    counter = 0;
    for (register int i = 1; i < nx - 1; i++)
      for (register int j = 1; j < ny - 1; j++) {
        vector[ns][i][j][1] += ghostZleftFace[counter];
        counter++;
      }

  }
}
/** insert the ghost cells Edge Z in the 3D phisical vector */
void parseEdgeZ(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge) {
  for (register int i = 1; i < (nz - 1); i++) {
    vector[nx - 1][ny - 1][i] = ghostXrightYrightZsameEdge[i - 1];
    vector[0][0][i] = ghostXleftYleftZsameEdge[i - 1];
    vector[nx - 1][0][i] = ghostXrightYleftZsameEdge[i - 1];
    vector[0][ny - 1][i] = ghostXleftYrightZsameEdge[i - 1];

  }
}
/** insert the ghost cells Edge Z in the 3D phisical vector */
void parseEdgeZ(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge) {
  for (register int i = 1; i < (nz - 1); i++) {
    vector[ns][nx - 1][ny - 1][i] = ghostXrightYrightZsameEdge[i - 1];
    vector[ns][0][0][i] = ghostXleftYleftZsameEdge[i - 1];
    vector[ns][nx - 1][0][i] = ghostXrightYleftZsameEdge[i - 1];
    vector[ns][0][ny - 1][i] = ghostXleftYrightZsameEdge[i - 1];

  }
}


/** insert the ghost cells Edge Z in the 3D phisical vector */
void addEdgeZ(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge, VirtualTopology3D * vct) {

  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[ns][nx - 2][ny - 2][i] += ghostXrightYrightZsameEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[ns][1][1][i] += ghostXleftYleftZsameEdge[i - 1];
  }
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[ns][nx - 2][1][i] += ghostXrightYleftZsameEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[ns][1][ny - 2][i] += ghostXleftYrightZsameEdge[i - 1];
  }

}
/** insert the ghost cells Edge Z in the 3D phisical vector */
void addEdgeZ(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZsameEdge, double *ghostXleftYleftZsameEdge, double *ghostXrightYleftZsameEdge, double *ghostXleftYrightZsameEdge, VirtualTopology3D * vct) {
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[nx - 2][ny - 2][i] += ghostXrightYrightZsameEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[1][1][i] += ghostXleftYleftZsameEdge[i - 1];
  }
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[nx - 2][1][i] += ghostXrightYleftZsameEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nz - 1); i++)
      vector[1][ny - 2][i] += ghostXleftYrightZsameEdge[i - 1];
  }
}
/** prepare ghost cell Edge Y for communication */
void parseEdgeY(int nx, int ny, int nz, double ***vector, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge) {
  for (register int i = 1; i < ny - 1; i++) {
    vector[nx - 1][i][nz - 1] = ghostXrightYsameZrightEdge[i - 1];
    vector[0][i][0] = ghostXleftYsameZleftEdge[i - 1];
    vector[0][i][nz - 1] = ghostXleftYsameZrightEdge[i - 1];
    vector[nx - 1][i][0] = ghostXrightYsameZleftEdge[i - 1];
  }
}
/** prepare ghost cell Edge Y for communication */
void parseEdgeY(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge) {
  for (register int i = 1; i < (ny - 1); i++) {
    vector[ns][nx - 1][i][nz - 1] = ghostXrightYsameZrightEdge[i - 1];
    vector[ns][0][i][0] = ghostXleftYsameZleftEdge[i - 1];
    vector[ns][0][i][nz - 1] = ghostXleftYsameZrightEdge[i - 1];
    vector[ns][nx - 1][i][0] = ghostXrightYsameZleftEdge[i - 1];
  }
}
/** add the ghost cell values Edge Y to the 3D phisical vector */
void addEdgeY(int nx, int ny, int nz, double ***vector, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge, VirtualTopology3D * vct) {
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][nz - 2] += ghostXrightYsameZrightEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[1][i][1] += ghostXleftYsameZleftEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[1][i][nz - 2] += ghostXleftYsameZrightEdge[i - 1];
  }
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][1] += ghostXrightYsameZleftEdge[i - 1];
  }
}
/** SPECIES: add the ghost cell values Edge Y to the 3D phisical vector */
void addEdgeY(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYsameZrightEdge, double *ghostXleftYsameZleftEdge, double *ghostXleftYsameZrightEdge, double *ghostXrightYsameZleftEdge, VirtualTopology3D * vct) {
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[ns][nx - 2][i][nz - 2] += ghostXrightYsameZrightEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[ns][1][i][1] += ghostXleftYsameZleftEdge[i - 1];
  }
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[ns][1][i][nz - 2] += ghostXleftYsameZrightEdge[i - 1];
  }
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (ny - 1); i++)
      vector[ns][nx - 2][i][1] += ghostXrightYsameZleftEdge[i - 1];
  }

}
/** insert the ghost cells Edge X in the 3D physical vector */
void parseEdgeX(int nx, int ny, int nz, double ***vector, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge) {
  for (register int i = 1; i < (nx - 1); i++) {
    vector[i][ny - 1][nz - 1] = ghostXsameYrightZrightEdge[i - 1];
    vector[i][0][0] = ghostXsameYleftZleftEdge[i - 1];
    vector[i][0][nz - 1] = ghostXsameYleftZrightEdge[i - 1];
    vector[i][ny - 1][0] = ghostXsameYrightZleftEdge[i - 1];

  }
}
/** insert the ghost cells Edge X in the 3D phisical vector */
void parseEdgeX(int nx, int ny, int nz, double ****vector, int ns, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge) {
  for (register int i = 1; i < (nx - 1); i++) {
    vector[ns][i][ny - 1][nz - 1] = ghostXsameYrightZrightEdge[i - 1];
    vector[ns][i][0][0] = ghostXsameYleftZleftEdge[i - 1];
    vector[ns][i][0][nz - 1] = ghostXsameYleftZrightEdge[i - 1];
    vector[ns][i][ny - 1][0] = ghostXsameYrightZleftEdge[i - 1];
  }
}
/** add the ghost values Edge X to the 3D phisical vector */
void addEdgeX(int nx, int ny, int nz, double ***vector, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, VirtualTopology3D * vct) {
  if (vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][nz - 2] += ghostXsameYrightZrightEdge[i - 1];
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[i][1][1] += ghostXsameYleftZleftEdge[i - 1];
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[i][1][nz - 2] += ghostXsameYleftZrightEdge[i - 1];
  }
  if (vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][1] += ghostXsameYrightZleftEdge[i - 1];
  }
}
/** add the ghost values Edge X to the 3D phisical vector */
void addEdgeX(int nx, int ny, int nz, double ****vector, int ns, double *ghostXsameYrightZrightEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrightEdge, double *ghostXsameYrightZleftEdge, VirtualTopology3D * vct) {
  if (vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[ns][i][ny - 2][nz - 2] += ghostXsameYrightZrightEdge[i - 1];
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[ns][i][1][1] += ghostXsameYleftZleftEdge[i - 1];
  }
  if (vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[ns][i][1][nz - 2] += ghostXsameYleftZrightEdge[i - 1];
  }
  if (vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL) {
    for (register int i = 1; i < (nx - 1); i++)
      vector[ns][i][ny - 2][1] += ghostXsameYrightZleftEdge[i - 1];
  }
}
// /////////////////
// ////////////////
// PARSE CORNERS
// ///////////////
// //////////////

/** insert the ghost cells Edge X in the 3D phisical vector */
void parseCorner(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner) {
  vector[nx - 1][ny - 1][nz - 1] = *ghostXrightYrightZrightCorner;
  vector[0][ny - 1][nz - 1] = *ghostXleftYrightZrightCorner;
  vector[nx - 1][0][nz - 1] = *ghostXrightYleftZrightCorner;
  vector[0][0][nz - 1] = *ghostXleftYleftZrightCorner;
  vector[nx - 1][ny - 1][0] = *ghostXrightYrightZleftCorner;
  vector[0][ny - 1][0] = *ghostXleftYrightZleftCorner;
  vector[nx - 1][0][0] = *ghostXrightYleftZleftCorner;
  vector[0][0][0] = *ghostXleftYleftZleftCorner;
}
/** insert the ghost cells Edge X in the 3D physical vector */
void parseCorner(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner) {
  vector[ns][nx - 1][ny - 1][nz - 1] = *ghostXrightYrightZrightCorner;
  vector[ns][0][ny - 1][nz - 1] = *ghostXleftYrightZrightCorner;
  vector[ns][nx - 1][0][nz - 1] = *ghostXrightYleftZrightCorner;
  vector[ns][0][0][nz - 1] = *ghostXleftYleftZrightCorner;
  vector[ns][nx - 1][ny - 1][0] = *ghostXrightYrightZleftCorner;
  vector[ns][0][ny - 1][0] = *ghostXleftYrightZleftCorner;
  vector[ns][nx - 1][0][0] = *ghostXrightYleftZleftCorner;
  vector[ns][0][0][0] = *ghostXleftYleftZleftCorner;
}
/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ***vector, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner, VirtualTopology3D * vct) {
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[nx - 2][ny - 2][nz - 2] += *ghostXrightYrightZrightCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[1][ny - 2][nz - 2] += *ghostXleftYrightZrightCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[nx - 2][1][nz - 2] += *ghostXrightYleftZrightCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[1][1][nz - 2] += *ghostXleftYleftZrightCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[nx - 2][ny - 2][1] += *ghostXrightYrightZleftCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[1][ny - 2][1] += *ghostXleftYrightZleftCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[nx - 2][1][1] += *ghostXrightYleftZleftCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[1][1][1] += *ghostXleftYleftZleftCorner;

}
/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ****vector, int ns, double *ghostXrightYrightZrightCorner, double *ghostXleftYrightZrightCorner, double *ghostXrightYleftZrightCorner, double *ghostXleftYleftZrightCorner, double *ghostXrightYrightZleftCorner, double *ghostXleftYrightZleftCorner, double *ghostXrightYleftZleftCorner, double *ghostXleftYleftZleftCorner, VirtualTopology3D * vct) {
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[ns][nx - 2][ny - 2][nz - 2] += *ghostXrightYrightZrightCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[ns][1][ny - 2][nz - 2] += *ghostXleftYrightZrightCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[ns][nx - 2][1][nz - 2] += *ghostXrightYleftZrightCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZright_neighbor_P() != MPI_PROC_NULL)
    vector[ns][1][1][nz - 2] += *ghostXleftYleftZrightCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[ns][nx - 2][ny - 2][1] += *ghostXrightYrightZleftCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYright_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[ns][1][ny - 2][1] += *ghostXleftYrightZleftCorner;
  if (vct->getXright_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[ns][nx - 2][1][1] += *ghostXrightYleftZleftCorner;
  if (vct->getXleft_neighbor_P() != MPI_PROC_NULL && vct->getYleft_neighbor_P() != MPI_PROC_NULL && vct->getZleft_neighbor_P() != MPI_PROC_NULL)
    vector[ns][1][1][1] += *ghostXleftYleftZleftCorner;
}
