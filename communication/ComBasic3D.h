/***************************************************************************
  ComBasic.h  -  Library to handle Basic Communication
  ------------------
 ***************************************************************************/

#ifndef ComBasic_H
#define ComBasic_H

#include <mpi.h>
#include "ComParser3D.h"


/** communicate particles along a direction **/
inline void communicateParticlesDIR(int buffer_size, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *b_right, double *b_left){
    MPI_Status status;
    double *LEN = new double[3];
    LEN[0] = XLEN, LEN[1] = YLEN, LEN[2] = ZLEN;
    switch (DIR){
        case 0:
            myrank = (int) floor((double)myrank/(ZLEN*YLEN));
            break;
        case 1:
            myrank = (int) floor((double)myrank/ZLEN);
            break;
    }
    if (myrank%2==0){
        // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
        if (right_neighbor != MPI_PROC_NULL){
            if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_right[0],buffer_size,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, MPI_COMM_WORLD, &status);
            else
                swapBuffer(buffer_size,b_left,b_right);
        }
    } else {
        // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
        if (left_neighbor != MPI_PROC_NULL){
            if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_left[0],buffer_size,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, MPI_COMM_WORLD, &status);
            else
                swapBuffer(buffer_size,b_left,b_right);
        }
    }
    if (myrank%2==1){
        // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
        if (right_neighbor != MPI_PROC_NULL){
            if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_right[0],buffer_size,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, MPI_COMM_WORLD, &status);

        }
    } else  {
        // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
        if (left_neighbor != MPI_PROC_NULL){
            if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_left[0],buffer_size,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, MPI_COMM_WORLD, &status);
        }
    }
    delete[] LEN;
}
/** communicate ghost along a direction **/
inline void communicateGhostFace(int b_len, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightFace, double *ghostLeftFace){

    MPI_Status status;
    double *LEN = new double[3];
    int rankF;
    LEN[0] = XLEN, LEN[1] = YLEN, LEN[2] = ZLEN;
    switch (DIR){
        case 0:
            rankF = (int) floor((double)myrank/(ZLEN*YLEN));
            break;
        case 1:
            rankF = (int) floor((double)myrank/ZLEN);
            break;
        case 2:
            rankF = myrank;
            break;
    }
    if (rankF%2==0 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)       //
        MPI_Sendrecv_replace(&ghostRightFace[0],b_len,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, MPI_COMM_WORLD, &status);
    else if (rankF%2== 1 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)           //
        MPI_Sendrecv_replace(&ghostLeftFace[0],b_len,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, MPI_COMM_WORLD, &status);

    if (rankF%2==1 && right_neighbor != MPI_PROC_NULL &&  LEN[DIR] > 1)             //
        MPI_Sendrecv_replace(&ghostRightFace[0],b_len,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, MPI_COMM_WORLD, &status);
    else if (rankF%2==0 && left_neighbor != MPI_PROC_NULL  && LEN[DIR] > 1 )   //
        MPI_Sendrecv_replace(&ghostLeftFace[0],b_len,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, MPI_COMM_WORLD, &status);
    // just swap the buffer if you have just a1 processor in 1 direction
    if (LEN[DIR] == 1 && right_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
        swapBuffer(b_len,ghostLeftFace,ghostRightFace);


    delete[] LEN;
}
/** communicate ghost edge along a direction; there are 6 Diagonal directions through which we exchange Ghost Edges :
  0 = from   XrightYrightZsame to YleftZleftZsame; we exchange Z edge
  1 = from   XrightYleftZsame to XleftYrightZsame; we exchange Z edge
  2 = from   XrightYsameZright to XleftYsameZsame; we exchange Y edge
  3 = from   XrightYsameZleft to XleftYsameZright; we exchange Y edge
  4 = from   XsameYrightZright to XsameYleftZleft; we exchange X edge
  5 = from   XsameYrightZleft to XsameYleftZright; we exchange X edge
  */
inline void communicateGhostEdge(int b_len, int myrank, int right_neighborD, int left_neighborD, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightEdge, double *ghostLeftEdge){
    MPI_Status status;
    int rankE;
    bool comNotDone = true;
    if (XLEN > 1 && YLEN > 1 && ZLEN > 1){
        if (DIR==0 || DIR==1) // Z same 
            rankE = (int) floor((double)myrank/(ZLEN*YLEN));
        else if (DIR==2 || DIR==3) // Y same
            rankE = (int) floor((double)myrank/(ZLEN*XLEN));
        else if (DIR==4 || DIR==5) // X same
            rankE = (int) floor((double)myrank/(ZLEN));
    }
    if (XLEN == 1 && YLEN > 1 && ZLEN > 1){
        if (DIR==0 || DIR==1) // Z same 
            rankE = (int) floor((double)myrank/ZLEN);
        else if (DIR==2 || DIR==3) // Y same
            rankE = myrank;
        else if (DIR==4 || DIR==5) // X same
            rankE = myrank; // not important for deadlock should just swap
    }
    if (XLEN > 1 && YLEN == 1 && ZLEN > 1){
        if (DIR==0 || DIR==1) // Z same 
            rankE = (int) floor((double)myrank/(ZLEN));
        else if (DIR==2 || DIR==3) // Y same
            rankE = myrank;
        else if (DIR==4 || DIR==5) // X same
            rankE = (int) floor((double)myrank/(ZLEN)); 
    }
    if (XLEN > 1 && YLEN > 1 && ZLEN == 1){
        if (DIR==0 || DIR==1) // Z same 
            rankE = myrank; // should not matter
        else if (DIR==2 || DIR==3) // Y same
            rankE = (int) floor((double)myrank/(YLEN));
        else if (DIR==4 || DIR==5) // X same
            rankE = (int) floor((double)myrank/(YLEN));
    } else {
        rankE = myrank;
    }


    // swap 
    if (right_neighborD == myrank && left_neighborD == myrank){
        swapBuffer(b_len,ghostLeftEdge,ghostRightEdge);
        comNotDone = false;
    }
    //If processors are available in diagonal comunicate in Diagonal: use rankE to avoid deadlocks
    if (rankE%2==0 && right_neighborD != MPI_PROC_NULL && comNotDone){    //
        MPI_Sendrecv_replace(&ghostRightEdge[0],b_len,MPI_DOUBLE,right_neighborD,1,right_neighborD,1, MPI_COMM_WORLD, &status);
    }
    else if (rankE%2== 1 && left_neighborD != MPI_PROC_NULL && comNotDone){           //
        MPI_Sendrecv_replace(&ghostLeftEdge[0],b_len,MPI_DOUBLE,left_neighborD,1,left_neighborD,1, MPI_COMM_WORLD, &status);
    }

    if (rankE%2==1 && right_neighborD != MPI_PROC_NULL && comNotDone){             //
        MPI_Sendrecv_replace(&ghostRightEdge[0],b_len,MPI_DOUBLE,right_neighborD,1,right_neighborD,1, MPI_COMM_WORLD, &status);

    }
    else if (rankE%2==0 && left_neighborD != MPI_PROC_NULL && comNotDone){
        MPI_Sendrecv_replace(&ghostLeftEdge[0],b_len,MPI_DOUBLE,left_neighborD,1,left_neighborD,1, MPI_COMM_WORLD, &status);

    }


}
/** Communicate ghost corners along a direction; there are 4 Diagonal directions through which we exchange Ghost Corners:
  0 =  from XrightYrightZright to XleftYleftZleft
  1 =  from XrightYleftZright  to XleftYrightZleft
  2 =  from XleftYrightZright  to XrightYleftZleft
  3 =  from XleftYleftZright   to XrightYrightZleft
  */
inline void communicateGhostCorner(int myrank, int right_neighborC, int left_neighborC, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightCorner, double *ghostLeftCorner){
    MPI_Status status;
    int rankC;
    bool comNotDone = true;
    if (XLEN > 1 && YLEN > 1 && ZLEN > 1){
        rankC = (int) floor((double)(myrank/(YLEN*ZLEN)) );
    } else if (XLEN == 1 && YLEN > 1 && ZLEN > 1){
        rankC = (int) floor((double)myrank/ZLEN);
    } else if (XLEN > 1 && YLEN == 1 && ZLEN > 1){
        rankC = (int) floor((double)myrank/(ZLEN));
    } else if (XLEN > 1 && YLEN > 1 && ZLEN == 1){
        rankC = (int) floor((double)myrank/(YLEN));
    } else {
        rankC = myrank;
    }


    // commmunicate
    if (right_neighborC == myrank && left_neighborC == myrank){
        swapBuffer(1,ghostLeftCorner,ghostRightCorner);
        comNotDone = false;
    }
    // if it's possible communicate corners
    if (rankC%2==0 && right_neighborC != MPI_PROC_NULL && comNotDone){     //
        MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, MPI_COMM_WORLD, &status);

    }
    else if (rankC%2== 1 && left_neighborC != MPI_PROC_NULL && comNotDone){           //
        MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, MPI_COMM_WORLD, &status);

    }

    if (rankC%2==1 && right_neighborC != MPI_PROC_NULL && comNotDone){             //
        MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, MPI_COMM_WORLD, &status);

    } else if (rankC%2==0 && left_neighborC != MPI_PROC_NULL && comNotDone){   //
        MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, MPI_COMM_WORLD, &status);

    }





}

















#endif
