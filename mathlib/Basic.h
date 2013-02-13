/*******************************************************************************************
  Basic.h  -  Basic operations 
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/
#ifndef Basic_H
#define Basic_H

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


/**
 *  
 * Basic operations defined. This library provides methods to calculate:
 *
 * - dot product of two vectors
 * - square norm of a vector
 * - norm of a vector
 * - difference of two vector
 * - scalar-vector product
 * - vector1 = vector1 + alfa*vector2
 * - vector1 = beta*vector1 + alfa*vector2
 * - opposite of a vector
 *
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */


/** method to calculate the parallel dot product with vect1, vect2 having the ghost cells*/
inline double dotP(double *vect1, double *vect2,int n){
    double result = 0;
    double local_result = 0;
    for (register int i=0; i< n; i++)
        local_result += vect1[i]*vect2[i];
    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return(result);

}
/** method to calculate dot product */
inline double dot(double *vect1, double *vect2, int n){
    double result = 0;
    for (int i=0; i<n; i++)
        result += vect1[i]*vect2[i];
    return(result);
}
/** method to calculate the square norm of a vector */
inline double norm2(double **vect, int nx, int ny){
    double result = 0;
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            result += vect[i][j]*vect[i][j];
    return(result);
}
/** method to calculate the square norm of a vector */
inline double norm2(double ***vect, int nx, int ny){
    double result = 0;
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            result += vect[i][j][0]*vect[i][j][0];
    return(result);
}
/** method to calculate the square norm of a vector */
inline double norm2(double *vect, int nx){
    double result = 0;
    for (int i=0; i<nx; i++)
        result += vect[i]*vect[i];
    return(result);
}



/** method to calculate the parallel dot product */
inline double norm2P(double ***vect, int nx, int ny, int nz){
    double result = 0;
    double local_result = 0;
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                local_result += vect[i][j][k]*vect[i][j][k];

    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return(result);
}
/** method to calculate the parallel norm of a vector on different processors with the ghost cell */
inline double norm2P(double *vect, int n){
    double result = 0;
    double local_result = 0;
    for (int i=0; i< n; i++)
        local_result += vect[i]*vect[i];
    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return(result);
}
/** method to calculate the parallel norm of a vector on different processors with the gost cell*/
inline double normP(double *vect, int n){
    double result = 0.0;
    double local_result = 0.0;
    for (register int i=0; i<n; i++)
        local_result += vect[i]*vect[i];


    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return(sqrt(result));

}
/** method to calculate the difference of two vectors*/
inline void sub(double *res, double *vect1, double *vect2, int n){
    for (register int i=0; i<n; i++)
        res[i] = vect1[i] - vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double *vect1, double *vect2,int n){
    for (register int i=0; i<n; i++)
        vect1[i] += vect2[i];


}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ***vect2, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] += vect2[i][j][k];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ***vect2, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0] += vect2[i][j][0];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ****vect2, int nx, int ny, int nz,int ns){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] += vect2[ns][i][j][k];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ****vect2, int nx, int ny,int ns){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0] += vect2[ns][i][j][0];
}
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(double ***vect1, double ***vect2, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] -= vect2[i][j][k];


}

/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(double ***vect1, double ***vect2, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0] -= vect2[i][j][0];


}


/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 + gamma*vector3 + delta*vector4 */
inline void sum4(double ***vect1, double alfa, double ***vect2, double beta, double ***vect3, double gamma, double ***vect4, double delta, double ***vect5, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] = alfa*(vect2[i][j][k] + beta*vect3[i][j][k] + gamma*vect4[i][j][k] + delta*vect5[i][j][k]);

}
/** method to calculate the scalar-vector product */
inline void scale(double* vect, double alfa, int n){
    for (register int i=0; i<n; i++)
        vect[i]*=alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double***  vect, double alfa, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect[i][j][0] *=alfa;
}


/** method to calculate the scalar-vector product */
inline void scale(double*** vect, double alfa, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect[i][j][k] *=alfa;
}
/** method to calculate the scalar product */
inline void scale(double vect[][2][2], double alfa, int nx, int ny, int nz){
    for (int i=0; i < nx; i++)
        for (int j=0; j < ny ; j++)
            for (int k =0; k < nz; k++)
                vect[i][j][k] *= alfa;
}
/** method to calculate the scalar-vector product */
inline void scale(double***  vect1, double***  vect2,double alfa, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] = vect2[i][j][k]*alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double***  vect1, double***  vect2,double alfa, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0] = vect2[i][j][0]*alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double *vect1, double *vect2,double alfa, int n){
    for (register int i=0; i<n; i++)
        vect1[i] = vect2[i]*alfa;
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double*** vect1, double*** vect2,int nx,int ny,int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] = vect1[i][j][k] +  alfa*vect2[i][j][k];
}
/** add scale for weights */
inline void addscale(double alfa, double vect1[][2][2], double vect2[][2][2],int nx, int ny, int nz){
    for (int i=0; i < nx;i++)
        for (int j=0;j < ny;j++)
            for (int k=0;k < nz; k++)
                vect1[i][j][k] = vect1[i][j][k] + alfa*vect2[i][j][k];

}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double*** vect1, double*** vect2,int nx,int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0]+= alfa*vect2[i][j][0];
}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double*  vect1, double*  vect2, int n){
    for (register int i=0; i < n; i++)
        vect1[i] += alfa*vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
inline void addscale(double alfa,double beta, double*  vect1, double*  vect2, int n){
    for (register int i=0; i < n; i++)
        vect1[i] = vect1[i]*beta + alfa*vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, double *** vect1, double *** vect2,int nx, int ny, int nz){

    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++){
                vect1[i][j][k] = beta*vect1[i][j][k]  + alfa*vect2[i][j][k];
            }

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, double ***vect1, double ***vect2,int nx, int ny){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            vect1[i][j][0] = beta*vect1[i][j][0] + alfa*vect2[i][j][0];

}


/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
inline void scaleandsum(double ***vect1, double alfa,double beta, double ***vect2, double ***vect3, int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = alfa*vect2[i][j][k] + beta*vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2 depending on species*/
inline void scaleandsum(double ***vect1, double alfa,double beta, double ****vect2, double ***vect3,int ns, int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = alfa*vect2[ns][i][j][k] + beta*vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on species*/
inline void prod(double ***vect1, double alfa,double ****vect2, int ns, double ***vect3, int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = alfa*vect2[ns][i][j][k]*vect3[i][j][k];

}
/** method to calculate vect1 = vect2/alfa */
inline void div(double ***vect1, double alfa, double ***vect2,int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = vect2[i][j][k]/alfa;

}
inline void prod6(double ***vect1,double ***vect2,double ***vect3,double ***vect4,double ***vect5,double ***vect6,double ***vect7,int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = vect2[i][j][k]*vect3[i][j][k] + vect4[i][j][k]*vect5[i][j][k] + vect6[i][j][k]*vect7[i][j][k];
}
/** method used for calculating PI */
inline void proddiv(double ***vect1, double ***vect2, double alfa, double ***vect3, double ***vect4, double ***vect5, double ***vect6, double beta, double ***vect7, double ***vect8, double gamma, double ***vect9,int nx, int ny, int nz){
    for (register int i=0; i< nx; i++)
        for (register int j=0; j< ny; j++)
            for (register int k=0; k< nz; k++)
                vect1[i][j][k] = (vect2[i][j][k] + alfa*(vect3[i][j][k]*vect4[i][j][k]  - vect5[i][j][k]*vect6[i][j][k]) + beta*vect7[i][j][k]*vect8[i][j][k])/(1+gamma*vect9[i][j][k]);

    // questo mi convince veramente poco!!!!!!!!!!!!!! CAZZO!!!!!!!!!!!!!!!!!!
    //  ***vect1++ = (***vect2++ + alfa*((***vect3++)*(***vect4++)  - (***vect5++)*(***vect6++)) + beta*(***vect7++)*(***vect8++))/(1+gamma*(***vect9++));
}  
/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx, int ny, int nz){
    for (register int i=0; i<nx;i++)
        for (register int j=0; j<ny;j++)
            for (register int k=0; k<nz;k++)
                vect[i][j][k] = -vect[i][j][k];


}

/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx, int ny){
    for (register int i=0; i<nx;i++)
        for (register int j=0; j<ny;j++)
            vect[i][j][0] = -vect[i][j][0];
}
/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx){
    for (register int i=0; i<nx;i++)
        vect[i][0][0] = -vect[i][0][0];
}
/** method to calculate the opposite of a vector */
inline void neg(double *vect, int n){
    for (register int i=0; i<n; i++)
        vect[i] = -vect[i];


}
/** method to set equal two vectors */
inline void eq(double ***vect1, double ***vect2, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[i][j][k] = vect2[i][j][k];

}
/** method to set equal two vectors */
inline void eq(double ***vect1, double ***vect2, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[i][j][0] = vect2[i][j][0];

}

/** method to set equal two vectors */
inline void eq(double ****vect1, double ***vect2, int nx, int ny, int is){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect1[is][i][j][0] = vect2[i][j][0];

}
/** method to set equal two vectors */
inline void eq(double ****vect1, double ***vect2, int nx, int ny, int nz, int is){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect1[is][i][j][k] = vect2[i][j][k];

}

inline void eq(double *vect1, double *vect2, int n){
    for (register int i=0; i < n; i++)
        vect1[i] = vect2[i];
}
/** method to set a vector to a Value */
inline void eqValue(double value, double ***vect, int nx, int ny, int nz){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            for (register int k=0; k<nz; k++)
                vect[i][j][k] = value;

}
inline void eqValue(double value, double vect[][2][2], int nx, int ny, int nz){
    for(int i=0; i < nx;i++)
        for(int j=0; j < ny; j++)
            for(int k=0;k < nz; k++)
                vect[i][j][k] = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double ***vect, int nx, int ny){
    for (register int i=0; i<nx; i++)
        for (register int j=0; j<ny; j++)
            vect[i][j][0]= value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double ***vect, int nx){
    for (register int i=0; i<nx; i++)
        vect[i][0][0]= value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double *vect, int n){
    for (register int i=0; i<n; i++)
        vect[i] = value;
}
/** method to put a column in a matrix 2D */
inline void putColumn(double **Matrix, double *vect, int column, int n){
    for (int i=0; i < n; i++)
        Matrix[i][column] = vect[i];

}
/** method to get a column in a matrix 2D */
inline void getColumn(double *vect, double **Matrix, int column, int n){
    for (int i=0; i < n; i++)
        vect[i] =  Matrix[i][column];
}
/** RIFAI QUESTA PARTE questo e' la tomba delle performance*/
inline void MODULO(double *x, double L){
    *x = *x - floor(*x/L)*L;

}
/** method to calculate the epsilon machine */
inline double eps(){
    double eps;
    int i=1;
    double num = 1;
    double newsum=1;
    double oldsum=1;
    while (true){
        num = num/(2*i);
        newsum+=num;
        if (newsum==oldsum)
            break;
        oldsum = newsum;
        i++ ;
    }
    eps =  num*2;
    return(eps);
}
/** method to get rid of the ghost cells */
inline void getRidGhost(double** out, double** in, int nx, int ny){
    for (register int i=1; i<nx-1; i++)
        for (register int j=1; j<ny-1; j++)
            out[i-1][j-1]=in[i][j];  
}    

#endif
