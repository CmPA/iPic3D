/***************************************************************************
  BcParticles.h  -  Library to manage boundary conditions for particles
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef BcParticles_H
#define BcParticles_H
/** set the boundary VirtualTopology3D3Dcondition  for particle in 2D
  <ul>
  <li>bc = 1 perfect mirror </li>
  <li>bc = 2 riemission     </li>
  </ul>*/
inline void BCpart(double *x,double *u,double Lx,double ut,int bcFaceXright, int bcFaceXleft){
    if (*x > Lx){
        switch(bcFaceXright){
            case 1:   // perfect mirror
                *x = 2*Lx -*x;
                *u = - *u;
                break;
            case 2:   // riemmission
                double harvest, prob,theta;
                *x = 2*Lx -*x;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = - fabs(ut*prob*cos(theta));



        }
    } else if (*x < 0) {
        switch(bcFaceXleft){
            case 1:   // perfect mirror
                MODULO(x,Lx);
                *u = -*u;
                break;
            case 2:   // riemmission
                MODULO(x,Lx);
                double harvest, prob,theta;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = fabs(ut*prob*cos(theta));
                break;

        }

    }

}
/** set the boundary condition  for particle in 2D
  <ul>
  <li>bc = 1 perfect mirror </li>
  <li>bc = 2 riemission     </li>
  </ul>*/
inline void BCpart(double *x, double *u, double *v, double *w, double Lx,double ut, double vt, double wt,  int bcFaceXright, int bcFaceXleft){
    if (*x > Lx){
        switch(bcFaceXright){
            case 1:   // perfect mirror
                *x = 2*Lx - *x;
                *u = - *u;
                break;
            case 2:   // riemmission
                double harvest, prob,theta;
                *x = 2*Lx - *x;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = - fabs(ut*prob*cos(theta));
                // v
                *v = vt*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta); 


        }
    } else if (*x < 0) {
        switch(bcFaceXleft){
            case 1:   // perfect mirror
                *x = -*x;
                *u = -*u;
                break;
            case 2:   // riemmission
                *x = -*x;
                double harvest, prob,theta;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = fabs(ut*prob*cos(theta));
                // v
                *v = vt*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta);   
                break;

        }

    }

}
/** set the boundary condition on boundaries for particle in 3D*/
inline void BCpart(double *x, double *y, double *z, double *u, double *v, double *w, double Lx, double Ly, double Lz, double ut, double vt, double wt,  int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft,int bcFaceZright,int bcFaceZleft, VirtualTopology3D *vct){
    if (*x > Lx && vct->getXright_neighbor()==MPI_PROC_NULL){
        switch(bcFaceXright){
            case 1:   // perfect mirror
                MODULO(x,Lx);
                *u = - *u;
                break;
            case 2:   // riemmission
                double harvest, prob,theta;
                MODULO(x,Lx);
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = - fabs(ut*prob*cos(theta));
                // v
                *v = vt*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta);   
                break;



        }
    }
    if (*x < 0 && vct->getXleft_neighbor()==MPI_PROC_NULL){
        switch(bcFaceXleft){
            case 1:   // perfect mirror
                MODULO(x,Lx);
                *u = -*u;
                break;
            case 2:   // riemmission
                MODULO(x,Lx);
                double harvest, prob,theta;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = fabs(ut*prob*cos(theta));
                // v
                *v = vt*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta);   
                break;



        }
    }
    if (*y > Ly && vct->getYright_neighbor()==MPI_PROC_NULL){
        switch(bcFaceYright){
            case 1:   // perfect mirror
                MODULO(y,Ly);
                *v = -*v;
                break;
            case 2:   // riemmission
                MODULO(y,Ly);
                double harvest, prob,theta;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = ut*prob*cos(theta);
                // v
                *v = -fabs(vt*prob*sin(theta));
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta); 
                break;



        }
    }
    if (*y < 0 && vct->getYleft_neighbor()==MPI_PROC_NULL){
        switch(bcFaceYleft){
            case 1:   // perfect mirror
                MODULO(y,Ly);
                *v = -*v;
                break;
            case 2:   // riemmission
                MODULO(y,Ly);
                double harvest, prob,theta;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *u = ut*prob*cos(theta);
                // v
                *v = fabs(vt*prob*sin(theta));
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                *w = wt*prob*cos(theta);
                break;



        }
    }
    /*        if (*z > Lz && vct->getZright_neighbor()==MPI_PROC_NULL){
              switch(bcFaceZright){
              case 1:   // perfect mirror
              MODULO(z,Lz);
     *w= -*w;
     break;
     case 2:   // riemmission
     MODULO(z,Lz);
     double harvest, prob,theta;
    // u
    harvest =   rand()/(double)RAND_MAX;
    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
    harvest =   rand()/(double)RAND_MAX;
    theta = 2.0*M_PI*harvest;
     *u = ut*prob*cos(theta);
    // v
     *v = vt*prob*sin(theta);
    // w
    harvest =   rand()/(double)RAND_MAX;
    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
    harvest =   rand()/(double)RAND_MAX;
    theta = 2.0*M_PI*harvest;
     *w = -fabs(wt*prob*cos(theta));
     break;



     }
     }

     if (*z < 0 && vct->getZleft_neighbor()==MPI_PROC_NULL){
     switch(bcFaceZleft){
     case 1:   // perfect mirror
     MODULO(z,Lz);
     *w= -*w;               
     break;
     case 2:   // riemmission
     MODULO(z,Lz);
     double harvest, prob,theta;
    // u
    harvest =   rand()/(double)RAND_MAX;
    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
    harvest =   rand()/(double)RAND_MAX;
    theta = 2.0*M_PI*harvest;
     *u = ut*prob*cos(theta);
    // v
     *v = vt*prob*sin(theta);
    // w
    harvest =   rand()/(double)RAND_MAX;
    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
    harvest =   rand()/(double)RAND_MAX;
    theta = 2.0*M_PI*harvest;
     *w = fabs(wt*prob*cos(theta));
     break;



     }
     }
     */





}

#endif

