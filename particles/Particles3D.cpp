/*******************************************************************************************
  Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/


#include <mpi.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include "asserts.h"

#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "CollectiveIO.h"
#include "Collective.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "Field.h"
#include "MPIdata.h"
#include "ipicdefs.h"
#include "TimeTasks.h"

#include "Particles3D.h"

#include "debug.h"
#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
// particles processed together
#define P_SAME_TIME 2

/**
 * 
 * Class for particles of the same species
 * @date Fri Jun 4 2009
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles3D::Particles3D() {
  // see allocate(int species, CollectiveIO* col, VirtualTopology3D* vct, Grid* grid)

}
/** deallocate particles */
Particles3D::~Particles3D() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
}

/** particles are uniformly distributed with zero velocity   */
void Particles3D::uniform_background(Grid * grid, Field * EMf) {
  int counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              u[counter] = 0.0;
              v[counter] = 0.0;
              w[counter] = 0.0;
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
              if (TrackParticleID)
                ParticleID[counter] = counter * (long long) pow(10.0, BirthRank[1]) + BirthRank[0];
              counter++;
            }


  cout << "Velocity Maxwellian Distribution " << endl;
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
  <ul>
  <li> dim = 0 --> constant velocity on X direction </li>
  <li> dim = 1 --> constant velocity on Y direction </li>
  <li> dim = 2 --> constant velocity on Z direction </li>
  </ul>

*/
void Particles3D::constantVelocity(double vel, int dim, Grid * grid, Field * EMf) {
  switch (dim) {
    case 0:
      for (int i = 0; i < nop; i++)
        u[i] = vel, v[i] = 0.0, w[i] = 0.0;
      break;
    case 1:
      for (int i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = vel, w[i] = 0.0;
      break;
    case 2:
      for (int i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = 0.0, w[i] = vel;
      break;

  }

}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles3D::alt_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {
}

#ifdef BATSRUS
/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::MaxwellianFromFluid(Grid* grid,Field* EMf,VirtualTopology3D* vct,Collective *col, int is){

  /*
   * Constuctiong the distrebution function from a Fluid model
   */

  // loop over grid cells and set position, velociy and charge of all particles indexed by counter
  // there are multiple (27 or so) particles per grid cell.
  int i,j,k,counter=0;
  for (i=1; i< grid->getNXC()-1;i++)
    for (j=1; j< grid->getNYC()-1;j++)
      for (k=1; k< grid->getNZC()-1;k++)
        MaxwellianFromFluidCell(grid,col,is, i,j,k,counter,x,y,z,q,u,v,w,ParticleID);
}

void Particles3D::MaxwellianFromFluidCell(Grid* grid, Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, long long* ParticleID)
{
  /*
   * grid           : local grid object (in)
   * col            : collective (global) object (in)
   * is             : species index (in)
   * i,j,k          : grid cell index on proc (in)
   * ip             : particle number counter (inout)
   * x,y,z          : particle position (out)
   * q              : particle charge (out)
   * vx,vy,vz       : particle velocity (out)
   * ParticleID     : particle tracking ID (out)
   */

  double harvest;
  double prob, theta;

  // loop over particles inside grid cell i,j,k
  for (int ii=0; ii < npcelx; ii++)
    for (int jj=0; jj < npcely; jj++)
      for (int kk=0; kk < npcelz; kk++){
        // Assign particle positions: uniformly spaced. x_cellnode + dx_particle*(0.5+index_particle)
        x[ip] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);
        y[ip] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
        z[ip] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
        // q = charge
        q[ip] =  (qom/fabs(qom))*(col->getFluidRhoCenter(i,j,k,is)/npcel)*(1.0/grid->getInvVOL());
        // u = X velocity
        harvest =   rand()/(double)RAND_MAX;
        prob  = sqrt(-2.0*log(1.0-.999999*harvest));
        harvest =   rand()/(double)RAND_MAX;
        theta = 2.0*M_PI*harvest;
        u[ip] = col->getFluidUx(i,j,k,is) + col->getFluidUthx(i,j,k,is)*prob*cos(theta);
        // v = Y velocity
        v[ip] = col->getFluidUy(i,j,k,is) + col->getFluidUthy(i,j,k,is)*prob*sin(theta);
        // w = Z velocity
        harvest =   rand()/(double)RAND_MAX;
        prob  = sqrt(-2.0*log(1.0-.999999*harvest));
        harvest =   rand()/(double)RAND_MAX;
        theta = 2.0*M_PI*harvest;
        w[ip] = col->getFluidUz(i,j,k,is) + col->getFluidUthz(i,j,k,is)*prob*cos(theta);
        if (TrackParticleID)
          ParticleID[ip]= ip*(long long)pow(10.0,BirthRank[1])+BirthRank[0];
        ip++ ;
      }
}
#endif

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  double harvest;
  double prob, theta, sign;
  int counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);  // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = u0 + uth * prob * cos(theta);
              // v
              v[counter] = v0 + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = w0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (long long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}

/** Force Free initialization (JxB=0) for particles */
void Particles3D::force_free(Grid * grid, Field * EMf, VirtualTopology3D * vct) {


  double harvest, prob, theta;
  int counter = 0;
  double shaperx, shapery, shaperz;
  double flvx = 1.0, flvy = 1.0, flvz = 1.0;


  /* initialize random generator */
  srand(vct->getCartesian_rank() + 1 + ns);
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              flvx = 1.0;
              flvy = 1.0;
              flvz = 1.0;
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / invVOL);
              shaperx = tanh((y[counter] - Ly / 2) / delta) / cosh((y[counter] - Ly / 2) / delta) / delta;
              shaperz = 1.0 / (cosh((y[counter] - Ly / 2) / delta) * cosh((y[counter] - Ly / 2) / delta)) / delta;
              shapery = shapery;
              // new drift velocity to satisfy JxB=0
              flvx = u0 * flvx * shaperx;
              flvz = w0 * flvz * shaperz;
              flvy = v0 * flvy * shapery;
              u[counter] = c;
              v[counter] = c;
              w[counter] = c;
              while ((fabs(u[counter]) >= c) | (fabs(v[counter]) >= c) | (fabs(w[counter]) >= c)) {
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                u[counter] = flvx + uth * prob * cos(theta);
                // v
                v[counter] = flvy + vth * prob * sin(theta);
                // w
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                w[counter] = flvz + wth * prob * cos(theta);
              }
              if (TrackParticleID)
                ParticleID[counter] = counter * (long long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }

}

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid * grid) {

  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);
  jx_mod *= alpha;
  jy_mod *= alpha;
  jz_mod *= alpha;
  for (int i = 0; i < nop; i++) {
    u[i] += jx_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jx_phase);
    v[i] += jy_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jy_phase);
    w[i] += jz_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jz_phase);
  }
}


/** explicit mover */
void Particles3D::mover_explicit(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  // to be implemented

}
/** mover with a Predictor-Corrector scheme */
void Particles3D::mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  #pragma omp master
  { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
  const pfloat dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
  #pragma omp for schedule(static)
  // why does single precision make no difference in execution speed?
  //#pragma simd vectorlength(VECTOR_WIDTH)
  for (int pidx = 0; pidx < nop; pidx++) {
    // copy the particle
    const pfloat xorig = x[pidx];
    const pfloat yorig = y[pidx];
    const pfloat zorig = z[pidx];
    const pfloat uorig = u[pidx];
    const pfloat vorig = v[pidx];
    const pfloat worig = w[pidx];
    pfloat xavg = xorig;
    pfloat yavg = yorig;
    pfloat zavg = zorig;
    pfloat uavg;
    pfloat vavg;
    pfloat wavg;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {
      // interpolation G-->P
      const pfloat ixd = floor((xavg - xstart) * inv_dx);
      const pfloat iyd = floor((yavg - ystart) * inv_dy);
      const pfloat izd = floor((zavg - zstart) * inv_dz);
      // interface of index to right of cell
      int ix = 2 + int(ixd);
      int iy = 2 + int(iyd);
      int iz = 2 + int(izd);

      // use field data of closest cell in domain
      //
      if (ix < 1) ix = 1;
      if (iy < 1) iy = 1;
      if (iz < 1) iz = 1;
      if (ix > nxc) ix = nxc;
      if (iy > nyc) iy = nyc;
      if (iz > nzc) iz = nzc;
      // index of cell of particle;
      //const int cx = ix - 1;
      //const int cy = iy - 1;
      //const int cz = iz - 1;

      const pfloat xi0   = xavg - grid->get_pfloat_XN(ix-1);
      const pfloat eta0  = yavg - grid->get_pfloat_YN(iy-1);
      const pfloat zeta0 = zavg - grid->get_pfloat_ZN(iz-1);
      const pfloat xi1   = grid->get_pfloat_XN(ix) - xavg;
      const pfloat eta1  = grid->get_pfloat_YN(iy) - yavg;
      const pfloat zeta1 = grid->get_pfloat_ZN(iz) - zavg;

      pfloat Exl = 0.0;
      pfloat Eyl = 0.0;
      pfloat Ezl = 0.0;
      pfloat Bxl = 0.0;
      pfloat Byl = 0.0;
      pfloat Bzl = 0.0;

      pfloat weights[8];
      const pfloat weight0 = invVOL*xi0;
      const pfloat weight1 = invVOL*xi1;
      const pfloat weight00 = weight0*eta0;
      const pfloat weight01 = weight0*eta1;
      const pfloat weight10 = weight1*eta0;
      const pfloat weight11 = weight1*eta1;
      weights[0] = weight00*zeta0; // weight000
      weights[1] = weight00*zeta1; // weight001
      weights[2] = weight01*zeta0; // weight010
      weights[3] = weight01*zeta1; // weight011
      weights[4] = weight10*zeta0; // weight100
      weights[5] = weight10*zeta1; // weight101
      weights[6] = weight11*zeta0; // weight110
      weights[7] = weight11*zeta1; // weight111
      //weights[0] = xi0 * eta0 * zeta0 * qi * invVOL; // weight000
      //weights[1] = xi0 * eta0 * zeta1 * qi * invVOL; // weight001
      //weights[2] = xi0 * eta1 * zeta0 * qi * invVOL; // weight010
      //weights[3] = xi0 * eta1 * zeta1 * qi * invVOL; // weight011
      //weights[4] = xi1 * eta0 * zeta0 * qi * invVOL; // weight100
      //weights[5] = xi1 * eta0 * zeta1 * qi * invVOL; // weight101
      //weights[6] = xi1 * eta1 * zeta0 * qi * invVOL; // weight110
      //weights[7] = xi1 * eta1 * zeta1 * qi * invVOL; // weight111

      // creating these aliases seems to accelerate this method by about 30%
      // on the Xeon host, processor, suggesting deficiency in the optimizer.
      //
      arr1_pfloat_get field_components[8];
      field_components[0] = fieldForPcls[ix  ][iy  ][iz  ]; // field000
      field_components[1] = fieldForPcls[ix  ][iy  ][iz-1]; // field001
      field_components[2] = fieldForPcls[ix  ][iy-1][iz  ]; // field010
      field_components[3] = fieldForPcls[ix  ][iy-1][iz-1]; // field011
      field_components[4] = fieldForPcls[ix-1][iy  ][iz  ]; // field100
      field_components[5] = fieldForPcls[ix-1][iy  ][iz-1]; // field101
      field_components[6] = fieldForPcls[ix-1][iy-1][iz  ]; // field110
      field_components[7] = fieldForPcls[ix-1][iy-1][iz-1]; // field111

      for(int c=0; c<8; c++)
      {
        Bxl += weights[c] * field_components[c][0];
        Byl += weights[c] * field_components[c][1];
        Bzl += weights[c] * field_components[c][2];
        Exl += weights[c] * field_components[c][3];
        Eyl += weights[c] * field_components[c][4];
        Ezl += weights[c] * field_components[c][5];
      }
      const double Omx = qdto2mc*Bxl;
      const double Omy = qdto2mc*Byl;
      const double Omz = qdto2mc*Bzl;

      // end interpolation
      const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const pfloat denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const pfloat ut = uorig + qdto2mc * Exl;
      const pfloat vt = vorig + qdto2mc * Eyl;
      const pfloat wt = worig + qdto2mc * Ezl;
      //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation 
      uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
      // update average position
      xavg = xorig + uavg * dto2;
      yavg = yorig + vavg * dto2;
      zavg = zorig + wavg * dto2;
    }                           // end of iteration
    // update the final position and velocity
    x[pidx] = xorig + uavg * dt;
    y[pidx] = yorig + vavg * dt;
    z[pidx] = zorig + wavg * dt;
    u[pidx] = 2.0 * uavg - uorig;
    v[pidx] = 2.0 * vavg - vorig;
    w[pidx] = 2.0 * wavg - worig;
  }                             // END OF ALL THE PARTICLES
  #pragma omp master
  { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
}

/** mover with a Predictor-Corrector scheme */
void Particles3D::mover_PC_vectorized(
  Grid * grid, VirtualTopology3D * vct, Field * EMf)
{
  assert_eq(nxc,nxn-1);
  assert_eq(nyc,nyn-1);
  assert_eq(nzc,nzn-1);
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  // initialize average positions
  #pragma omp for schedule(static)
  for(int pidx = 0; pidx < nop; pidx++)
  {
    xavg[pidx] = x[pidx];
    yavg[pidx] = y[pidx];
    zavg[pidx] = z[pidx];
  }

  const pfloat dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
  for(int niter=1; niter<=NiterMover; niter++)
  {
    // sort particles based on the time-averaged position
    if(niter>1) // on first iteration already was sorted to sum moments
    {
      #pragma omp master
      {
        timeTasks_begin_task(TimeTasks::MOVER_PCL_SORTING);
        sort_particles_serial(xavg, yavg, zavg, grid,vct);
        timeTasks_end_task(TimeTasks::MOVER_PCL_SORTING);
      }
      #pragma omp barrier
    }

    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
    // move particles in parallel
    //
    // iterate over mesh cells
    //const int ncells=nxc*nyc*nzc;
    //int *numpcls_in_bucket_1d = &numpcls_in_bucket[0][0][0];
    //int *bucket_offset_1d = &bucket_offset[0][0][0];
    int serial_pidx = 0;
    #pragma omp for collapse(2) // schedule(static)
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    //for(int cell=0; cell<ncells; cell++)
    {
      // interface to the right of cell
      const int ix = cx+1;
      const int iy = cy+1;
      const int iz = cz+1;

      arr1_pfloat_get field_components[8];
      field_components[0] = fieldForPcls[ix][iy][iz]; // field000
      field_components[1] = fieldForPcls[ix][iy][cz]; // field001
      field_components[2] = fieldForPcls[ix][cy][iz]; // field010
      field_components[3] = fieldForPcls[ix][cy][cz]; // field011
      field_components[4] = fieldForPcls[cx][iy][iz]; // field100
      field_components[5] = fieldForPcls[cx][iy][cz]; // field101
      field_components[6] = fieldForPcls[cx][cy][iz]; // field110
      field_components[7] = fieldForPcls[cx][cy][cz]; // field111

      // push all particles in mesh cell
      //
      //const int numpcls_in_cell = numpcls_in_bucket_1d[cell];
      const int numpcls_in_cell = get_numpcls_in_bucket(cx,cy,cz);
      const int bucket_offset = get_bucket_offset(cx,cy,cz);
      const int bucket_end = bucket_offset+numpcls_in_cell;
      // this should vectorize, but could be faster if particle
      // data were aligned.
      #pragma simd
      //for(int pidx=bucket_offset_1d[cell]; pidx<numpcls_in_cell; pidx++)
      for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
      {
        // serial case: check that pidx is correct
        //assert_eq(pidx,serial_pidx++);
        // confirm that particle is in correct cell
        if(true)
        {
          int cx_,cy_,cz_;
          get_safe_cell_for_pos(cx_,cy_,cz_,xavg[pidx],yavg[pidx],zavg[pidx]);
          //if((cx_!=cx)
          // ||(cy_!=cy)
          // ||(cz_!=cz))
          //{
          //  dprintf("\n\t cx =%d, cy =%d, cz =%d"
          //          "\n\t cx_=%d, cy_=%d, cz_=%d"
          //          "\n\t x=%g, y=%g, z_=%g",
          //          cx,cy,cz,
          //          cx_,cy_,cz_,
          //          xavg[pidx], yavg[pidx], zavg[pidx]);
          //}
          assert_eq(cx_,cx);
          assert_eq(cy_,cy);
          assert_eq(cz_,cz);
        }

        // copy the particle
        const pfloat xorig = x[pidx];
        const pfloat yorig = y[pidx];
        const pfloat zorig = z[pidx];
        const pfloat uorig = u[pidx];
        const pfloat vorig = v[pidx];
        const pfloat worig = w[pidx];

        // compute weights for field components
        //
        double weights[8];
        const double abs_xpos = xavg[pidx];
        const double abs_ypos = yavg[pidx];
        const double abs_zpos = zavg[pidx];
        // xstart marks start of domain excluding ghosts
        const double rel_xpos = abs_xpos - xstart;
        const double rel_ypos = abs_ypos - ystart;
        const double rel_zpos = abs_zpos - zstart;
        // cell position minus 1 (due to ghost cells)
        const double cxm1_pos = rel_xpos * inv_dx;
        const double cym1_pos = rel_ypos * inv_dy;
        const double czm1_pos = rel_zpos * inv_dz;
        // index of interface to right of cell
        const int ix = cx + 1;
        const int iy = cy + 1;
        const int iz = cz + 1;
        // fraction of the distance from the right of the cell
        const double w1x = cx - cxm1_pos;
        const double w1y = cy - cym1_pos;
        const double w1z = cz - czm1_pos;
        // fraction of distance from the left
        const double w0x = 1-w1x;
        const double w0y = 1-w1y;
        const double w0z = 1-w1z;
        const double weight00 = w0x*w0y;
        const double weight01 = w0x*w1y;
        const double weight10 = w1x*w0y;
        const double weight11 = w1x*w1y;
        weights[0] = weight00*w0z; // weight000
        weights[1] = weight00*w1z; // weight001
        weights[2] = weight01*w0z; // weight010
        weights[3] = weight01*w1z; // weight011
        weights[4] = weight10*w0z; // weight100
        weights[5] = weight10*w1z; // weight101
        weights[6] = weight11*w0z; // weight110
        weights[7] = weight11*w1z; // weight111

        pfloat Exl = 0.0;
        pfloat Eyl = 0.0;
        pfloat Ezl = 0.0;
        pfloat Bxl = 0.0;
        pfloat Byl = 0.0;
        pfloat Bzl = 0.0;
        for(int c=0; c<8; c++)
        {
          Bxl += weights[c] * field_components[c][0];
          Byl += weights[c] * field_components[c][1];
          Bzl += weights[c] * field_components[c][2];
          Exl += weights[c] * field_components[c][3];
          Eyl += weights[c] * field_components[c][4];
          Ezl += weights[c] * field_components[c][5];
        }
        const double Omx = qdto2mc*Bxl;
        const double Omy = qdto2mc*Byl;
        const double Omz = qdto2mc*Bzl;

        // end interpolation
        const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
        const pfloat denom = 1.0 / (1.0 + omsq);
        // solve the position equation
        const pfloat ut = uorig + qdto2mc * Exl;
        const pfloat vt = vorig + qdto2mc * Eyl;
        const pfloat wt = worig + qdto2mc * Ezl;
        //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
        const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
        // solve the velocity equation 
        const pfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
        const pfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
        const pfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
        // update average position
        xavg[pidx] = xorig + uavg * dto2;
        yavg[pidx] = yorig + vavg * dto2;
        zavg[pidx] = zorig + wavg * dto2;

        // if it is the last iteration, update the position and velocity
        // (hopefully this will not compromise vectorization...)
        if(niter==NiterMover)
        {
          x[pidx] = xorig + uavg * dt;
          y[pidx] = yorig + vavg * dt;
          z[pidx] = zorig + wavg * dt;
          u[pidx] = 2.0 * uavg - uorig;
          v[pidx] = 2.0 * vavg - vorig;
          w[pidx] = 2.0 * wavg - worig;
        }
      }
    }
    #pragma omp master
    { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
  }
}

/** communicate particle after moving them */
int Particles3D::communicate_particles(VirtualTopology3D * vct)
{
  timeTasks_set_communicating(); // communicating until end of scope
  const int avail = communicate(vct);
  if (avail < 0)
    return (-1);
  MPI_Barrier(MPI_COMM_WORLD);
  // communicate again if particles are not in the correct domain
  while (isMessagingDone(vct) > 0) {
    // COMMUNICATION
    const int avail = communicate(vct);
    if (avail < 0)
      return (-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return 0; // exit successfully
}

/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  return (0);
}

int Particles3D::particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf){

  if (vct->getCartesian_rank()==0){
    cout << "*** Repopulator species " << ns << " ***" << endl;
  }
  double  FourPI =16*atan(1.0);
  int avail;
  int store_nop=nop;

  ////////////////////////
  // INJECTION FROM XLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft == 2){ // use Field topology in this case
    int particles_index=0;
    int nplast = nop-1;

    while (particles_index < nplast+1) {
      if (x[particles_index] < 3.0*dx ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }

    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< 4;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];


                particles_index++ ;
              }
    nop = particles_index;
  }

  store_nop = nop;

  ////////////////////////
  // INJECTION FROM YLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getYleft_neighbor() == MPI_PROC_NULL  && bcPfaceYleft == 2)
  {
    int particles_index=0;
    int nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] < 3.0*dy ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< 4;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getZleft_neighbor() == MPI_PROC_NULL  && bcPfaceZleft == 2)
  {
    int particles_index=0;
    int nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] < 3.0*dz ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< 4;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM XRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getXright_neighbor() == MPI_PROC_NULL  && bcPfaceXright == 2){
    int particles_index=0;
    int nplast = nop-1;
    while (particles_index < nplast+1) {
      if (x[particles_index] > (Lx-3.0*dx) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=(grid->getNXC()-4); i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM YRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getYright_neighbor() == MPI_PROC_NULL  && bcPfaceYright == 2)
  {
    int particles_index=0;
    int nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] > (Ly-3.0*dy) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=(grid->getNYC()-4); j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getZright_neighbor() == MPI_PROC_NULL  && bcPfaceZright == 2)
  {
    int particles_index=0;
    int nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] > (Lz-3.0*dz) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=(grid->getNZC()-4); k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                q[particles_index] =  (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(long long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  if (vct->getCartesian_rank()==0){
    cout << "*** number of particles " << nop << " ***" << endl;
  }

  //********************//
  // COMMUNICATION
  // *******************//
  avail = communicate(vct);
  if (avail < 0) return(-1);

  MPI_Barrier(MPI_COMM_WORLD);

  // communicate again if particles are not in the correct domain
  while(isMessagingDone(vct) >0){
    // COMMUNICATION
    avail = communicate(vct);
    if (avail < 0)
      return(-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return(0); // exit succcesfully (hopefully)
}

/** apply a linear perturbation to particle distribution */
void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  double value1 = 0.0, value2 = 0.0, max_value = 0.0, min_value = 0.0, phi, n;
  int counter = 0, total_generated = 0;
  bool rejected;
  double harvest, prob, theta;
  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  double integral = 0.0;

  alpha = deltaBoB * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0) + EMf->getBz(1, 1, 0) * EMf->getBz(1, 1, 0)) / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;



  // find the maximum value of f=1+delta_f/f0
  for (register double vpar = -2 * uth; vpar <= 2 * uth; vpar += 0.0005)
    for (register double vperp = 1e-10; vperp <= 2 * vth; vperp += 0.0005)
      for (register double X = xstart; X <= xend; X += 2 * grid->getDX())
        for (register double Y = ystart; Y <= yend; Y += 2 * grid->getDY()) {
          value1 = 1 + delta_f(vpar, vperp, 0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(vpar, vperp);

          if (value1 > max_value)
            max_value = value1;


        }



  max_value *= 3.2;
  phi = 1.48409;
  n = 2.948687;                 // security factor...
  if (ns == 1) {
    max_value *= 3.0;
    phi = -1.65858;
    n = 2.917946;
  }                             // security factor...
  cout << "max-value=" << max_value << " min-value=" << min_value << endl;

  /* initialize random generator */
  srand(vct->getCartesian_rank() + 2);

  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int ii = 0; ii < npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))); ii++)
        for (int jj = 0; jj < npcely; jj++) {
          x[counter] = (ii + .5) * (dx / (npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))))) + grid->getXN(i, j, 0);
          y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, 0);
          q[counter] = (qom / fabs(qom)) * ((0.19635) / npcel) * (1.0 / invVOL);

          // apply rejection method in velocity space
          rejected = true;
          while (rejected) {
            total_generated++;
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            // u
            u[counter] = u0 + uth * prob * cos(theta);
            // v
            v[counter] = v0 + vth * prob * sin(theta);
            // w
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            w[counter] = w0 + wth * prob * cos(theta);

            // test: if rand < (1+delta_f/f0)/max_value --> accepted
            if (rand() / (double) RAND_MAX <= (1 + delta_f(u[counter], v[counter], w[counter], x[counter], y[counter], kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(u[counter], sqrt(v[counter] * v[counter] + w[counter] * w[counter]))) / max_value)
              rejected = false;

          }
          if (TrackParticleID)
            ParticleID[counter] = counter * (long long) pow(10.0, BirthRank[1]) + BirthRank[0];
          counter++;
        }
  nop = counter + 1;
  // if (vct->getCartesian_rank()==0)
  cout << "Rejection method: " << (counter + 1) / double (total_generated) * 100 << " % of particles are accepted for species " << ns << " counter=" << counter << endl;
}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double theta, Field * EMf) {
  const complex < double >I(0.0, 1.0);
  const double vperp = sqrt(v * v + w * w);
  const double vpar = u;
  const double kpar = kx;
  double kperp;
  if (ky == 0.0)                // because this formula is not valid for exactly parallel
    kperp = 1e-9;
  else
    kperp = ky;
  const double om_c = qom / c * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0)) / 2 / M_PI;
  const double phi = atan2(w, v);
  const double lambda = kperp * vperp / om_c;
  const complex < double >omega(omega_re, omega_i);

  const int lmax = 5;           // sum from -lmax to lmax

  double bessel_Jn_array[lmax + 2];
  double bessel_Jn_prime_array[lmax + 1];
  complex < double >a1[2 * lmax + 1], a2[2 * lmax + 1], a3[2 * lmax + 1];
  complex < double >factor, deltaf;

  // rotation of x,y
  double temp;
  temp = x;
  x = x * cos(theta) - y * sin(theta);
  y = temp * sin(theta) + y * cos(theta);


  /** for compilation issues comment this part: PUT in the math stuff */
  // calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
  factor = (kpar * vperp / omega * df0_dvpar(vpar, vperp) + (1.0 - (kpar * vpar / omega)) * df0_dvperp(vpar, vperp));
  for (register int l = -lmax; l < 0; l++) {  // negative index
    a1[l + lmax] = factor / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * 0.5 * pow(-1.0, -l) * (bessel_Jn_array[-l - 1] - bessel_Jn_array[-l + 1]);
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * pow(-1.0, -l) * bessel_Jn_array[-l];
  }

  for (register int l = 0; l < lmax + 1; l++) { // positive index
    a1[l + lmax] = factor / lambda * bessel_Jn_array[l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * bessel_Jn_prime_array[l];
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * bessel_Jn_array[l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * bessel_Jn_array[l];
  }

  deltaf = (0.0, 0.0);
  for (register int l = -lmax; l < lmax + 1; l++) {
    deltaf += (a3[l + lmax] * Ex_mod * exp(I * Ex_phase) + a1[l + lmax] * Ey_mod * exp(I * Ey_phase) + a2[l + lmax] * Ez_mod * exp(I * Ez_phase)) / (kpar * vpar + l * om_c - omega) * exp(-I * phi * (double) l);
  }
  deltaf *= I * qom * exp(I * lambda * sin(phi)) * exp(I * (2 * M_PI * kx * x + 2 * M_PI * ky * y));

  return (real(deltaf));
}

double Particles3D::df0_dvpar(double vpar, double vperp) {
  double result;
  result = -2 * (vpar - u0) / uth / uth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::df0_dvperp(double vpar, double vperp) {
  double result;
  result = -2 * (vperp) / vth / vth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::f0(double vpar, double vperp) {
  double result;
  result = exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

void Particles3D::RotatePlaneXY(double theta) {
  double temp, temp2;
  for (register int s = 0; s < nop; s++) {
    temp = u[s];
    temp2 = v[s];
    u[s] = temp * cos(theta) + v[s] * sin(theta);
    v[s] = -temp * sin(theta) + temp2 * cos(theta);
  }
}

/*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
double Particles3D::deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center){

  int np_current = 0;
  int nplast     = nop-1;

  while (np_current < nplast+1){

    double xd = x[np_current] - x_center;
    double yd = y[np_current] - y_center;
    double zd = z[np_current] - z_center;

    if ( (xd*xd+yd*yd+zd*zd) < R*R ){
      Q_removed += q[np_current];
      del_pack(np_current,&nplast);

    } else {
      np_current++;
    }
  }
  nop = nplast +1;
  return(Q_removed);
}

