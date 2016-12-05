/*******************************************************************************************
  Particles3D.h  -  Class for particles of the same species, in a 2D space and 3 component velocity
  -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
 ********************************************************************************************/

#ifndef Part2D_H
#define Part2D_H

#include "Particles3Dcomm.h"
#include "TimeTasks.h"

/**
 * 
 * Class for particles of the same species, in a 2D space and 3 component velocity
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3D:public Particles3Dcomm {
  public:
    /** constructor */
    Particles3D();
    /** destructor */
    ~Particles3D();
    /** Initial condition: uniform in space and motionless */
    void uniform_background(Grid * grid, Field * EMf);
    /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
      <ul>
      <li> dim = 0 --> constant velocity on X direction </li>
      <li> dim = 1 --> constant velocity on Y direction </li>
      </ul>
      */
    void constantVelocity(double vel, int dim, Grid * grid, Field * EMf);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space and maxwellian with reversal across Y in velocity */
    void maxwellian_reversed(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: non uniform in space with maxwellian in thermal equilibrium outside a domain */
    void maxwellian_whistler(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space with Kappa */
    void kappa(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void MaxwellianFromFields(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Force Free initialization (JxB=0) for particles */
    void force_free(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void alt_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Linear_perturbation */
    void linear_perturbation(double deltaBX, double kx, double ky, double theta, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /**Add a periodic perturbation in velocity exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
    void AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid * grid);
    /** Linear delta f for bi-maxwellian plasma */
    double delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_ampl, double Ex_phase, double Ey_ampl, double Ey_phase, double Ez_ampl, double Ez_phase, double theta, Field * EMf);
    /** Derivative of f0 wrt vpar */
    double df0_dvpar(double vpar, double vperp);
    /** Derivative of f0 wrt vperp */
    double df0_dvperp(double vpar, double vperp);
    /** Equilibrium bi-maxwellian f0 */
    double f0(double vpar, double vperp);
    /** Rotate velocities in plane XY of angle theta */
    void RotatePlaneXY(double theta);
    /** mover with the esplicit non relativistic scheme */
    void mover_explicit(Grid * grid, VirtualTopology3D * vct, Field * EMf);

    void get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext);
    void get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez);
    void get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]);
    /** mover with a Predictor-Corrector Scheme */
    int mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    int mover_PC_sub(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** relativistic mover with a Predictor-Corrector scheme */
    int mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** particle repopulator */
    int particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is);
    /** interpolation Particle->Grid only charge density, current */
    void interpP2G_notP(Field * EMf, Grid * grid, VirtualTopology3D * vct);
    /** interpolation Particle->Grid only for pressure tensor */
    void interpP2G_onlyP(Field * EMf, Grid * grid, VirtualTopology3D * vct);
    /*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
    double deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center);

    /*! Initial condition: given a fluid model (BATSRUS) */
    void MaxwellianFromFluid(Grid* grid,Field* EMf,VirtualTopology3D* vct,Collective *col, int is);
    /*! Initiate dist. func. for a single cell form a fluid model (BATSRUS) */
    void MaxwellianFromFluidCell(Grid* grid, Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, unsigned long* ParticleID);

};


#endif
