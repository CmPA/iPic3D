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
    void initMaxwellJuttner(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void initMaxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** Initial condition: uniform in space with anisotropic Kappa */
    void initAnisotropicKappa(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** KAW turbulence setup */
    void KAWTurbulencePert(Grid * grid, Field * EMf, VirtualTopology3D * vct, double B0x, double mime, double TiTe, bool symmetric);  

    /** particle injector */
    int injector_rand_box(Grid* grid,VirtualTopology3D* vct, Field* EMf, double x_center_inject, double y_center_inject, double z_center_inject, double L_inject);
    /** particle injector monoenergetic*/
    int injector_rand_box_mono(Grid* grid,VirtualTopology3D* vct, Field* EMf);
    /** Derivative of f0 wrt vpar */
    double df0_dvpar(double vpar, double vperp);
    /** Derivative of f0 wrt vperp */
    double df0_dvperp(double vpar, double vperp);
    /** Equilibrium bi-maxwellian f0 */
    double f0(double vpar, double vperp);
    /** Rotate velocities in plane XY of angle theta */
    void RotatePlaneXY(double theta);
    /** Compute overlap of intervals */
    double interval_overlap(double xd0, double xd1, double xp0, double xp1, double& x_start_overlap);
    /** mover with the esplicit non relativistic scheme */
    void mover_explicit(Grid * grid, VirtualTopology3D * vct, Field * EMf);

    void get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext);
    void get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez);
    void get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]);
    /** mover with a Predictor-Corrector Scheme  woith older code versions*/
    int mover_PC_old(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** mover with a Predictor-Corrector Scheme */
    int mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    int mover_PC_sub(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** mover with a Predictor-Corrector Scheme for 2D cylindrical symmetry*/
    int mover_PC_sub_cyl(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** relativistic mover with a Boris-like scheme */
    int mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf);
    /** relativistic mover with the Celested3D scheme */
    int mover_relativistic_celeste(Grid * grid, VirtualTopology3D * vct, Field * EMf);

    /** particle repopulator */
    int particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is);
    /** particle reflector */
    int particle_reflector(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is);
    /** interpolation Particle->Grid only charge density, current */
    void interpP2G_notP(Field * EMf, Grid * grid, VirtualTopology3D * vct);
    /** interpolation Particle->Grid only for pressure tensor */
    void interpP2G_onlyP(Field * EMf, Grid * grid, VirtualTopology3D * vct);
    /*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
    double deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center);
	 /** Delete the particles outside the sphere with radius R and center x_center y_center */
	 double deleteParticlesOutsideSphere(double R, double x_center, double y_center, double z_center);
	/** Delete the particles outside the cube with dimension L */
	 double deleteParticlesOutsideBox(double L);
	 /** Delete particles in the outer shell towards the center */
	double deleteParticlesOuterFrame(double multx, double multy, double multz);
	/** Reflect particles in the outer shell towards the center */
	double ReturnToCenterOuterFrame(double multx, double multy, double multz);
	/** Reflect particles in outside a circle */
	double ReturnToCenterCircle();
	/** Reflect particles in the outer shell towards the center after regenerating their speed form intial temperature*/
	double ReturnRegeneratedToCenterOuterFrame(double multx, double multy, double multz);
	/** Initial condition: localised in a box and maxwellian in velocity */
    void dual_spark_plug(Grid* grid,Field* EMf, VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center);
	 /** record the Flux through the FluxLoops **/
	 void recordFlux(double oldX, double oldY, double oldZ, double newX, double newY, double newZ, int ptcl);
    /*! Initial condition: given a fluid model (BATSRUS) */
    void MaxwellianFromFluid(Grid* grid,Field* EMf,VirtualTopology3D* vct,Collective *col, int is);
    /*! Initiate dist. func. for a single cell form a fluid model (BATSRUS) */
    void MaxwellianFromFluidCell(Grid* grid, Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, unsigned long* ParticleID);

};


#endif
