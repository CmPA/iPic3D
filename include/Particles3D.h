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
    /** Initial condition: uniform in space and maxwellian in velocity */
    void maxwellian(Grid * grid, VirtualTopology3D * vct);
    void twostream1D(Grid * grid, VirtualTopology3D * vct);

    void get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext);
    void get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez);
    void get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]);

    int mover_relativistic_pos(Grid*, VirtualTopology3D*);
    //int mover_relativistic_mom_EM(Grid*, VirtualTopology3D*, double*** Ex, double*** Ey, double*** Ez);
    int mover_relativistic_mom_ES(Grid*, VirtualTopology3D*, double*** Ex, double*** Ey, double*** Ez);
    //void interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct);
    void mover_relativistic_vel(Grid*, VirtualTopology3D*);
};


#endif
