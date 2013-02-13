/***************************************************************************
  ESfield3D.h  -  Electrostatic field definition
  -------------------
begin             : Wed Jul 14 2004
copyright         : (C) 2004 Los Alamos National Laboratory
developers        : Stefano Markidis, Giovanni Lapenta
email             : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef ESfield3D_H
#define ESfield3D_H


#include <iostream>

#include <math.h>
#include <mpi.h>


#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace.h"
#include "../solvers/CG.h"
#include "../solvers/GMRES.h"

using std::cout;
using std::cerr;
using std::endl;



/**
 *  Electrostatic field and sources defined for each local grid 
 *   
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */

class ESfield3D : public Field{
    public:
        /** constructor */
        ESfield3D(CollectiveIO *col,Grid *grid);
        /** destructor */
        ~ESfield3D();
        /** set to 0 all the densities fields */
        void setZeroDensities();
        /** calculate error in norm 2 for a test case */
        double calculateError(Grid *grid);
        /** Image of Poisson Solver */
        void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct);
        /** Calculate Electric field with the implicit solver */
        void calculateE(Grid *grid, VirtualTopology *vct);
        /** calculate PHI in the electrostatic limit with the Poisson equation */
        void calculatePHI(Grid *grid, VirtualTopology *vct);
        /** communicate ghost for grid -> Particles interpolation */
        void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, 	int bcFaceZleft,VirtualTopology *vct);
        /** add an amount of charge density to charge density field at node X,Y,Z */
        void addRho(double ***weight, int X, int Y, int Z, int ns);
        /** set to 0 all the densities fields */
        void setZeroDensities();
        /** get PHI(X,Y,Z)  */
        double &getPHI(int indexX, int indexY, int indexZ) const;
        /** get Ex(X,Y,Z)  */
        double &getEx(int indexX, int indexY, int indexZ) const;
        /** get Ey(X,Y,Z)  */
        double &getEy(int indexX, int indexY, int indexZ) const;
        /** get Ez(X,Y,Z)  */
        double &getEz(int indexX, int indexY, int indexZ) const;
        /** get Bx(X,Y,Z)  */
        double &getBx(int indexX, int indexY, int indexZ) const;
        /** get By(X,Y,Z)  */
        double &getBy(int indexX, int indexY, int indexZ) const;
        /** get Bz(X,Y,Z)  */
        double &getBz(int indexX, int indexY, int indexZ) const;
        /** get rhoc(X,Y,Z) */
        double*** getRHOC();
        /** get PHI(X,Y,Z) */
        double*** getPHI();
        /** add an amount of charge density to charge density field at node X,Y,Z */
        void addRho(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of current density - direction X to current density field at node X,Y,Z */
        void addJx(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of current density - direction Y to current density field at node X,Y,Z */
        void addJy(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of current density - direction Z to current density field at node X,Y,Z */
        void addJz(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
        void addPxx(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
        void addPxy(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
        void addPxz(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
        void addPyy(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
        void addPyz(double ***weight, int X, int Y, int Z, int ns);
        /** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
        void addPzz(double ***weight, int X, int Y, int Z, int ns);
        /** print electromagnetic fields info */
        void print(void) const;
    private:
        /** free space dielectric */
        double eps0;
        /** magnetic permeability */
        double mu0;
        /** light speed */
        double c;
        /** time step */
        double dt;
        /** decentering parameter */
        double th;
        /** delt = c*th*dt */
        double delt;
        /** number of particles species */
        int ns;
        /** charge to mass ratio array for different species */
        double *qom;


        // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
        /** number of cells - X direction, including + 2 (guard cells) */
        int nxc;
        /** number of nodes - X direction, including + 2 extra nodes for guard cells */
        int nxn;
        /** number of cell - Y direction, including + 2 (guard cells) */
        int nyc;
        /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
        int nyn;
        /** number of cell - Z direction, including + 2 (guard cells) */
        int nzc;
        /** number of nodes - Z direction, including + 2 extra nodes for guard cells */
        int nzn;
        /** local grid boundaries coordinate  */
        double xStart, xEnd, yStart, yEnd, zStart, zEnd;
        /** grid spacing */
        double dx, dy, dz, invVOL;
        /** simulation box length - X direction   */
        double Lx;
        /** simulation box length - Y direction   */
        double Ly;
        /** simulation box length - Z direction   */
        double Lz;

        /** PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes*/
        double ***PHI;
        /** Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
        double ***Ex;
        /** Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
        double ***Ey;
        /** Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
        double ***Ez;
        /*********************************************************************************
        /*************                SOURCES                                           **
        /********************************************************************************/

        /** rhoc: charge density (indexX, indexY, indexZ), defined on central points between nodes */
        double ***rhoc;
        /** rhon: charge density (indexX, indexY, indexZ), defined on nodes */
        double ***rhon;
        /** rhon: charge density (indexX, indexY, indexZ, #species), defined on nodes */
        double ****rhons;


        /** Field Boundary Condition
          0 =  Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain
          1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
          2 =  Periodic condition
          */

        /** Boundary Condition Electrostatic Potential: FaceXright */
        int bcPHIfaceXright;
        /** Boundary Condition Electrostatic Potential:FaceXleft */
        int bcPHIfaceXleft;
        /** Boundary Condition Electrostatic Potential:FaceYright */
        int bcPHIfaceYright;
        /** Boundary Condition Electrostatic Potential:FaceYleft */
        int bcPHIfaceYleft;
        /** Boundary Condition Electrostatic Potential:FaceYright */
        int bcPHIfaceZright;
        /** Boundary Condition Electrostatic Potential:FaceYleft */
        int bcPHIfaceZleft;



};
inline ESfield3D::ESfield3D(CollectiveIO *col,Grid *grid){
    nxc = grid->getNXC();
    nxn = grid->getNXN();
    nyc = grid->getNYC();
    nyn = grid->getNYN();
    nzc = grid->getNZC();
    nzn = grid->getNZN();
    dx = grid->getDX();
    dy = grid ->getDY();
    dz = grid->getDZ();
    invVOL = 1/(dx*dy*dz);
    xStart = grid->getXstart();
    xEnd   = grid->getXend();
    yStart = grid->getYstart();
    yEnd   = grid->getYend();
    zStart = grid->getZstart();
    zEnd   = grid->getZend();
    Lx = col->getLx();
    Ly = col->getLy();
    Lz = col->getLz();
    ns  = col->getNs();
    eps0 = col->getEps0();
    mu0 = col->getMu0();
    c = col->getC();
    dt = col->getDt();
    th = col->getTh();
    delt = c*th*dt;
    qom = new double[ns];
    for (int i=0; i < ns;i++)
        qom[i] = col->getQOM(i);
    // boundary conditions
    bcPHIfaceXright = col->getBcPHIfaceXright();
    bcPHIfaceXleft = col->getBcPHIfaceXleft();
    bcPHIfaceYright = col->getBcPHIfaceYright();
    bcPHIfaceYleft = col->getBcPHIfaceYleft();
    bcPHIfaceZright = col->getBcPHIfaceZright();
    bcPHIfaceZleft = col->getBcPHIfaceZleft();
    // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
    Ex = newArr3(double,nxn,nyn,nzn);
    Ey = newArr3(double,nxn,nyn,nzn);
    Ez = newArr3(double,nxn,nyn,nzn);
    // involving species
    rhons = newArr4(double,nxn,nyn,nzn,ns);
    // arrays allocation: central points ---> the first central point has index 1, the last has index nxn-2!
    PHI    = newArr3(double,nxc,nyc,nzc);
    rhoc   = newArr3(double,nxc,nyc,nzc);

}
/** destructor: deallocate arrays*/
inline ESfield3D::~ESfield3D(){
    // nodes
    delArr3(Ex,nxn,nyn);
    delArr3(Ey,nxn,nyn);
    delArr3(Ez,nxn,nyn);
    // nodes and species
    delArr4(rhons,nxn,nyn,nzn);
    // central points
    delArr3(PHI,nxc,nyc);
    delArr3(rhoc,nxc,nyc);

}
/** set to 0 all the densities fields */
inline  void ESfield3D::setZeroDensities(){
    for (register int i=0; i < nxn*nyn*nzn; i++)
        ***rhon++ = 0.0;
    for (register int i=0; i < nxc*nyc*nzc; i++)
        ***rhoc++ = 0.0;
    for (register int i=0; i < nxn*nyn*nzn*ns; i++){
        ****rhons++ = 0.0;                      

    }
    /** calculate PHI in the electrostatic limit with the Poisson equation */
    inline void ESfield3D::calculatePHI(Grid *grid, VirtualTopology *vct){
        // solve the equation Laplacian PHI = drho
        double *xkrylov = new double[(nxc-2)*(nyc-2)*(nzc-2)];
        double *sorg    = new double[(nxc-2)*(nyc-2)*(nzc-2)];
        // know term
        phys2solver(sorg,rhoc,nxc,nyc,nzc);
        CG(xkrylov,(nxc-2)*(nyc-2)*(nzc-2),sorg,100,10E-4,&Field::PoissonImage,grid,vct,this);
        //GMRES(xkrylov,(nxc-2)*(nyc-2)*(nzc-2), sorg, 100, 10E-4, 10,&Field::PoissonImage,grid,vct,this);
        solver2phys(PHI,xkrylov,nxc,nyc,nzc);

    }
    /** Image of Poisson Solver */
    inline void ESfield3D::PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct){
        // allocate  2 three dimensional service vectors
        double ***temp = newArr3(double,nxc,nyc,nzc);
        double ***im   = newArr3(double,nxc,nyc,nzc);
        // move from krylov space to physical space and communicate ghost cells
        solver2phys(temp,vector,nxc,nyc,nzc);
        communicateGhost(nxc,nyc,nzc,temp,0,0,0,0,0,0,dx,dy,dz,vct);
        // calculate the laplacian
        grid->lapC2C(im,temp,vct);
        // move from physical space to krylov space, in this case ghost cells don't count
        phys2solver(image,im,nxc,nyc,nzc);
        // deallocate temporary array and objects
        delArr3(temp,nxc,nyc);
        delArr3(im,nxc,nyc);


    }
    /** communicate ghost for grid -> Particles interpolation */
    inline void ESfield3D::communicateGhostP2G(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft,VirtualTopology *vct){
        communicateInterp(nxn,nyn,nzn,rhon,0,0,0,0,0,0,dx,dy,dz,vct);
        communicateInterp(nxn,nyn,nzn,ns,rhons,0,0,0,0,0,0,dx,dy,dz,vct);
    }
    /** add an amount of charge density to charge density field at node X,Y,Z */
    inline void ESfield3D::addRho(double ***weight, int X, int Y, int Z, int ns){
        for (int i=0; i < 2; i++)
            for (int j=0; j < 2; j++)
                for(int k=0; k < 2; k++){
                    rhon[X +(-1 + i)][Y+(-1 + j)][Z + (-1 +k)] += weight[i][j][k]*invVOL;
                    rhons[X +(-1 + i)][Y+(-1 + j)][Z+(-1 +k)][ns] += weight[i][j][k]*invVOL;
                }
    }
    /** get PHI(X,Y,Z)  */
    inline double &ESfield3D::getPHI(int indexX, int indexY, int indexZ) const{
        return(PHI[indexX][indexY][indexZ]);
    }
    /** get Ex(X,Y,Z)  */
    inline double &ESfield3D::getEx(int indexX, int indexY, int indexZ) const{
        return(Ex[indexX][indexY][indexZ]);
    }
    /** get Ey(X,Y,Z)  */
    inline double &ESfield3D::getEy(int indexX, int indexY, int indexZ) const{
        return(Ey[indexX][indexY][indexZ]);
    }
    /** get Ez(X,Y,Z)  */
    inline double &ESfield3D::getEz(int indexX, int indexY, int indexZ) const{
        return(0.0);
    }
    /** get Bx(X,Y,Z)  */
    inline double &ESfield3D::getBx(int indexX, int indexY, int indexZ) const{
        return(0.0;
                }
                /**  get By(X,Y,Z) */
                inline double &ESfield3D::getBy(int indexX, int indexY, int indexZ) const{
                return(0.0);
                }
                /**  get Bz(X,Y,Z) */
                inline double &ESfield3D::getBz(int indexX, int indexY, int indexZ) const{
                return(0.0);
                }
                /** get rhoc(X,Y,Z) */
                inline double*** ESfield3D::getRHOC(){
                return(rhoc);
                }
                /** get PHI(X,Y,Z) */
                inline double*** ESfield3D::getPHI(){
                return(PHI);
                }
#endif
