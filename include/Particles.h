/*******************************************************************************************
  Particles.h  -  Abstract class for particles of the same species on one processor
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef Particles_H
#define Particles_H

#include "CollectiveIO.h"
#include "VirtualTopology3D.h"
#include "Grid.h"
#include "Field.h"

/**
 * 
 * Abstract class for particles of the same species
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

#include "Particles3Dcomm.h"
//class Particles {
//public:
//  /** allocate particles */
//  virtual void allocate(int species, CollectiveIO * col, VirtualTopology3D * vct, Grid * grid) = 0;
//  /** interpolation Particle -> grid */
//  virtual void interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct) = 0;
//
//
//  /** get X-position array for all the particles */
//  virtual double *getXall() const = 0;
//  /** get Y-position array for all the particles */
//  virtual double *getYall() const = 0;
//  /** get Z-position array for all the particles */
//  virtual double *getZall() const = 0;
//  /** get u (X-velocity) array for all the particles */
//  virtual double *getUall() const = 0;
//  /** get v (Y-velocity) array for all the particles */
//  virtual double *getVall() const = 0;
//  /** get w (Z-velocity) array for all the particles */
//  virtual double *getWall() const = 0;
//  /** get ID array for all the particles */
//  virtual long long *getParticleIDall() const = 0;
//  /**get charge of particle array */
//  virtual double *getQall() const = 0;
//  /** get X-position of particle with label indexPart */
//  virtual double getX(int indexPart) const = 0;
//  /** get Y-position of particle with label indexPart */
//  virtual double getY(int indexPart) const = 0;
//  /** get Z-position of particle with label indexPart */
//  virtual double getZ(int indexPart) const = 0;
//  /** get u (X-velocity) of particle with label indexPart */
//  virtual double getU(int indexPart) const = 0;
//  /** get v (Y-velocity) of particle with label indexPart */
//  virtual double getV(int indexPart) const = 0;
//  /** get w (Z-velocity) of particle with label indexPart */
//  virtual double getW(int indexPart) const = 0;
//  /** get ID of particle with label indexPart */
//  virtual long long getParticleID(int indexPart) const = 0;
//  /**get charge of particle with label indexPart */
//  virtual double getQ(int indexPart) const = 0;
//  /** get the number of particles of this subdomain */
//  virtual int getNOP() const = 0;
//  /** return the Kinetic energy */
//  virtual double getKe() = 0;
//  /** return the maximum kinetic energy */
//  virtual double getMaxVelocity() = 0;
//  /** return energy distribution*/
//  virtual long long *getVelocityDistribution(int nBins, double maxVel) = 0;
//  /** retturn the momentum */
//  virtual double getP() = 0;
//  /** Print particles info: positions, velocities */
//  virtual void Print(VirtualTopology3D * ptVCT) const = 0;
//  /** Print the number of particles of this subdomain */
//  virtual void PrintNp(VirtualTopology3D * ptVCT) const = 0;
//
//};
#endif
