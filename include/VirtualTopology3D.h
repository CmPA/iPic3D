/***************************************************************************
  VirtualTopology.h - Abstract Base class for virtual process topologies
  -------------------
begin                : May 2008
copyright            : KUL Leuven
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef VirtualTopology3D_H
#define VirtualTopology3D_H

#include "mpi.h"
/* mlmd: need to include Collective.h for new setup_vctopology */
#include "Collective.h"

/**
 *  Abstract base class for virtual process topologies
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
class VirtualTopology3D {
public:
  /** Find the neighbors in the new communicator  */
  /* pre-mlmd
     virtual void setup_vctopology(MPI_Comm comm_old) = 0; */
  virtual void setup_vctopology(MPI_Comm comm_old, Collective *col) = 0;
  /** Print topology info */
  virtual void Print() = 0;
  /** Print the mapping of topology */
  virtual void PrintMapping() = 0;
  /** these give info on the current grid **/
  /** get XLEN */
  virtual int getXLEN() = 0;
  /** get YLEN */
  virtual int getYLEN() = 0;
  /** get ZLEN */
  virtual int getZLEN() = 0;
  /** these give info on thr grid system **/
  /** get XLEN */
  virtual int getXLEN(int N) = 0;
  /** get YLEN */
  virtual int getYLEN(int N) = 0;
  /** get ZLEN */
  virtual int getZLEN(int N) = 0;
  /** get nprocs */
  virtual int getNprocs() = 0;
  /** get periodicity on boundaries - DIRECTION X*/
  virtual bool getPERIODICX() = 0;
  virtual bool getPERIODICX_P() = 0;
  /** get periodicity on boundaries - DIRECTION Y*/
  virtual bool getPERIODICY() = 0;
  virtual bool getPERIODICY_P() = 0;
  /** get periodicity on boundaries - DIRECTION Z*/
  virtual bool getPERIODICZ() = 0;
  virtual bool getPERIODICZ_P() = 0;
  /** get the cartesian rank of the process */
  virtual int getCartesian_rank() = 0;
  /** get the total number of process */
  virtual int getNproc() = 0;
  /** get the cartesian rank of XLEFT neighbor */
  virtual int getXleft_neighbor() = 0;
  /** get the cartesian rank of XRIGHT neighbor */
  virtual int getXright_neighbor() = 0;
  /** get the cartesian rank of YLEFT neighbor */
  virtual int getYleft_neighbor() = 0;
  /** get the cartesian rank of YRIGHT neighbor */
  virtual int getYright_neighbor() = 0;
  /** get the cartesian rank of ZLEFT neighbor */
  virtual int getZleft_neighbor() = 0;
  /** get the cartesian rank of ZRIGHT neighbor */
  virtual int getZright_neighbor() = 0;
  /** get the PARTICLES cartesian rank of XLEFT neighbor */
  virtual int getXleft_neighbor_P() = 0;
  /** get the PARTICLES cartesian rank of XRIGHT neighbor */
  virtual int getXright_neighbor_P() = 0;
  /** get the PARTICLES cartesian rank of YLEFT neighbor */
  virtual int getYleft_neighbor_P() = 0;
  /** get the PARTICLES cartesian rank of YRIGHT neighbor */
  virtual int getYright_neighbor_P() = 0;
  /** get the PARTICLES cartesian rank of ZLEFT neighbor */
  virtual int getZleft_neighbor_P() = 0;
  /** get the PARTICLES cartesian rank of ZRIGHT neighbor */
  virtual int getZright_neighbor_P() = 0;
  /** get the coordinates in dir direction of process*/
  virtual int getCoordinates(int dir) = 0;
  /** get the coordinates of process*/
  virtual int *getCoordinates() = 0;
  /** get the number of processors in each direction*/
  virtual int *getDivisions() = 0;
  /** get Periodicity condition in dir direction */
  virtual int getPeriods(int dir) = 0;
  /** if cVERBOSE == true, print to the screen all the comunication */
  virtual bool getcVERBOSE() = 0;
  /** get the MPI communicator */
  virtual MPI_Comm getComm() = 0;

  /*! mlmd specific functions */
  /*! returns the non cartesian communicator at grid level */
  virtual MPI_Comm getCommGrid() = 0;
  /*! returns the number of the current grid */
  virtual int getNumGrid() = 0;
  /*! returns the communicator to the parent grid */
  virtual MPI_Comm getCommToParent() = 0;
  virtual MPI_Comm getCommToParent_BCGhost() = 0;
  virtual MPI_Comm getCommToParent_BCBuffer() = 0;
  virtual MPI_Comm getCommToParent_BCFix3B() = 0;
  virtual MPI_Comm getCommToParent_Proj() = 0;
  /*! returns the communicator to the parent grid for particles */
  virtual MPI_Comm getCommToParent_P(int is) = 0;
  /*! returns the number of children of a grid */ 
  virtual int getNumChildren() =0;
  /* return the communicator to the child grid n in the mlmd hierarchy */
  // for fields
  virtual MPI_Comm getCommToChild(int n) =0;
  virtual MPI_Comm getCommToChild_BCGhost(int n) =0;
  virtual MPI_Comm getCommToChild_BCBuffer(int n) =0;
  virtual MPI_Comm getCommToChild_BCFix3B(int n) =0;
  virtual MPI_Comm getCommToChild_Proj(int n) =0;
  // for particles
  virtual MPI_Comm getCommToChild_P(int n, int is) =0;

  virtual int getSystemWide_rank()=0;
  virtual int getRank_CommToParent()=0;
  virtual int getRank_CommToParent_P(int is)=0;
  virtual int getRank_CommToChildren(int nc)=0;

  virtual int getParentGridNum()=0;
  virtual int getChildGridNum(int n) =0;
  
  virtual int getMaxGridCoreN()=0;
  virtual int getMaxGridPer()=0;
  virtual int getMaxRF1()=0;

  virtual int getRank_CommToChildren_P(int nc, int is) =0;
  /* return the values of the cartesian coordinate lookup table  
     NB: N is the rank number in the inter-communicator (CommToParent or CommToChildren)*/
  virtual int getXcoord_CommToParent(int N)=0;
  virtual int getYcoord_CommToParent(int N)=0;
  virtual int getZcoord_CommToParent(int N)=0;
  /*! end mlmd specific functions */

  virtual MPI_Comm getCommField_XLeft() =0; 
  virtual MPI_Comm getCommField_XRight() =0; 
  virtual MPI_Comm getCommField_YLeft() =0; 
  virtual MPI_Comm getCommField_YRight() =0; 
  virtual MPI_Comm getCommField_ZLeft() =0; 
  virtual MPI_Comm getCommField_Zright() =0; 

};
#endif
