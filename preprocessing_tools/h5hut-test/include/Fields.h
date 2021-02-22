#ifndef __FIELDS_H__
#define __FIELDS_H__

#include <cmath>
#include "Alloc.h"

class Fields{
  public:
    Fields();
    ~Fields();
    void InitRandom(int, int, int, int);
    void InitCell(int, int, int, int);
    void InitNode(int, int, int, int);
    void InitSequentialCell(int, int, int, int, int);
    void InitSequentialNode(int, int, int, int, int);
    void Init(int, int, int, int);
    void FreeCell();
    void FreeNode();

    void   SetEn  (int i, int j, int k, double E0x, double E0y, double E0z) {Ex[i][j][k] = E0x; Ey[i][j][k] = E0y; Ez[i][j][k] = E0z;};
    void   SetBn  (int i, int j, int k, double B0x, double B0y, double B0z) {Bx[i][j][k] = B0x; By[i][j][k] = B0y; Bz[i][j][k] = B0z;};
    void   SetRhon(int s, int i, int j, int k, double rho0) {rhons[s][i][j][k] = rho0;};

    void   SetEc  (int i, int j, int k, double E0x, double E0y, double E0z) {Exc[i][j][k] = E0x; Eyc[i][j][k] = E0y; Ezc[i][j][k] = E0z;};
    void   SetBc  (int i, int j, int k, double B0x, double B0y, double B0z) {Bxc[i][j][k] = B0x; Byc[i][j][k] = B0y; Bzc[i][j][k] = B0z;};
    void   SetRhoc(int s, int i, int j, int k, double rho0) {rhocs[s][i][j][k] = rho0;};

    double ShowEx (int i, int j, int k) {return Ex [i][j][k];}
    double ShowExc(int i, int j, int k) {return Exc[i][j][k];}

    double ***getExc() { return(Exc);}
    double ***getEyc() { return(Eyc);}
    double ***getEzc() { return(Ezc);}
    double ***getBxc() { return(Bxc);}
    double ***getByc() { return(Byc);}
    double ***getBzc() { return(Bzc);}
    double ***& getRHOcs(int s) { return(rhocs[s]);}

    double ***getEx() { return(Ex);}
    double ***getEy() { return(Ey);}
    double ***getEz() { return(Ez);}
    double ***getBx() { return(Bx);}
    double ***getBy() { return(By);}
    double ***getBz() { return(Bz);}
    double ***& getRHOns(int s) { return(rhons[s]);}

  private:
    double ***Ex;
    double ***Ey;
    double ***Ez;
    double ***Bx;
    double ***By;
    double ***Bz;
    double ****rhons;

    double ***Exc;
    double ***Eyc;
    double ***Ezc;
    double ***Bxc;
    double ***Byc;
    double ***Bzc;
    double ****rhocs;
};

Fields::Fields(){
}

Fields::~Fields(){
}

void Fields::FreeNode(){
  int d1, d2, d3;

  delArr3<double>(Ex, d1, d2);
  delArr3<double>(Ey, d1, d2);
  delArr3<double>(Ez, d1, d2);
  delArr3<double>(Bx, d1, d2);
  delArr3<double>(By, d1, d2);
  delArr3<double>(Bz, d1, d2);
  delArr4<double>(rhons, d1, d2, d3);
}

void Fields::FreeCell(){
  int d1, d2, d3;

  delArr3<double>(Exc, d1, d2);
  delArr3<double>(Eyc, d1, d2);
  delArr3<double>(Ezc, d1, d2);
  delArr3<double>(Bxc, d1, d2);
  delArr3<double>(Byc, d1, d2);
  delArr3<double>(Bzc, d1, d2);
  delArr4<double>(rhocs, d1, d2, d3);
}

void Fields::Init(int nspec, int nxc, int nyc, int nzc){
  Exc   = newArr3(double, nxc, nyc, nzc);
  Eyc   = newArr3(double, nxc, nyc, nzc);
  Ezc   = newArr3(double, nxc, nyc, nzc);
  Bxc   = newArr3(double, nxc, nyc, nzc);
  Byc   = newArr3(double, nxc, nyc, nzc);
  Bzc   = newArr3(double, nxc, nyc, nzc);
  rhocs = newArr4(double, nspec, nxc, nyc, nzc);

  for (int i=0; i<nxc; i++){
    for (int j=0; j<nyc; j++){
      for (int k=0; k<nzc; k++){
        Exc[i][j][k] = 0.0;
        Eyc[i][j][k] = 0.0;
        Ezc[i][j][k] = 0.0;
        Bxc[i][j][k] = 0.0;
        Byc[i][j][k] = 0.0;
        Bzc[i][j][k] = 0.0;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nxc; i++){
      for (int j=0; j<nyc; j++){
        for (int k=0; k<nzc; k++){
          rhocs[s][i][j][k] = 0.0;
        }
      }
    }
  }

  Ex    = newArr3(double, nxc+1, nyc+1, nzc+1);
  Ey    = newArr3(double, nxc+1, nyc+1, nzc+1);
  Ez    = newArr3(double, nxc+1, nyc+1, nzc+1);
  Bx    = newArr3(double, nxc+1, nyc+1, nzc+1);
  By    = newArr3(double, nxc+1, nyc+1, nzc+1);
  Bz    = newArr3(double, nxc+1, nyc+1, nzc+1);
  rhons = newArr4(double, nspec, nxc+1, nyc+1, nzc+1);

  for (int i=0; i<nxc+1; i++){
    for (int j=0; j<nyc+1; j++){
      for (int k=0; k<nzc+1; k++){
        Ex[i][j][k] = 0.0;
        Ey[i][j][k] = 0.0;
        Ez[i][j][k] = 0.0;
        Bx[i][j][k] = 0.0;
        By[i][j][k] = 0.0;
        Bz[i][j][k] = 0.0;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nxc+1; i++){
      for (int j=0; j<nyc+1; j++){
        for (int k=0; k<nzc+1; k++){
          rhons[s][i][j][k] = 0.0;
        }
      }
    }
  }
}

void Fields::InitSequentialNode(int nspec, int nx, int ny, int nz, int tag){
  Ex    = newArr3(double, nx, ny, nz);
  Ey    = newArr3(double, nx, ny, nz);
  Ez    = newArr3(double, nx, ny, nz);
  Bx    = newArr3(double, nx, ny, nz);
  By    = newArr3(double, nx, ny, nz);
  Bz    = newArr3(double, nx, ny, nz);
  rhons = newArr4(double, nspec, nx, ny, nz);

  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        Ex[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        Ey[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        Ez[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        Bx[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        By[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        Bz[i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          rhons[s][i][j][k] = tag*1000 + k*ny*nx + j*nx + i;
        }
      }
    }
  }
}

void Fields::InitNode(int nspec, int nx, int ny, int nz){
  Ex    = newArr3(double, nx, ny, nz);
  Ey    = newArr3(double, nx, ny, nz);
  Ez    = newArr3(double, nx, ny, nz);
  Bx    = newArr3(double, nx, ny, nz);
  By    = newArr3(double, nx, ny, nz);
  Bz    = newArr3(double, nx, ny, nz);
  rhons = newArr4(double, nspec, nx, ny, nz);

  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        Ex[i][j][k] = 1.0;
        Ey[i][j][k] = 2.0;
        Ez[i][j][k] = 3.0;
        Bx[i][j][k] = 4.0;
        By[i][j][k] = 5.0;
        Bz[i][j][k] = 6.0;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          rhons[s][i][j][k] = 7.0+s;
        }
      }
    }
  }
}

void Fields::InitSequentialCell(int nspec, int nxc, int nyc, int nzc, int tag){
  Exc   = newArr3(double, nxc, nyc, nzc);
  Eyc   = newArr3(double, nxc, nyc, nzc);
  Ezc   = newArr3(double, nxc, nyc, nzc);
  Bxc   = newArr3(double, nxc, nyc, nzc);
  Byc   = newArr3(double, nxc, nyc, nzc);
  Bzc   = newArr3(double, nxc, nyc, nzc);
  rhocs = newArr4(double, nspec, nxc, nyc, nzc);

  for (int i=0; i<nxc; i++){
    for (int j=0; j<nyc; j++){
      for (int k=0; k<nzc; k++){
        Exc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        Eyc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        Ezc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        Bxc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        Byc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        Bzc[i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nxc; i++){
      for (int j=0; j<nyc; j++){
        for (int k=0; k<nzc; k++){
          rhocs[s][i][j][k] = tag*1000 + k*nyc*nxc + j*nxc + i;
        }
      }
    }
  }
}

void Fields::InitCell(int nspec, int nxc, int nyc, int nzc){
  Exc   = newArr3(double, nxc, nyc, nzc);
  Eyc   = newArr3(double, nxc, nyc, nzc);
  Ezc   = newArr3(double, nxc, nyc, nzc);
  Bxc   = newArr3(double, nxc, nyc, nzc);
  Byc   = newArr3(double, nxc, nyc, nzc);
  Bzc   = newArr3(double, nxc, nyc, nzc);
  rhocs = newArr4(double, nspec, nxc, nyc, nzc);

  for (int i=0; i<nxc; i++){
    for (int j=0; j<nyc; j++){
      for (int k=0; k<nzc; k++){
        Exc[i][j][k] = 0.0;
        Eyc[i][j][k] = 0.0;
        Ezc[i][j][k] = 0.0;
        Bxc[i][j][k] = 0.0;
        Byc[i][j][k] = 0.0;
        Bzc[i][j][k] = 0.0;
      }
    }
  }

  for (int s=0; s<nspec; s++){
    for (int i=0; i<nxc; i++){
      for (int j=0; j<nyc; j++){
        for (int k=0; k<nzc; k++){
          rhocs[s][i][j][k] = 0.0;
        }
      }
    }
  }
}

void Fields::InitRandom(int nspec, int nx, int ny, int nz){
  Ex    = newArr3(double, nx, ny, nz);
  Ey    = newArr3(double, nx, ny, nz);
  Ez    = newArr3(double, nx, ny, nz);
  Bx    = newArr3(double, nx, ny, nz);
  By    = newArr3(double, nx, ny, nz);
  Bz    = newArr3(double, nx, ny, nz);
  rhons = newArr4(double, nspec, nx, ny, nz);

  double LO = -10.0;
  double HI =   0.0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        if (i>10) {
          LO = 0.0;
          HI = 10.0;
        }
        Ex[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        Ey[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        Ez[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        Bx[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        By[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        Bz[i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
      }
    }
  }

  LO = 0.0;
  HI = 1.0;
  for (int s=0; s<nspec; s++){
    for (int i=0; i<nx; i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          rhons[s][i][j][k] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
        }
      }
    }
  }
}

#endif
