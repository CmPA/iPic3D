#ifndef _Particle_
#define _Particle_

// Depends on width of vector unit;
// need to be known at compile time.
//
#define AoS_PCLS_AT_A_TIME 2

namespace ParticleType
{
  enum Type
  {
    AoS = 0,
    SoA
  };
}

// intended to occupy 64 bytes
//
// particle for a specific species
class SpeciesParticle
{
  long long ID;
  double x[3];
  double u[3];
  double q;
 public:
  // accessors
  long long get_ID()const{ return ID; }
  double get_x(int i)const{ return x[i]; }
  double get_u(int i)const{ return u[i]; }
  double get_q()const{ return q; }
  void set_ID(long long in){ ID=in; }
  void set_x(int i, double in) { x[i] = in; }
  void set_u(int i, double in) { u[i] = in; }
  void set_q(double in) { q = in; }
  // alternative accessors
  double get_x()const{ return x[0]; }
  double get_y()const{ return x[1]; }
  double get_z()const{ return x[2]; }
  double get_u()const{ return u[0]; }
  double get_v()const{ return u[1]; }
  double get_w()const{ return u[2]; }
  void set_x(double in){ x[0]=in; }
  void set_y(double in){ x[1]=in; }
  void set_z(double in){ x[2]=in; }
  void set_u(double in){ u[0]=in; }
  void set_v(double in){ u[1]=in; }
  void set_w(double in){ u[2]=in; }
  void set(long long _ID,
    double _x, double _y, double _z,
    double _u, double _v, double _w,
    double _q)
  {
    ID = _ID;
    x[0] = _x; x[1] = _y; x[2] = _z;
    u[0] = _u; u[1] = _v; u[2] = _w;
    q = _q;
  }
};

// intended to occupy 64 bytes
//
// to be used when sorting with every particle advance
struct CellParticle
{
  long long ID; // 8 bytes
  int cx[3]; // mesh cell
  float fx[3]; // mesh cell position (fraction)
  float u[3];
  float fxavg[3]; // for implicit push
  float q; // float m would be better for stitching to MHD for dusty plasma
  float qom; // for dusty plasma
 public:
  // accessors
  //
  // read access
  long long get_ID()const{ return ID; }
  float get_fx()const{ return fx[0]; }
  float get_fy()const{ return fx[1]; }
  float get_fz()const{ return fx[2]; }
  float get_u()const{ return u[0]; }
  float get_v()const{ return u[1]; }
  float get_w()const{ return u[2]; }
  float get_q()const{ return q; }
  void set_ID(long long in){ ID=in; }
  // write access
  void set_u(float in){ u[0]=in; }
  void set_v(float in){ u[1]=in; }
  void set_w(float in){ u[2]=in; }

  void init(const SpeciesParticle& pcl,
    double cxstart[3], // starting position of cell coordinates
    float dx_inv[3],
    float _qom)
  {
    ID = pcl.get_ID();
    // position in mesh coordinates
    //

    float xpos[3];
    for(int i=0;i<3;i++)
    {
      float xpos = (pcl.get_x(i)-cxstart[i])*dx_inv[i];
      float cxpos = floor(xpos);
      cx[i] = int(cxpos);
      fxavg[i] = fx[i] = cxpos - cx[i];
      u[i] = pcl.get_u(i);
    }
    q = pcl.get_q();
    qom = _qom;
  }
};

#endif
