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
// species particle for second-order-accuracy implicit advance
class ISpcl
{
  long long ID;
  double x[3];
  float hdx[3]; // xavg = x + hdx
  float u[3];
  double q;
 public:
  // accessors
  long long get_ID()const{ return ID; }
  double get_x(int i)const{ return x[i]; }
  double get_hdx(int i)const{ return hdx[i]; }
  float get_u(int i)const{ return u[i]; }
  double get_q()const{ return q; }
  void set_ID(long long in){ ID=in; }
  void set_x(int i, double in) { x[i] = in; }
  void set_hdx(int i, float in) { hdx[i] = in; }
  void set_u(int i, double in) { u[i] = in; }
  void set_q(double in) { q = in; }
};

// intended to occupy 64 bytes
//
// dust particle for second-order-accuracy implicit advance
class IDpcl
{
  long long ID;
  double x[3]; // could replace with cell index and float x
  float hdx[3]; // xavg = x + hdx
  float u[3];
  float qom; // charge to mass ratio of particle
  float m; // mass of particle
};

#endif
