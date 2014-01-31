#ifndef _Parameters_h_
#define _Parameters_h_

// namespace provides a more flexible, succinct singleton via "using Parameters"
//
namespace Parameters
{
  enum Enum
  {
    SoA=0, // struct of arrays
    AoS, // array of structs
    // for moments type
    AoSvec,
    SoAvec,
    // for mover type
    SoA_vec_onesort,
    AoS_vec_onesort,
    SoA_vec_resort,
    AoS_vec_resort,
  };

  void init_parameters();

  bool get_USING_AOS();
  bool get_SORTING_SOA();
  bool get_SORTING_PARTICLES();
  // for resorting particles with each iteration of mover
  bool get_RESORTING_PARTICLES();
  inline bool get_USING_XAVG() { return get_RESORTING_PARTICLES(); }
  bool get_VECTORIZE_MOMENTS();
  //bool get_VECTORIZE_MOVER();
  Enum get_MOVER_TYPE();
  Enum get_MOMENTS_TYPE();
}
#endif
