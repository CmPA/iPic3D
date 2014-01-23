#ifndef _Parameters_h_
#define _Parameters_h_

// namespace provides a more flexible, succinct singleton via "using Parameters"
//
namespace Parameters
{
  enum MoverType
  {
    SoA=0,
    AoS,
    SoAvec_onesort,
    AoSvec_onesort,
    SoAvec_resort,
    AoSvec_resort,
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
  MoverType get_MOVER_TYPE();
}
#endif
