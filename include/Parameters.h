#ifndef _Parameters_h_
#define _Parameters_h_

// namespace provides a more flexible, succinct singleton via "using Parameters"
//
namespace Parameters
{
  void init_parameters();

  bool get_SORTING_PARTICLES();
  bool get_VECTORIZE_MOMENTS();
  bool get_VECTORIZE_MOVER();
  bool get_USING_XAVG();
}
#endif
