#include "Parameters.h"

using namespace Parameters;

static bool SORTING_PARTICLES;

void Parameters::init_parameters()
{
  SORTING_PARTICLES = get_VECTORIZE_MOMENTS() || get_VECTORIZE_MOVER();
}

//bool Parameters::get_SORTING_PARTICLES() { return SORTING_PARTICLES; }
bool Parameters::get_SORTING_PARTICLES() { return true; }
bool Parameters::get_VECTORIZE_MOMENTS() { return false; }
bool Parameters::get_VECTORIZE_MOVER() { return false; }
// this must also return true if we communicate particles per iteration
//bool Parameters::get_USING_XAVG() { return get_VECTORIZE_MOVER(); }
bool Parameters::get_USING_XAVG() { return get_SORTING_PARTICLES(); }
