#ifndef injIFields
#define injIFields

#include "Alloc.h"

struct injInfoFields {
  double ***ExITemp;
  double ***EyITemp;
  double ***EzITemp;
  double ***BxITemp;
  double ***ByITemp;
  double ***BzITemp;
  int Nxsize_store;
  int Nysize_store;

    injInfoFields(int Nxsize, int Nysize, int Nzsize);
   ~injInfoFields();
};


struct injInfoParticles {
  double ***VxITemp;
  double ***VyITemp;
  double ***VzITemp;

  double ***VthxITemp;
  double ***VthyITemp;
  double ***VthzITemp;

  int ***Npcelx_array;
  int ***Npcely_array;
  int ***Npcelz_array;

  double ***QITemp;             // weight of injected particles: make it variable, keeping the amount of injected particles.

  double ***RoITemp;
  int Nxsize_store;
  int Nysize_store;

    injInfoParticles(int Nxsize, int Nysize, int Nzsize);
   ~injInfoParticles();
};

#endif
