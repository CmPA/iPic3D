
#include "BCStructure.h"

injInfoFields::injInfoFields(int Nxsize, int Nysize, int Nzsize) {
  ExITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  EyITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  EzITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  BxITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  ByITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  BzITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  Nxsize_store = Nxsize;
  Nysize_store = Nysize;
};
injInfoFields::~injInfoFields() {
  delArr3(ExITemp, Nxsize_store, Nysize_store);
  delArr3(EyITemp, Nxsize_store, Nysize_store);
  delArr3(EzITemp, Nxsize_store, Nysize_store);
  delArr3(BxITemp, Nxsize_store, Nysize_store);
  delArr3(ByITemp, Nxsize_store, Nysize_store);
  delArr3(BzITemp, Nxsize_store, Nysize_store);
};

injInfoParticles::injInfoParticles(int Nxsize, int Nysize, int Nzsize) {
  VxITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  VyITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  VzITemp = newArr3(double, Nxsize, Nysize, Nzsize);

  VthxITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  VthyITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  VthzITemp = newArr3(double, Nxsize, Nysize, Nzsize);

  Npcelx_array = newArr3(int, Nxsize, Nysize, Nzsize);
  Npcely_array = newArr3(int, Nxsize, Nysize, Nzsize);
  Npcelz_array = newArr3(int, Nxsize, Nysize, Nzsize);

  QITemp = newArr3(double, Nxsize, Nysize, Nzsize);
  RoITemp = newArr3(double, Nxsize, Nysize, Nzsize);

  Nxsize_store = Nxsize;
  Nysize_store = Nysize;
};

injInfoParticles::~injInfoParticles() {
  delArr3(VxITemp, Nxsize_store, Nysize_store);
  delArr3(VyITemp, Nxsize_store, Nysize_store);

  delArr3(VzITemp, Nxsize_store, Nysize_store);

  delArr3(VthxITemp, Nxsize_store, Nysize_store);
  delArr3(VthyITemp, Nxsize_store, Nysize_store);
  delArr3(VthzITemp, Nxsize_store, Nysize_store);

  delArr3(Npcelx_array, Nxsize_store, Nysize_store);
  delArr3(Npcely_array, Nxsize_store, Nysize_store);
  delArr3(Npcelz_array, Nxsize_store, Nysize_store);

  delArr3(QITemp, Nxsize_store, Nysize_store);
  delArr3(RoITemp, Nxsize_store, Nysize_store);
};
