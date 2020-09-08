#include "EBox.h"
#include <iomanip>
using std::showpoint;
/*! constructor */
EBox::EBox(Collective * col) {

  dt    = col->getDt();
  th    = col->getTh();
  UEB_0 = col->getUEB_0();
  REB_0 = col->getREB_0();

  // this variable is to be used when restarting an EB simulation starting
  // from a non EB one
  // if it is false, R= R_0 + U_0*dt*cycles
  // if it is true, R= R_0
  bool Restart_From_REB_0= col->getRestart_From_REB_0();
  
  int restart_or_solinit= col->getrestart_or_solinit();
  if (restart_or_solinit ==0 or (restart_or_solinit ==1 and Restart_From_REB_0 )) {
    //so it becomes REB_0 + dt*UEB_0*th at the first update
    R_nth  = REB_0- dt*UEB_0 + dt*UEB_0*th; 
    //so it becomes REB_0 at tge first update
    R  = REB_0- dt*UEB_0;
  }
  else // restart, either for hdf5 or h5hut
    {
      R= REB_0 + col->getLast_cycle()*dt*UEB_0 ; // check: you do not have to remove dt*UEB_0
      R_nth= R+    dt*UEB_0*th;
    }
  
}

/*! destructor */
EBox::~EBox(){

}

/*! update expanding both parameters: to be done once at the beginning of the cycle */
void EBox::UpdateEbParameter(){
  

  /* EVERY TIME STEP, IT HAS TO BE UPDATED OF A DT!!!       
     (AT INIT, SET TO R_nth  = REB_0- dt*UEB_0 + dt*UEB_0*th, so it becomes REB_0 + dt*UEB_0*th
     at the first update*/
  R_nth = R_nth+ dt*UEB_0;

  R = R+ dt*UEB_0;

  cout.precision(11);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) {

    cout << "Grid intermediate position updated to Rnth=" << R_nth << " in EBox" << endl;
    cout << "REB_0: " << REB_0 << " UEB_0: " << UEB_0 << " th: " << th << " dt: " << dt << endl;
   
  }

}


/** output **/
double EBox::getUEB_0() {return UEB_0;}
double EBox::getREB_0() {return REB_0;}
double EBox::getR_nth() {return R_nth;}
double EBox::getR() {return R;}
