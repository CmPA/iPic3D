/*************************************************************************** 
 EBox.h  -  Abstract Base class for expanding box        
 -------------------                   
begin                : Fri Jun 8 2018   
copyright            : (C) ???
developer            : Maria Elena Innocenti
email                : maria.elena.innocenti@jpl.nasa.gov
***************************************************************************/

#ifndef EBox_H
#define EBox_H

#include "Collective.h"

class EBox                                            
{
public:
  /** constructor */
  EBox(Collective * col);
  /** destructor */
  ~EBox();
  /** update expanding box parameters, 
      to be executed at the beginning of each cycle **/
  void UpdateEbParameter();

  /** output **/
  double getUEB_0();
  double getREB_0();
  double getR_nth();
  double getR();
private:
  double dt;
  double th;
  // velocity of the plasma parcel, from inputfile
  double UEB_0;
  // initial distance from the Sun, from inputfile
  double REB_0;

  /* R^{n+theta}= R0+U0 t = Rn+ U0*th dt */
  double R_nth;

  double R;
};
#endif
