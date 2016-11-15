Module ArGasModule
     Use TypeModule
     Use Constants
     Use SigmaModule
   !  This module defines the gas properties.
   Integer(4),Private,parameter :: NSpecyAr=1_4
   !Type(GasPhysical),parameter :: ArGas=(GasPhysical(2,1,1,6.63d-26,0.d0,0.d0)) 

    !  This module defines the particles properties.   
   !Type(OneSpecy),parameter :: ArSpecy(1:NSpecyAr)=(OneSpecy('Ar+',1,1,6.63d-26))
    
   !  This module defines the Collision properties. 
   Type(ReactionBundle),save :: ArReaction(0:NSpecyAr)

   Data ArReaction(0)%NReaction /3_4/
        Data ArReaction(0)%Reaction(1:3) /&
                      ReactionOne(0,11,b'01',0.d0),&      ! Elastic "e + Ar ->   e + Ar" 
                      ReactionOne(0,12,b'01',11.5d0),&  ! Inelastic "e + Ar -> e + Ar*" 11.5eV
                      ReactionOne(0,13,b'10',15.8d0)/     ! Ionisation "e + Ar -> 2e + Ar+"15.8eV
    Data ArReaction(1)%NReaction /2_4/
        Data ArReaction(1)%Reaction(1:2) /&
                      ReactionOne(1,101,b'10',0.d0),&    !elastic "Ar+ + Ar -> Ar+ + Ar"
                      ReactionOne(1,103,b'10',0.d0)/       !charge exchange "Ar+ + Ar -> Ar + Ar+"
          
   Type(SigmaNormalized),save :: ArSigma(0:NSpecyAr)
   
   contains
     Subroutine  ArSigmaInitialization()
        Implicit none
            !Call SigmaInit(0_4,ArSigma(0))
                  Call FuncSigmaNormalization(0_4,ArSigma(0),SigmaEArElastic)
                  Call FuncSigmaNormalization(0_4,ArSigma(0),SigmaEArInelastic)
                  Call FuncSigmaNormalization(0_4,ArSigma(0),SigmaEArIonisation)
            !Call SigmaInit(1_4,ArSigma(1))
                  Call FuncSigmaNormalization(1_4,ArSigma(1),SigmaArArElastic)
                  Call FuncSigmaNormalization(1_4,ArSigma(1),SigmaArArChargeExchange)
          Return        
     End Subroutine ArSigmaInitialization

    ! Crosssection (m^2) for elastic "e + Ar -> e + Ar" collision as Elastic(energy);
      Function SigmaEArElastic(energy)
        implicit none
        real(8) :: energy
        real(8) ::  SigmaEArElastic
            SigmaEArElastic=1.d-20*(dabs(6.0d0/(1.0d0+(energy/0.1d0)+(energy/0.6d0)**2)**3.3-  &
               1.1d0*energy**1.4d0/(1.0d0+(energy/1.50d1)**1.2)/dsqrt(1.0d0+(energy/5.5d0)**2.5+  &
               (energy/60.0d0)**4.1))+0.05d0/(1.0d0+(energy/1.0d1))**2.0+  &
               1d-2*energy**3.0/(1.0d0+(energy/1.2d1)**6.0))
        return
      end Function
    ! Crosssection (m^2) for inelastic "e + Ar -> e + Ar*" collision as Inelastic(energy);
      Function SigmaEArInelastic(energy)
        implicit none
        real(8) :: energy
        real(8) ::  SigmaEArInelastic
        if (energy<=11.5d0) then
                SigmaEArInelastic=0.d0
           else if (energy>11.5d0) then
                SigmaEArInelastic=1.d-20*(0.034d0*(energy-11.5d0)**1.1*(1.0d0+(energy/15.0d0)**2.8)/(1.0d0    &
                           +(energy/23.0d0)**5.5)+0.023d0*(energy-11.5d0)/(1.0d0+(energy/80.0d0)**1.9))
       end if
       return
      end Function
	! Crosssection (m^2) for ionisation "e + Ar -> 2e + Ar+" collision as Ionisation(energy);  
	  Function SigmaEArIonisation(energy)
    	implicit none
    	real(8) :: energy
        real(8) :: SigmaEArIonisation
        if (energy<=15.8d0) then
                SigmaEArIonisation=0.d0
           else if (energy>15.8d0) then
                SigmaEArIonisation=1.d-20*(970.d0*(energy-15.8d0)/(70.0+energy)**2.d0 &
			                 + 0.06d0*(energy-15.8d0)**2.d0*dexp(-energy/9.0d0))

       end if
	   return
	  end Function
	 ! Crosssection (m^2) for charge exchange "Ar+ + Ar -> Ar + Ar+" as ChargeExchange(energy); 
	  Function SigmaArArChargeExchange(energy)
    	implicit none
    	real(8),intent(in) :: energy
        real(8) ::  SigmaArArChargeExchange
        if (energy > 100.d0)  SigmaArArChargeExchange = 2.0d-19 + 5.5d-19/(Dsqrt(100.d0)+1.0d-30)
	    if (energy < 4.0d0) then 
	          SigmaArArChargeExchange = -2.95d-19*Dsqrt(energy)+10.65d-19
	    else
	         SigmaArArChargeExchange = 2.0d-19 + 5.5d-19/(Dsqrt(energy)+1.0d-30)
	    end if
	  return
	  end Function SigmaArArChargeExchange
	! Crosssection (m^2) for elastic "Ar+ + Ar -> Ar+ + Ar" IonElastic(double eps);
	  Function SigmaArArElastic(energy)
    	implicit none
    	real(8),intent(in) :: energy
        real(8) ::  SigmaArArElastic
        if (energy > 100.d0)  SigmaArArElastic = 1.8d-19 + 4.0d-19/(Dsqrt(100.d0)+1.0d-30)
	    if (energy < 4.d0) then
	         SigmaArArElastic = -2.0d-19*Dsqrt(energy)+7.8d-19
	    else
		     SigmaArArElastic = 1.8d-19 + 4.0d-19/(Dsqrt(energy)+1.0d-30)
	    end if
	  return 
	  end function SigmaArArElastic
end module ArGasModule
