Module ParticleModule
      ! This module is the operations for the particles.
      Use TypeModule
      Use Constants 
     Implicit none
     Contains
           Subroutine CalEnergy(Mass,InputParticle,Energy) 
              Implicit none
              Real(8),intent(in):: Mass 
              Type(ParticleOne),intent(in) :: InputParticle
              Real(8),intent(out) ::Energy
              Energy=InputParticle%Vx*InputParticle%Vx+InputParticle%Vy*InputParticle%Vy
              Energy=0.5d0*Mass*(Energy+InputParticle%Vz*InputParticle%Vz)
             return 
           end  Subroutine  CalEnergy
           
           Subroutine MometumEnergyChange(Mass,PN,PO)
                     Implicit none
                     Type(ParticleOne),intent(inout) :: PN
                     Type(ParticleOne),intent(in),optional :: PO
                     Real(8) :: Mass
                     Real(8) :: EnergyNew,EnergyOld
                     If (Present(PO)) then
                       Call CalEnergy(Mass,PN,EnergyNew)
                       Call CalEnergy(Mass,PO,EnergyOld)
                           
                       PN%dP_x=Mass*(PN%Vx-PO%Vx)
                       PN%dP_y=Mass*(PN%Vy-PO%Vy)
                       PN%dP_z=Mass*(PN%Vz-PO%Vz)
                       PN%dE=EnergyNew-EnergyOld
                     Else
                       Call CalEnergy(Mass,PN,EnergyNew)
                       PN%dP_x=Mass*PN%Vx
                       PN%dP_y=Mass*PN%Vy
                       PN%dP_z=Mass*PN%Vz
                       PN%dE=EnergyNew
                    ENd if
                       
                    return  
            End subroutine MometumEnergyChange
 End  Module ParticleModule

