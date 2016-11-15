!Module ElectronModule
! Use TypeModule
!    Implicit none
!        Type(OneSpecy),parameter :: Electron=(OneSpecy('Electron',1,-1,9.1095d-31))
!        Type(Gas),parameter :: ElectronGas=(Gas(1,0,9.1095d-31,3.d0*11605.d0,0.d0,0.d0)) 
! End Module ElectronModule

Module ElectronMCCModule
    Use MCCModule
    Use ParticleModule 
    Implicit none
    contains 
   !  ReactionType is definition of the reactions:
   !   -1-The injected particle are removed because the generated particle is no longer included in simulations . 
   !   0-no collision
   !   1~99 for electrons  (1~10 Isotropic : 1-Elastic; 2-Exitation; 3- Ionization; 4-Attachment; 5-Dissociation)
   !                                    (11-20 Anisotropic Ar: 11-Elastic; 12-Exitation; 13- Ionization;)
    subroutine SelectCollisionElectron(InputParticle,InputGas,InputReaction,NPar,OuputParticle)          
         Implicit none
         Type(MCCParticleOne),intent(inout) :: InputParticle
         Type(ReactionOne),intent(in) :: InputReaction
         Type(Gas),intent(in) :: InputGas
         Integer(4),intent(inout)  :: NPar
         Type(ParticleOne) :: OuputParticle(3)
         Integer(4) :: NSpecy,ParticleIndex(2_4)
         !Type(ParticleBundleSmall),intent(Out) :: OuputParticle(0:InputGas%NSpecy)
         Type(ParticleOne) :: TempParticle(2_4) 
         NSpecy=InputGas%NSpecy  
             Select case(InputReaction%ReactionType)
             Case(0_4) ! No collisions occur.
                          !Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                 NPar=1
                 OuputParticle(1)=InputParticle%PhaseSpace
                 OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                 OuputParticle(1)%CollisionType=-1
                  !Case(1_4)  ! Isotropic Elastic collisions.
                  !        Call ElectronPostVelocity(InputParticle,IsotropicCosKai)
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  !Case(2_4)  ! Isotropic Excitation collisions occur.
                  !        Call  ExcitationElectron(InputParticle,InputReaction%Threshold,IsotropicCosKai)
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  !Case(3_4)  ! Isotropic Ionization collisions.
                  !        Call IonizationElectron(InputParticle,InputGas,TempParticle,InputReaction%Threshold,IsotropicCosKai,IsotropicEnergy)
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  !        Call AddParticleSmall(TempParticle(1),OuputParticle(InputReaction%Reactant))
                  !        Call MCCParticleIndex(NSpecy,InputReaction%Resultant,1_4,ParticleIndex(1))
                  !         TempParticle(2)%Ax=-TempParticle(2)%Ax*InputParticle%Mass/InputGas%MGas
                  !          TempParticle(2)%Ay=-TempParticle(2)%Ay*InputParticle%Mass/InputGas%MGas
                  !           TempParticle(2)%Az=-TempParticle(2)%Az*InputParticle%Mass/InputGas%MGas
                  !        Call AddParticleSmall(TempParticle(2),OuputParticle(ParticleIndex(1)))
                  !Case(4_4)   ! Attachment collisions.
                  !        TempParticle(1)=InputParticle%PhaseSpace
                  !        Call Maxwellian(InputGas,TempParticle(1))
                  !        Call MCCParticleIndex(NSpecy,InputReaction%Resultant,1_4,ParticleIndex(1))
                  !        Call AddParticleSmall(TempParticle(1),OuputParticle(ParticleIndex(1)))  
                  !Case(5_4)
                  !        Call DissociationElectron(InputParticle,InputGas,TempParticle,InputReaction%Threshold,IsotropicCosKai,IsotropicEnergy)
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  !        Call MCCParticleIndex(NSpecy,InputReaction%Resultant,2_4,ParticleIndex(1:2))
                  !        Call AddParticleSmall(TempParticle(1),OuputParticle(ParticleIndex(1)))
                  !        Call AddParticleSmall(TempParticle(2),OuputParticle(ParticleIndex(2)))  
                   Case(11_4) !Anisotropic Elastic collisions for Ar Gas
                           Call ElectronPostVelocity(InputParticle,EArCosKai) 
                           NPar=1
                           OuputParticle(1)=InputParticle%PhaseSpace
                           OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                           OuputParticle(1)%CollisionType=0
                   Case(12_4) !!Anisotropic Excitation collisions for Ar Gas
                            Call  ExcitationElectron(InputParticle,InputReaction%Threshold,EArCosKai)
                            NPar=1
                            OuputParticle(1)=InputParticle%PhaseSpace
                            OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                            OuputParticle(1)%CollisionType=1
                   Case(13_4) !!Anisotropic Ionization collisions for Ar Gas
                            Call  IonizationElectron(InputParticle,InputGas,TempParticle,InputReaction%Threshold,EArCosKai) 
                            NPar=3
                            OuputParticle(1)=InputParticle%PhaseSpace
                            OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                            OuputParticle(1)%CollisionType=2
                            
                            TempParticle(1)%idp=666
                            OuputParticle(2)=TempParticle(1)
                            !OuputParticle(2)%sp_in=0
                            OuputParticle(2)%sp_out=OuputParticle(1)%sp_in
                            OuputParticle(2)%CollisionType=2
                            Call MometumEnergyChange(InputParticle%Mass,OuputParticle(2))
                            
                            Call MCCParticleIndex(NSpecy,InputReaction%Resultant,1_4,ParticleIndex(1))
                            TempParticle(2)%Sp_out=ParticleIndex(1)
                            TempParticle(2)%qp=-1.d0*TempParticle(2)%qp
                            TempParticle(2)%idp=666
                            !Write(*,*)  TempParticle(2)%A1,TempParticle(2)%A2,TempParticle(1)%A1,TempParticle(1)%A2
                            OuputParticle(3)=TempParticle(2)
                            !OuputParticle(2)%sp_in=1
                            OuputParticle(2)%CollisionType=2
                            Call MometumEnergyChange(InputParticle%Mass,OuputParticle(3))
                     End select
                            Call MometumEnergyChange(InputParticle%Mass,OuputParticle(1),InputParticle%Particle_Old)
        return
      contains
                subroutine ExcitationElectron(InputParticle,Threshold,CosTheta)
                    Implicit none
                    Type(MCCParticleOne),intent(inout) :: InputParticle
                     Real(8),intent(in) :: Threshold
                     Real(8),external :: CosTheta
                    Real(8) :: VFactor
                    VFactor=Dsqrt(1.d0-Threshold/InputParticle%Energy)
                    Call UpdateVelocityMCC(VFactor,InputParticle)
                    Call ElectronPostVelocity(InputParticle,CosTheta)
                    return
                end subroutine ExcitationElectron
               
                subroutine IonizationElectron(InputParticle,InputGas,OutputParticle,Threshold,CosTheta)
               ! This is the  subroutine for mainly the ionization collision of electrons. After the  collision, a new electron  and a new ion
                ! will be created. 
                    Implicit none
                    Type(MCCParticleOne),intent(inout) :: InputParticle
                    Type(Gas),intent(in) :: InputGas
                    Type(MCCParticleOne) :: ParticleTemp
                    Type(ParticleOne) ::  OutputParticle(2)
                    Real(8),intent(in) :: Threshold
                    Real(8),external :: CosTheta,CreateEnergy
                    Real(8) :: VFactor,EnergyTemp,EnergyOld,EnergyNew
                    ! Particle data are temperally stored.
                    ParticleTemp=InputParticle
                    ! Injected electron's velocity is determinated below.
                    
                    EnergyTemp=InputParticle%Energy-Threshold
                    Call DRandom(R)
                    EnergyOld=R*EnergyTemp
                    VFactor=Dsqrt(EnergyOld/InputParticle%Energy)
                    Call UpdateVelocityMCC(VFactor,InputParticle)
                    Call ElectronPostVelocity(InputParticle,CosTheta)
                    
                    ! Created electron's velocity is determinated below. Energy is conserved.
                     EnergyNew=EnergyTemp-EnergyOld
                     VFactor=Dsqrt(EnergyNew/ParticleTemp%Energy)
                     Call UpdateVelocityMCC(VFactor,ParticleTemp)
                     Call ElectronPostVelocity(ParticleTemp,CosTheta)
                     OutputParticle(1)=ParticleTemp%PhaseSpace 
                     OutputParticle(2)=InputParticle%PhaseSpace
                    ! Created ion's velocity is determinated below.
                     Call Maxwellian(InputGas,OutputParticle(2))
                    return
                End subroutine IonizationElectron

                Subroutine DissociationElectron(InputParticle,InputGas,OutputParticle,Threshold,CosTheta,CreateEnergy)
                    Implicit none
                    Type(MCCParticleOne),intent(inout) :: InputParticle
                    Type(Gas),intent(in) :: InputGas
                    Type(MCCParticleOne) :: ParticleTemp
                    Type(ParticleOne) ::  OutputParticle(2)
                    Real(8),intent(in) :: Threshold
                    Real(8),external :: CosTheta,CreateEnergy 
                    Real(8) :: VFactor,EnergyTemp
                    EnergyTemp=InputParticle%Energy-Threshold 
                    VFactor=Dsqrt(EnergyTemp/InputParticle%Energy)
                    Call UpdateVelocityMCC(VFactor,InputParticle)
                    Call ElectronPostVelocity(InputParticle,CosTheta)
                    OutputParticle(1)=InputParticle%PhaseSpace
                    OutputParticle(2)=InputParticle%PhaseSpace  
                    Call Maxwellian(InputGas,OutputParticle(1)) 
                    Call Maxwellian(InputGas,OutputParticle(2)) 
                    return
                End subroutine DissociationElectron

                Subroutine ElectronPostVelocity(InputParticle,CosKai)
                  Use Constants
                ! Note:Here the velocity of the target before the collision has be set to zero to save time.
                   implicit none
                   Type(MCCParticleOne),intent(inout) :: InputParticle
                   Real(8),external ::  CosKai
                   real(8) :: V
                   real(8) :: Vx,Vy,Vz,Vper,hx,hy,hz
                   real(8) :: AFai,CosFai,SinFai
                   real(8) ::  CosTheta,FcosTheta,SinTheta
                   Vx=InputParticle%PhaseSpace%Vx
                   Vy=InputParticle%PhaseSpace%Vy
                   Vz=InputParticle%PhaseSpace%Vz
                   V=InputParticle%V
                   !MassRatio=1.d0
                   
                   CosTheta=CosKai(InputParticle%Energy)
                   FcosTheta=1.d0-CosTheta
                   SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)

                   Call DRandom(R)
                   AFai=2.d0*Pi*R
                   CosFai=DCos(AFai)
                   SinFai=DSin(AFai)

                   Vper=Dsqrt(Vy*Vy+Vz*Vz)
                   If(Abs(Vper)<MinReal)  then
                           Call RandomVelocity(V,Vx,Vy,Vz)
                   else
                           hx=Vper*CosFai
                           hy=-(Vx*Vy*CosFai+V*Vz*SinFai)/Vper
                           hz=-(Vx*Vz*CosFai-V*Vy*SinFai)/Vper
                           Vx=Vx-(Vx*FcosTheta+hx*SinTheta)
                           Vy=Vy-(Vy*FcosTheta+hy*SinTheta)
                           Vz=Vz-(Vz*FcosTheta+hz*SinTheta)
                           InputParticle%PhaseSpace%Vx=Vx
                           InputParticle%PhaseSpace%Vy=Vy
                           InputParticle%PhaseSpace%Vz=Vz 
                      End If 
                   return
                end subroutine ElectronPostVelocity   
       end subroutine SelectCollisionElectron  
End   Module ElectronMCCModule