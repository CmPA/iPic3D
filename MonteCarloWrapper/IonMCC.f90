Module IonMCC
     Use MCCModule
     Use Constants
     Use ParticleModule 
     Implicit none
     contains
     subroutine SelectCollisionIon(InputParticle,InputGas,InputReaction,NPar,OuputParticle)          
         Implicit none
         Type(MCCParticleOne),intent(inout) :: InputParticle
         Type(ReactionOne),intent(in) :: InputReaction
         Type(Gas),intent(in) :: InputGas
         Integer(4) :: NSpecy,ParticleIndex(2_4)
         Integer(4),intent(out)  :: NPar
         Type(ParticleOne) :: OuputParticle(3)
         Type(ParticleOne) :: TempParticle(2_4) 
         NSpecy=InputGas%NSpecy
         Select case (InputReaction%ReactionType)
                  Case(0_4)
                          NPar=1
                          OuputParticle(1)=InputParticle%PhaseSpace
                          OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                          OuputParticle(1)%CollisionType=-1
                  Case(101_4)
                          Call PostVelocityIonElastic(InputParticle,InputGas,IsotropicCosKai)
                          NPar=1
                          OuputParticle(1)=InputParticle%PhaseSpace
                          OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                          OuputParticle(1)%CollisionType=4
                  Case(102_4)
                          !Call PostVelocityIonElastic(InputParticle,InputGas,ReactiveCosKai)
                          !Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  Case(103_4)  ! ChargeExchange
                          NPar=1
                          OuputParticle(1)=InputParticle%GasParticle
                          OuputParticle(1)%sp_out=OuputParticle(1)%sp_in
                          OuputParticle(1)%CollisionType=3
                          !Call AddParticleSmall(InputParticle%GasParticle,OuputParticle(InputReaction%Reactant))
                  !Case(111_4)
                  !        Call  PostVelocityIonReactive(InputParticle,InputGas,InputReaction%Threshold)
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(InputReaction%Reactant))
                  !Case(112_4)
                  !        Call  PostVelocityIonReactive(InputParticle,InputGas,InputReaction%Threshold)
                  !        Call MCCParticleIndex(NSpecy,InputReaction%Resultant,1_4,ParticleIndex(1))
                  !        Call AddParticleSmall(InputParticle%PhaseSpace,OuputParticle(ParticleIndex(1)))
                  Case(-1_4)
                  end select
                  Call MometumEnergyChange(InputParticle%Mass,OuputParticle(1),InputParticle%Particle_Old)
        return
        contains
    
        subroutine PostVelocityIonElastic(InputParticle,InputGas,CosKai)
               Implicit none
               Type(MCCParticleOne),intent(inout) :: InputParticle
               Type(Gas),intent(in) :: InputGas
               Real(8),external ::  CosKai
               real(8) :: V,MassRatio
               real(8) :: g,gx,gy,gz,gper,gtemp
               real(8) :: Vx,Vy,Vz,hx,hy,hz
               real(8) :: AFai,CosFai,SinFai
               real(8) :: CosTheta,FCosTheta,SinTheta
               MassRatio=InputGas%MGas/(InputGas%MGas+InputParticle%Mass) 

               Vx=InputParticle%PhaseSpace%Vx
               Vy=InputParticle%PhaseSpace%Vy
               Vz=InputParticle%PhaseSpace%Vz 
               gx=Vx-InputParticle%GasParticle%Vx
               gy=Vy-InputParticle%GasParticle%Vy
               gz=Vz-InputParticle%GasParticle%Vz
               gtemp=gy*gy+gz*gz
               gper=Dsqrt(gtemp)
               g=Dsqrt(gtemp+gx*gx)
                If(Abs(gper)<MinReal)  then
                           Call RandomVelocity(V,Vx,Vy,Vz)
               else
                       CosTheta=CosKai(InputParticle%Energy)
                       FcosTheta=1.d0-cosTheta
                       SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)

                       Call DRandom(R)
                       AFai=2.d0*Pi*R
                       CosFai=DCos(AFai)
                       SinFai=DSin(AFai)
                       hx=gper*CosFai
                       hy=-(gx*gy*CosFai+g*gz*SinFai)/gper
                       hz=-(gx*gz*CosFai-g*gy*SinFai)/gper
                       Vx=Vx-MassRatio*(gx*FcosTheta+hx*SinTheta)
                       Vy=Vy-MassRatio*(gy*FcosTheta+hy*SinTheta)
                       Vz=Vz-MassRatio*(gz*FcosTheta+hz*SinTheta)
                       InputParticle%PhaseSpace%Vx=Vx
                       InputParticle%PhaseSpace%Vy=Vy
                       InputParticle%PhaseSpace%Vz=Vz
                 End if 
               return
            end subroutine PostVelocityIonElastic
            
            Subroutine PostVelocityIonReactive(InputParticle,InputGas,Threshold)
            ! Note:Here the velocity of the target before the collision has be set to zero to save time.
               implicit none
               Type(MCCParticleOne),intent(inout) :: InputParticle
               Type(Gas),intent(in) :: InputGas
               Real(8),intent(in) :: Threshold
               real(8) :: MassRatioA,MassRatioB,Miu
               real(8) :: g,gx,gy,gz,Gafter
               real(8) :: Vx,Vy,Vz
               MassRatioA=InputGas%MGas/(InputGas%MGas+InputParticle%Mass)
               MassRatioB=1.d0-MassRatioA
               Miu=MassRatioA*InputParticle%Mass 
                Vx=InputParticle%PhaseSpace%Vx
               Vy=InputParticle%PhaseSpace%Vy
               Vz=InputParticle%PhaseSpace%Vz 
               gx=Vx-InputParticle%GasParticle%Vx
               gy=Vy-InputParticle%GasParticle%Vy
               gz=Vz-InputParticle%GasParticle%Vz
               g=Dsqrt(gx*gx+gy*gy+gz*gz)
               Gafter=Dsqrt(g*g-2.d0*Threshold/Miu)
               Call RandomVelocity(Gafter,gx,gy,gz)  
               InputParticle%PhaseSpace%Vx=MassRatioB*Vx+MassRatioA*InputParticle%GasParticle%Vx+MassRatioB*gx
               InputParticle%PhaseSpace%Vy=MassRatioB*Vy+MassRatioA*InputParticle%GasParticle%Vy+MassRatioB*gy
               InputParticle%PhaseSpace%Vz=MassRatioB*Vz+MassRatioA*InputParticle%GasParticle%Vz+MassRatioB*gz 
               return
            end subroutine PostVelocityIonReactive
  end  subroutine SelectCollisionIon  
end Module IonMCC


