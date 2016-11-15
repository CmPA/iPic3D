Module MCCModule
    Use TypeModule
    USE Constants
    Use  EnergyKai
    Implicit none
    !Integer(4),parameter :: NMCParMax=20000_4
  contains
   Subroutine UpdateVelocityMCC(VFactor,InputParticle)
        Implicit none
        Type(MCCParticleOne),intent(inout) :: InputParticle
         Real(8) :: VFactor
         InputParticle%V=VFactor*InputParticle%V
         InputParticle%PhaseSpace%Vx=VFactor*InputParticle%PhaseSpace%Vx
         InputParticle%PhaseSpace%Vy=VFactor*InputParticle%PhaseSpace%Vy
         InputParticle%PhaseSpace%Vz=VFactor*InputParticle%PhaseSpace%Vz
         return
   End  subroutine UpdateVelocityMCC

    subroutine Maxwellian(InputGas,InputParticle)
       implicit none
       Type(Gas),intent(in) :: InputGas
       Type(ParticleOne),intent(out) :: InputParticle
       real(8) ::  Mass,Temprature 
       real(8) :: V,Beta,FuncA,FuncB
       real(8) :: Theta,CosTheta,SinTheta,Fai 
       Mass=InputGas%MGas
       Temprature=InputGas%TGas
       Beta=1.d0/(kB*Temprature)
       FuncA=1.d0
       FuncB=0.d0
      do while(FuncA>FuncB)
         Call DRandom(R)
              FuncA=R*R
         Call DRandom(R)
              FuncB=-exp*R*Dlog(R)
      end do
      V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)
      Call RandomVelocity(V,InputParticle%Vx,InputParticle%Vy,InputParticle%Vz)
     return 
   end subroutine Maxwellian
  
   subroutine RandomVelocity(V,Vx,Vy,Vz)
       implicit none
       Real(8),intent(in) ::  V
       Real(8),intent(out) ::  Vx,Vy,Vz  
       Real(8) :: Fai,CosFai,SinFai
       Real(8) :: Theta,CosTheta,FcosTheta,SinTheta
       Call DRandom(R)
       CosTheta=IsotropicCosKai(Theta)
        SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)
        Call DRandom(R)
        Fai=2.d0*PI*R
        Vx=V*CosTheta
        Vy=V*SinTheta*DCos(Fai)
        Vz=V*SinTheta*Dsin(Fai)
       return
   end subroutine RandomVelocity   
  
Subroutine MCCParticleIndex(NSpecy,Resultant,NIndex,ParticleIndex)
            Implicit none
            Integer(4),Intent(in) ::  NSpecy,Resultant,NIndex
            Integer(4),Intent(Out) ::  ParticleIndex(NIndex)
            Integer(4) :: i,j
            j=1 
            do i=0,NSpecy
                 If (Btest(Resultant,i)) then
                      ParticleIndex(j)=i
                      ! If(j=NIndex)  exit
                      j=j+1
                End if
             End do        
            return 
    End subroutine MCCParticleIndex
End Module MCCModule   

