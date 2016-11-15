Module MCCPublic
     Use TypeModule
     Use Constants
     Use ParticleModule
     Use  MCCModule
     Use ElectronMCCModule
    Use  IonMCC
     Implicit none
     Integer(4),parameter :: NSpecyMCCTemp=6_4
     Type(MCCParticleOne),save,private :: TempMCCParticle
     

   Contains
   subroutine  MCCBundleInit(InputGas,NSpecy,InputSpecy,InputReactionBundle,InputSigma,OutputMCCBundle,dt)
             Use SigmaModule
             Implicit none
             Type(Gas),intent(in) ::  InputGas
              Integer(4),intent(in) ::  NSpecy
             Type(OneSpecy),intent(in) :: InputSpecy(0:NSpecy)
             Type(ReactionBundle),intent(in) :: InputReactionBundle(0:1)
             Type(SigmaNormalized),intent(in) :: InputSigma(0:1)
             Type(MCCBundle),intent(out) ::  OutputMCCBundle(0:NSpecy)
             Real(8),intent(in) :: dt 
             Integer(4) ::  i,j,SPIndex
             Real(8) :: Miu,SigmaT
             
             do i=0,NSpecy
                  OutputMCCBundle(i)%Model=InputSpecy(i)%Model
                  OutputMCCBundle(i)%Mass=InputSpecy(i)%Mass
                  OutputMCCBundle(i)%VFactor=InputSpecy(i)%VFactor
            End do

              !Ng=(InputGas%PGas)/(kB*InputGas%TGas)
             
             do i=0, NSpecy
                     SPIndex=InputSpecy(i)%Model
                     OutputMCCBundle(i)%MFactor=InputGas%MFactor
                     OutputMCCBundle(i)%EFactor=InputGas%EFactor
                     
                    OutputMCCBundle(i)%EnergyInterval=InputSigma(SPIndex)%EnergyInterval
                    OutputMCCBundle(i)%NReaction=InputSigma(SPIndex)%NReaction
                    OutputMCCBundle(i)%NSigma=InputSigma(SPIndex)%NSigma
                    OutputMCCBundle(i)%Reaction(0)%Reactant=InputReactionBundle(SPIndex)%Reaction(1)%Reactant
                    do j=1,InputReactionBundle(SPIndex)%NReaction
                            OutputMCCBundle(i)%Reaction(j)=InputReactionBundle(SPIndex)%Reaction(j)
                    End do
                    Miu=InputSpecy(i)%Mass*InputGas%MGas/(InputGas%MGas+InputSpecy(i)%Mass)
                    Call  ProbilityNonReactive(InputSigma(SPIndex)%NReaction,InputSigma(SPIndex)%NSigma,InputSigma(SPIndex)%Value,OutputMCCBundle(i)%Probility,OutputMCCBundle(i)%EnergyInterval,Miu,SigmaT)
                    OutputMCCBundle(i)%CollisionRatio=1.d0-DExp(-SigmaT*InputGas%NGas*dt)
              enddo 
              Write(*,*) OutputMCCBundle%CollisionRatio,'ratio'  
           return
         End subroutine MCCBundleInit 
         
          subroutine  ProbilityNonReactive(NReaction,NSigma,Sigma,Probility,EnergyInterval,Mass,SigmaT)
                implicit none
                Integer(4),intent(in) ::  NReaction,NSigma
                Real(8),intent(in) ::  Sigma(NSigma,NReaction),EnergyInterval,Mass
                Real(8),intent(out) ::  Probility(NReaction,NSigma),SigmaT
                Real(8) ::Energy,V,Max
                Integer(4) :: i,j
                do i=1,NSigma
                     do j=1,NReaction
                          Probility(j,i)=Sigma(i,j)
                     Enddo
                Enddo
                   
                do i=1,NSigma
                      do j=2,NReaction
                            Probility(j,i)=Probility(j,i)+Probility(j-1,i)
                      end do
                end do
                
                Max=0.d0
                do i=1,NSigma
                      Energy=dble(i)*EnergyInterval
                      V=DSqrt(2.d0*Energy*JtoeV/Mass)  
                      do j=1,NReaction
                            Probility(j,i)=Probility(j,i)*V
                            If  (Probility(j,i)>Max) Then
                                  Max=Probility(j,i)
                            End if
                      end do
                end do
                
                do i=1,NSigma
                     do j=1,NReaction
                          Probility(j,i)=Probility(j,i)/Max
                     Enddo
                Enddo
                SigmaT=Max
             return 
          End subroutine  ProbilityNonReactive

      subroutine SelectProbility(InputMCCParticle,InputMCCBundle)
           Implicit None 
           Type(MCCParticleOne),intent(inout) :: InputMCCParticle
           Type(MCCBundle),intent(in) :: InputMCCBundle
           Integer(4) :: i,N,Index,Upper,Center,Lower,NReaction
           Real(8) :: EnergyInterval,Energy,S1,S2
           Real(8) :: TempProbility(InputMCCBundle%NReaction)
            TempProbility=0.d0 
            EnergyInterval=InputMCCBundle%EnergyInterval
            Energy=InputMCCParticle%Energy
            NReaction=InputMCCBundle%NReaction
            N=Int(Energy/EnergyInterval)
            If(N<2) then
                 TempProbility=InputMCCBundle%Probility(1:NReaction)
            else if (N<InputMCCBundle%NSigma) then
                 Index=N*NReaction
                 do i=1,NReaction
                       Center=Index+i
                       Upper=Center+NReaction
                       Lower=Center-NReaction
                       If(Energy<=InputMCCBundle%Reaction(i)%Threshold) Then
                           TempProbility(i:NReaction)=-1.d0
                           exit       
                      else
                           S1=(Energy-(dble(N)*EnergyInterval))/EnergyInterval
                           S2=1.d0-S1 
                           TempProbility(i)=S1*InputMCCBundle%Probility(Lower)+S2*InputMCCBundle%Probility(Upper)
                      end If
                  end do
           else
                  Index=(InputMCCBundle%NSigma-1)*NReaction 
                  TempProbility=InputMCCBundle%Probility(Index+1:Index+NReaction)
           End if
           Call DRandom(R)
           InputMCCParticle%Index=0 
           do i=1, NReaction
               If(TempProbility(i)<MinReal) then  
                   InputMCCParticle%Index=0
                   exit 
              else if(R<TempProbility(i)) then
                   InputMCCParticle%Index=i
                  exit
               End If
           end do 
  End  subroutine SelectProbility
           
        Subroutine UpdateParticleMCCElectron(InputParticle,OutputParticle,Mass)
                     Implicit none
                     Type(ParticleOne),intent(in) :: InputParticle
                     Type(MCCParticleOne),intent(out) :: OutputParticle
                     Real(8),intent(in) :: Mass
                     Real(8) :: Vx,Vy,Vz,V2,V,Energy
                     OutputParticle%PhaseSpace=InputParticle
                     OutputParticle%Particle_Old=InputParticle
                     OutputParticle%Mass=Mass
                     Vx=InputParticle%Vx
                     Vy=InputParticle%Vy
                     Vz=InputParticle%Vz
                     V2=Vx*Vx+Vy*Vy+Vz*Vz
                     V=Dsqrt(V2)
                     Energy=0.5d0*Mass*V2/JtoeV
                     OutputParticle%V2=V2
                     OutputParticle%V=V
                     OutputParticle%Energy=Energy
                    return  
            End subroutine UpdateParticleMCCElectron
                
            Subroutine UpdateParticleMCCIon(InputParticle,OutputParticle,Mass,InputGas)
                     Implicit none
                     Type(ParticleOne),intent(in) :: InputParticle
                     Type(MCCParticleOne),intent(out) :: OutputParticle
                     Real(8),intent(in) :: Mass
                     Type(Gas),intent(in) :: InputGas
                     Real(8) :: Vx,Vy,Vz,V2,V,Energy,Miu
                       OutputParticle%GasParticle=InputParticle
                       OutputParticle%PhaseSpace=InputParticle
                       OutputParticle%Particle_Old=InputParticle
                       OutputParticle%Mass=Mass 
                       Miu=Mass*InputGas%MGas/(InputGas%MGas+Mass)
                       Call Maxwellian(InputGas,OutputParticle%GasParticle)
                       Vx=InputParticle%Vx-OutputParticle%GasParticle%Vx
                       Vy=InputParticle%Vy-OutputParticle%GasParticle%Vy
                       Vz=InputParticle%Vz-OutputParticle%GasParticle%Vz
                       V2=Vx*Vx+Vy*Vy+Vz*Vz
                       V=Dsqrt(V2)
                       Energy=0.5d0*Miu*V2/JtoeV
                       OutputParticle%V2=V2
                       OutputParticle%V=V
                       OutputParticle%Energy=Energy
                    return  
            End subroutine UpdateParticleMCCIon
            

            
          
            
end module MCCPublic