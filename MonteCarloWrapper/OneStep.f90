Module OneStepModule
      Use MCCPublic
      Implicit none
      ! This section defines the particles.
                  Integer(4),parameter :: NSpecyMax=5_4
                  !Type(ParticleBundle),private ::  ParticleGlobal(0:NSpecyMax)
                  Integer(4),save :: NSpecy=3_4
                  Type(Gas),save ::  GasGlobal
                  Type(MCCBundle),save,private ::  MCCBundleGlobal(0:NSpecyMax)
                  !Type(ParticleBundleSmall),save,private ::  MCCBefore(0:NSpecyMax),MCCAfter(0:NSpecyMax)
                  Type(OneSpecy),save :: SPGlobal(0:NSpecyMax)
                  
                  !Real(8) :: CollisionRatio(0:NSpecyMax)
      !  This section defines the MCC and gas properties.
           
    contains
    Subroutine MCCInitialization(Inputdt,InputNSpecy,CollisionRatio) bind(C,name= "MCCInitialization")
             Use GasModule
             Use MCCPublic
             Use ArGasModule
             Use ISO_C_BINDING
             Implicit none
             !Integer(4),intent(in) ::  InputNSpecy
             !Real(8),intent(in):: dt
             !Real(8),intent(out):: CollisionRatio(0:InputNSpecy)
             Real(C_DOUBLE),intent(in):: Inputdt
             Integer(C_int),intent(in) ::  InputNSpecy
             Real(C_DOUBLE),intent(out):: CollisionRatio(0:InputNSpecy)
             !Integer(C_int),intent(in) :: LengthName
             !Character(kind=C_CHAR,Len=LengthName),intent(in) :: Filename
             !Character(99) :: Filename
             Integer(4) ::  i
             Real(8) :: dt
             dt=Inputdt
             !Write(*,*) "InitFortranBefore", InputNSpecy,dt,Inputdt,CollisionRatio(0)
             Call ArSigmaInitialization()
             
             Call GasInit(NSpecyMax,SPGlobal,GasGlobal)
             !Write(*,*) "Re",ArReaction(0),"s22s",ArSigma(0)%Value(1:6)
             NSpecy=InputNSpecy 
             Call MCCBundleInit(GasGlobal,NSpecy,SPGlobal,ArReaction,ArSigma,MCCBundleGlobal,dt)  
             !Write(*,*) "GE",GasGlobal,"b22b",MCCBundleGlobal(0)%Probility(1:6)
             do i=0,NSpecy
                    CollisionRatio(i)=MCCBundleGlobal(i)%CollisionRatio
             end do
             !Write(*,*) "InitFortranAfter", InputNSpecy,dt,Inputdt,CollisionRatio(0),MCCBundleGlobal(0)%CollisionRatio
        return  
    End Subroutine MCCInitialization
    
    Subroutine MCC(NParBefore,MCCBefore,NParAfter,MCCAfter) bind(C,name= "MCC")
        Implicit none
        ! j is the Particle Index;
        Integer(C_int),intent(in) :: NParBefore
        Integer(C_int),intent(inout) :: NParAfter
        Type(ParticleOne),intent(inout) :: MCCBefore(NParBefore)
        Type(ParticleOne),intent(out) :: MCCAfter(NParAfter)
        Integer(4) :: NTP
        Type(ParticleOne) :: TempParticle(3)
        Type(MCCParticleOne) :: TempMCCParticle
        Integer(4) :: i,j,SP,Index
        
        !Type(Gas),intent(in) :: InputGas
        !Type(MCCBundle),intent(in) :: InputMCCBundle
        NParAfter=0
        !Write(*,*) "before colliding",NParBefore,MCCBefore(1),MCCBefore(1)%Vx
               do i=1,NParBefore
                       SP=MCCBefore(i)%SP_in
 
                       MCCBefore(i)%Vx=MCCBefore(i)%Vx*MCCBundleGlobal(SP)%VFactor
                       MCCBefore(i)%Vy=MCCBefore(i)%Vy*MCCBundleGlobal(SP)%VFactor
                       MCCBefore(i)%Vz=MCCBefore(i)%Vz*MCCBundleGlobal(SP)%VFactor
                       
                       select case (MCCBundleGlobal(SP)%Model)
                         Case(0_4)
                             Call UpdateParticleMCCElectron(MCCBefore(i),TempMCCParticle,MCCBundleGlobal(SP)%Mass)
                             Call SelectProbility(TempMCCParticle,MCCBundleGlobal(SP))
                             Index=TempMCCParticle%Index
                             Call  SelectCollisionElectron(TempMCCParticle,GasGlobal,MCCBundleGlobal(SP)%Reaction(Index),NTP,TempParticle)
                         Case(1_4)
                              Call UpdateParticleMCCIon(MCCBefore(i),TempMCCParticle,MCCBundleGlobal(SP)%Mass,GasGlobal) 
                              Call SelectProbility(TempMCCParticle,MCCBundleGlobal(SP))
                              Index=TempMCCParticle%Index
                              Call  SelectCollisionIon(TempMCCParticle,GasGlobal,MCCBundleGlobal(SP)%Reaction(Index),NTP,TempParticle)
                         end select
                         do j=1,NTP
                             NParAfter=NParAfter+1
                             TempParticle(j)%Vx=TempParticle(j)%Vx/MCCBundleGlobal(TempParticle(j)%sp_in)%VFactor
                             TempParticle(j)%Vy=TempParticle(j)%Vy/MCCBundleGlobal(TempParticle(j)%sp_in)%VFactor
                             TempParticle(j)%Vz=TempParticle(j)%Vz/MCCBundleGlobal(TempParticle(j)%sp_in)%VFactor
                             TempParticle(j)%dP_x=TempParticle(j)%dP_x/MCCBundleGlobal(TempParticle(j)%sp_in)%MFactor
                             TempParticle(j)%dP_y=TempParticle(j)%dP_y/MCCBundleGlobal(TempParticle(j)%sp_in)%MFactor
                             TempParticle(j)%dP_z=TempParticle(j)%dP_z/MCCBundleGlobal(TempParticle(j)%sp_in)%MFactor
                             TempParticle(j)%dE=TempParticle(j)%dE/MCCBundleGlobal(TempParticle(j)%sp_in)%EFactor
                             MCCAfter(NParAfter)=TempParticle(j)
                        end do
               End do
               !Write(*,*) "After colliding",NParAfter,MCCAfter(1),MCCAfter(1)%Vx
              return
     End  Subroutine MCC
End Module OneStepModule
