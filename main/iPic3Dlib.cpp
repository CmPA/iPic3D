
#include "iPic3D.h"

using namespace iPic3D;

int c_Solver::Init(int argc, char **argv) {
  // initialize MPI environment
  // nproOBcs = number of processors
  // myrank = rank of tha process*/
  mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;

  col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  verbose = col->getVerbose();
  restart_cycle = col->getRestartOutputCycle();
  SaveDirName = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart = col->getRestart_status();
  ns = col->getNs();            // get the number of particle species involved in simulation
  first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  // initialize the virtual cartesian topology 
  vct = new VCtopology3D(col);
  // Check if we can map the processes into a matrix ordering defined in Collective.cpp
  if (nprocs != vct->getNprocs()) {
    if (myrank == 0) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << "x" << vct->getZLEN() << " matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()" << endl;
      mpi->finalize_mpi();
      return (1);
    }
  }
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology(MPI_COMM_WORLD);
  // initialize the central cell index

#ifdef BATSRUS
  // set index offset for each processor
  col->setGlobalStartIndex(vct);
#endif

  nx0 = col->getNxc() / vct->getXLEN(); // get the number of cells in x for each processor
  ny0 = col->getNyc() / vct->getYLEN(); // get the number of cells in y for each processor
  nz0 = col->getNzc() / vct->getZLEN(); // get the number of cells in z for each processor
  // Print the initial settings to stdout and a file
  if (myrank == 0) {
    mpi->Print();
    vct->Print();
    col->Print();
    col->save();
  }
  // Create the local grid
  MPI_Barrier(MPI_COMM_WORLD);
  grid = new Grid3DCU(col, vct);  // Create the local grid
  EMf = new EMfields3D(col, grid);  // Create Electromagnetic Fields Object

  if (col->getSolInit()) {
    /* -------------------------------------------- */
    /* If using parallel H5hut IO read initial file */
    /* -------------------------------------------- */
    ReadFieldsH5hut(ns, false, EMf, col, vct, grid);

  }
  else {
    /* --------------------------------------------------------- */
    /* If using 'default' IO initialize fields depending on case */
    /* --------------------------------------------------------- */
    if      (col->getCase()=="GEMnoPert") EMf->initGEMnoPert(vct,grid,col);
    else if (col->getCase()=="ForceFree") EMf->initForceFree(vct,grid,col);
    else if (col->getCase()=="GEM")       EMf->initGEM(vct, grid,col);
    else if (col->getCase()=="BATSRUS")   EMf->initBATSRUS(vct,grid,col);
    else if (col->getCase()=="Dipole")    EMf->init(vct,grid,col);
    else if (col->getCase()=="ByPert")    EMf->initByPert(vct,grid,col);
    else if (col->getCase()=="ByPert_NoEq")    EMf->initByPert_NoEq(vct,grid,col);
    else if (col->getCase()=="BxPert")    EMf->initBxPert(vct,grid,col);
    else if (col->getCase()=="ExPert")    EMf->initExPert(vct,grid,col);
    else if (col->getCase()=="NPert")    EMf->initNPert(vct,grid,col);
    else if (col->getCase()=="HarrisDP_noPert")    EMf->initDoublePeriodicHarrisNoPerturbation(vct,grid,col);
    else if (col->getCase()=="FF")       EMf->initForceFreePert(vct,grid,col);
    else if (col->getCase()=="FF_kappa") EMf->initForceFreePert(vct,grid,col);
    else if (col->getCase()=="Asymmetric") EMf->initAsymmetric(vct,grid,col);
    // this is Papini's turbulence init for Alfredo, perturbation in the xy plane
    else if (col->getCase()=="alfredo")   EMf->alfredo_turbulence(vct,grid,col);
    // this is Papini's turbulence init for Alfredo, perturbation in the yz plane
    else if (col->getCase()=="alfredo_yz")   EMf->alfredo_turbulence_yz(vct,grid,col);
    else if (col->getCase()=="DoubleGEM")   EMf->initDoubleGEM(vct,grid,col);
    else if (col->getCase()=="HarrisDP_lowerGEM")   EMf->initDP_lowerGEMPerturbed_upperUnperturbed(vct,grid,col);
    else {
      if (myrank==0) {
        cout << " =========================================================== " << endl;
        cout << " WARNING: The case '" << col->getCase() << "' was not recognized. " << endl;
        cout << "          Runing simulation with the default initialization. " << endl;
        cout << " =========================================================== " << endl;
      }
      EMf->init(vct,grid,col);
    }
  }

  // OpenBC
  EMf->updateInfoFields(grid,vct,col);

  // Allocation of particles
  part = new Particles3D[ns];
  if (col->getSolInit()) {
    if (col->getPartInit()=="File") ReadPartclH5hut(ns, part, col, vct, grid);
    else {
      if (myrank==0) cout << "WARNING: Particle drift velocity from ExB " << endl;
      for (int i = 0; i < ns; i++){
        part[i].allocate(i, 0, col, vct, grid);
        if (col->getPartInit()=="EixB") part[i].MaxwellianFromFields(grid, EMf, vct);
        else                            part[i].maxwellian(grid, EMf, vct);


      }
    }
  }
  else {
    for (int i = 0; i < ns; i++)
      part[i].allocate(i, 0, col, vct, grid);

    // Initial Condition for PARTICLES if you are not starting from RESTART
    if (restart == 0) {
      // wave = new Planewave(col, EMf, grid, vct);
      // wave->Wave_Rotated(part); // Single Plane Wave
      for (int i = 0; i < ns; i++)
        if      (col->getCase()=="ForceFree") part[i].force_free(grid,EMf,vct);
        else if (col->getCase()=="BATSRUS")   part[i].MaxwellianFromFluid(grid,EMf,vct,col,i);
        else if (col->getPartInit()=="MaxWhistler") part[i].maxwellian_WhistlerCurrent(grid, EMf, vct);
        else if (col->getCase()=="HarrisDP_noPert") part[i].maxwellian_HarrisDoublePeriodic(grid, EMf, vct);
        else if (col->getPartInit()=="Max_EBDeformed") part[i].maxwellian_EBDeformationWithoutInteraction(grid, EMf, vct);
        else if (col->getCase()=="FF") part[i].Maxwellianspacedist(col,grid,EMf,vct);
        else if (col->getCase()=="FF_kappa") part[i].Kappaspacedist(col,grid,EMf,vct);
        else if (col->getCase()=="Asymmetric") part[i].MaxwellianAsymmetric(col,grid,EMf,vct);
        else if (col->getPartInit()=="Kappa") part[i].kappa(grid,EMf,vct);
        // this is Papini's turbulence initialization for Alfredo, perturbation in the xy plane
	else if (col->getCase()=="alfredo")   part[i].alfredoturbulence(grid, EMf, vct,col);
        // this is Papini's turbulence initialization for Alfredo, perturbation in the yz plane
	else if (col->getCase()=="alfredo_yz")   part[i].alfredoturbulence_yz(grid, EMf, vct,col);
        else if (col->getCase()=="DoubleGEM")   part[i].maxwellian_DoubleGEM(grid, EMf, vct);
	else if (col->getCase()=="HarrisDP_lowerGEM")   part[i].maxwellian_DoubleGEM(grid, EMf, vct);
	else                                  part[i].maxwellian(grid, EMf, vct);
	
    }
  }


  /// initialization of OutputChunks (donw before setting up printing methods)
  outputChunk = new OutputChunks(col, grid);
  
  ///

  if (col->getWriteMethod() == "default") {
    // Initialize the output (simulation results and restart file)
    // PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    // myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
    hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
    if (outputChunk->PrintingChunks){
      // add the outputChunks simulation pointer
      hdf5_agent.set_OutputChunks_simulation_pointer(outputChunk);
    }
    for (int i = 0; i < ns; ++i)
      hdf5_agent.set_simulation_pointers_part(&part[i]);
    output_mgr.push_back(&hdf5_agent);  // Add the HDF5 output agent to the Output Manager's list
    if (myrank == 0 & restart < 2) {
      hdf5_agent.open(SaveDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
      hdf5_agent.open(RestartDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
    }
    // Restart
    num_proc << myrank;
    if (restart == 0) {           // new simulation from input file
      hdf5_agent.open(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
    else {                        // restart append the results to the previous simulation 
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
  }

  for (int i=0; i< ns; i++)
    if (part[i].GetTrackSpecies()){
      part[i].AssignParticlesID(vct);
    }

  Eenergy, Benergy, TOTenergy = 0.0, TOTmomentum = 0.0;
  Ke = new double[ns];
  qx = new double[ns];
  qy = new double[ns];
  qz = new double[ns];
  momentum = new double[ns];
  cq = SaveDirName + "/ConservedQuantities.txt";
  if (myrank == 0) {
    ofstream my_file(cq.c_str());
    my_file.close();
  }
  
  // // Distribution functions
  // nDistributionBins = 1000;
  // VelocityDist = new unsigned long[nDistributionBins];
  // ds = SaveDirName + "/DistributionFunctions.txt";
  // if (myrank == 0) {
  //   ofstream my_file(ds.c_str());
  //   my_file.close();
  // }

  cqsat = SaveDirName + "/VirtualSatelliteTraces" + num_proc.str() + ".txt";
  // if(myrank==0){
  ofstream my_file(cqsat.c_str(), fstream::binary);
  nsat = 3;
  for (int isat = 0; isat < nsat; isat++) {
    for (int jsat = 0; jsat < nsat; jsat++) {
      for (int ksat = 0; ksat < nsat; ksat++) {
        int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
        int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
        int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
        my_file << grid->getXC(index1, index2, index3) << "\t" << grid->getYC(index1, index2, index3) << "\t" << grid->getZC(index1, index2, index3) << endl;
      }}}
  my_file.close();

  Qremoved = new double[ns];

  my_clock = new Timing(myrank);

#ifdef EB
  ebox= new EBox(col);
  if (myrank==0){
    cout <<"Expanding box initialised: U_EB_0= " << ebox->getUEB_0() << ", R_EB_0=" << ebox->getREB_0() << endl;
  }

  cycle_reverseEBdir= col->getCycle_reverseEBdir();
#endif

  /* for VDF */
  /* whether to write VDF, vx, vy, vz */

  WriteVDF_xyz= col->getWriteVDF_xyz();

  if (WriteVDF_xyz){
  
    // Distribution functions
    nDistributionBins = 200;
    VelocityDist_dir = new long double[nDistributionBins];
    
    ds_vx= new string [ns];
    ds_vy= new string [ns];
    ds_vz= new string [ns];
    
    for (int i=0; i< ns; i++){

      std::ostringstream os;
      os << i;
      
      ds_vx[i] = SaveDirName + "/VDF_sp_" + os.str() + "_vx.txt";
      ds_vy[i] = SaveDirName + "/VDF_sp_" + os.str() + "_vy.txt";
      ds_vz[i] = SaveDirName + "/VDF_sp_" + os.str() + "_vz.txt";
    }
    if (myrank == 0) {
      for (int i=0; i< ns; i++){
	ofstream my_file1(ds_vx[i].c_str());
	my_file1.close();
	ofstream my_file2(ds_vy[i].c_str());
        my_file2.close();
	ofstream my_file3(ds_vz[i].c_str());
        my_file3.close();	
	}
    }
  }// end if (WriteVDF_xyz)

  if (outputChunk->PrintingChunks){ // here read if to do this

    //num_proc << myrank; 
    hdf5_agent.open(SaveDirName + "/Chunks_"+ num_proc.str() + ".hdf");
    output_mgr.output("printInfoCh", 0);
    hdf5_agent.close();
  }
  
  return 0;
}

void c_Solver::GatherMoments(){
  // timeTasks.resetCycle();
  // interpolation
  // timeTasks.start(TimeTasks::MOMENTS);

  EMf->updateInfoFields(grid,vct,col);
  EMf->setZeroDensities();                  // set to zero the densities

  for (int i = 0; i < ns; i++)
    part[i].interpP2G(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)

  EMf->sumOverSpecies(vct);                 // sum all over the species
  //
  // Fill with constant charge the planet
  if (col->getCase()=="Dipole") {
    EMf->ConstantChargePlanet(grid, vct, col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
  }

  // EMf->ConstantChargeOpenBC(grid, vct);     // Set a constant charge in the OpenBC boundaries

}

void c_Solver::UpdateCycleInfo(int cycle) {

  EMf->UpdateFext(cycle);
  if (myrank == 0) cout << " Fext = " << EMf->getFext() << endl;
  if (cycle == first_cycle) {
    if (col->getCase()=="Dipole") {
      EMf->SetDipole_2Bext(vct,grid,col);
      EMf->SetLambda(grid, cycle);
    }
  }

#ifdef EB

  double UEB0;
  double REB0;
  if (cycle_reverseEBdir >0){
    if ((cycle % cycle_reverseEBdir ==0) and (cycle != first_cycle)){

      if (myrank == 0) {
	cout << "Cycle " << cycle <<": I am reverting expansion "<< endl;
	cout << "Before reverting: U_0: " << ebox->getUEB_0() << ", R_0: " << ebox->getREB_0() << endl;
      }
      // reverse direction of R0 or U0 in EBox
      ebox->reverseEBdirection();
      // get updated values from EBox
      UEB0= ebox->getUEB_0();
      REB0= ebox->getREB_0();
      // set them in fields 
      EMf->set_UEB_0(UEB0);
      EMf->set_REB_0(REB0);
      // and particles
      for (int i = 0; i < ns; i++){
	part[i].set_UEB_0(UEB0);
	part[i].set_REB_0(REB0);
      }
      if (myrank == 0) {
	cout << "After reverting: U_0: " << ebox->getUEB_0() << ", R_0: " << ebox->getREB_0() << endl;
      }
    }// end if ((cycle % cycle_reverseEBdir ==0) and (cycle != first_cycle)){
  } // end if (cycle_reverseEBdir >0){
  
  // so it enters, as old magnetic field, in the En+theta calculation
  ebox->UpdateEbParameter();
  EMf->UpdateEBVectors(grid, ebox);

  // try damping E the first gyroperiod
  //EMf->SetLambda(grid, cycle);

#endif


}

void c_Solver::CalculateField(int cycle) {

  // timeTasks.resetCycle();
  // interpolation
  // timeTasks.start(TimeTasks::MOMENTS);

  EMf->interpDensitiesN2C(vct, grid);       // calculate densities on centers from nodes
  EMf->calculateHatFunctions(grid, vct);    // calculate the hat quantities for the implicit method
  MPI_Barrier(MPI_COMM_WORLD);
  // timeTasks.end(TimeTasks::MOMENTS);

  // MAXWELL'S SOLVER
  // timeTasks.start(TimeTasks::FIELDS);
  EMf->calculateE(grid, vct, col);               // calculate the E field
  
  // timeTasks.end(TimeTasks::FIELDS);

}

void c_Solver::CalculateBField() {
  /* --------------------- */
  /* Calculate the B field */
  /* --------------------- */

  // timeTasks.start(TimeTasks::BFIELD);
#ifdef EB
  EMf->calculateB_EB(grid, vct, col); 
#else
  EMf->calculateB(grid, vct, col);   // calculate the B field
#endif
  // timeTasks.end(TimeTasks::BFIELD);

  // print out total time for all tasks
  // timeTasks.print_cycle_times();
}

bool c_Solver::ParticlesMover() {

  /*  -------------- */
  /*  Particle mover */
  /*  -------------- */

  // timeTasks.start(TimeTasks::PARTICLES);
  for (int i = 0; i < ns; i++)  // move each species
  {
    // #pragma omp task inout(part[i]) in(grid) target_device(booster)
    
#ifdef EB
    mem_avail= part[i].mover_PC_EB(grid, vct, EMf, ebox); // no subcycling, Expanding Box
#else
    mem_avail = part[i].mover_PC_sub(grid, vct, EMf); // use the Predictor Corrector scheme 
#endif
  }
  // timeTasks.end(TimeTasks::PARTICLES);

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */

  InjectBoundaryParticles();

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  return (false);

}

void c_Solver::InjectBoundaryParticles(){

    /* --------------------------------------- */
    /* Remove particles from depopulation area */
    /* --------------------------------------- */

    if (col->getCase()=="Dipole") {
      for (int i=0; i < ns; i++)
          Qremoved[i] = part[i].deleteParticlesInsideSphere(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());

    }

    /* ------------------------------------------------------------------------ */
    /* Remove all old particles and inject new ones only in the injeciton faces */
    /* ------------------------------------------------------------------------ */

    for (int i=0; i < ns; i++)
      if (col->getRHOinject(i)>0.0){
        mem_avail = part[i].particle_repopulator(grid,vct,EMf,i);

  }

}

void c_Solver::WriteRestart(int cycle) {
  // write the RESTART file
  if (cycle % restart_cycle == 0 && cycle != first_cycle) {
    if (col->getWriteMethod() != "h5hut") {
      // without ,0 add to restart file
      writeRESTART(RestartDirName, myrank, cycle, ns, mpi, vct, col, grid, EMf, part, 0);
    }
  }

}

void c_Solver::WriteConserved(int cycle) {
  // write the conserved quantities
  // regarding the second part of the statement: I want to print only "real" first cycle,
  // not the first restart cycle
  if (cycle % col->getDiagnosticsOutputCycle() == 0 || cycle == 0) {
    Eenergy = EMf->getEenergy();
    Benergy = EMf->getBenergy();
    TOTenergy = 0.0;
    TOTmomentum = 0.0;
    for (int is = 0; is < ns; is++) {
      Ke[is] = part[is].getKe();
      qx[is] = part[is].getqx();
      qy[is] = part[is].getqy();
      qz[is] = part[is].getqz();
      TOTenergy += Ke[is];
      momentum[is] = part[is].getP();
      TOTmomentum += momentum[is];
    }
    if (myrank == 0) {
      ofstream my_file(cq.c_str(), fstream::app);
      my_file << cycle << "\t" << "\t" << (Eenergy + Benergy + TOTenergy) << "\t" << TOTmomentum << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTenergy  << "\t";
      for (int is = 0; is < ns; is++) {
	my_file << Ke[is] << "\t" ;
      }
      for (int is = 0; is < ns; is++) {
        my_file << qx[is] << "\t" ;
	my_file << qy[is] << "\t" ;
	my_file << qz[is] << "\t" ;
      }
      my_file << endl;
      my_file.close();
    }
    // // Velocity distribution
    // for (int is = 0; is < ns; is++) {
    //   double maxVel = part[is].getMaxVelocity();
    //   VelocityDist = part[is].getVelocityDistribution(nDistributionBins, maxVel);
    //   if (myrank == 0) {
    //     ofstream my_file(ds.c_str(), fstream::app);
    //     my_file << cycle << "\t" << is << "\t" << maxVel;
    //     for (int i = 0; i < nDistributionBins; i++)
    //       my_file << "\t" << VelocityDist[i];
    //     my_file << endl;
    //     my_file.close();
    //   }
    // }

  } // end if (cycle % col->getDiagnosticsOutputCycle() == 0 || cycle == 0) {

  if (cycle % col->getFieldOutputCycle() == 0 || cycle == 0) {  
    if (WriteVDF_xyz){
      double minVel, maxVel;
      for (int is = 0; is < ns; is++) {
	// vx
	minVel= part[is].getMinVelocity_x();
	maxVel= part[is].getMaxVelocity_x();
	
	VelocityDist_dir = part[is].getVelocityDistribution_x(nDistributionBins, minVel, maxVel);
	
	if (myrank == 0) {
	  ofstream my_file(ds_vx[is].c_str(), fstream::app);
	  my_file << cycle << "\t"  << minVel  << "\t" << maxVel << "\t" << nDistributionBins << endl;
	  for (int i = 0; i < nDistributionBins; i++){
	    my_file  << VelocityDist_dir[i] <<"\t";
	  }
	  my_file << endl;
	  my_file.close();
	} // end if (myrank == 0) {

	// vy
	minVel= part[is].getMinVelocity_y();
	maxVel= part[is].getMaxVelocity_y();
	
	VelocityDist_dir = part[is].getVelocityDistribution_y(nDistributionBins, minVel, maxVel);
	
	if (myrank == 0) {
	  ofstream my_file(ds_vy[is].c_str(), fstream::app);
	  my_file << cycle << "\t"  << minVel  << "\t" << maxVel << "\t" << nDistributionBins << endl;
	  for (int i = 0; i < nDistributionBins; i++)
	    my_file << "\t" << VelocityDist_dir[i];
	  my_file << endl;
	  my_file.close();
	} // end if (myrank == 0) {

	// vz
	minVel= part[is].getMinVelocity_z();
	maxVel= part[is].getMaxVelocity_z();
	
	VelocityDist_dir = part[is].getVelocityDistribution_z(nDistributionBins, minVel, maxVel);
	
	if (myrank == 0) {
	  ofstream my_file(ds_vz[is].c_str(), fstream::app);
	  my_file << cycle << "\t"  << minVel  << "\t" << maxVel << "\t" << nDistributionBins << endl;
	  for (int i = 0; i < nDistributionBins; i++)
	    my_file << "\t" << VelocityDist_dir[i];
	  my_file << endl;
	  my_file.close();
	} // end if (myrank == 0) {
	
	
      } // end for (int is = 0; is < ns; is++) {
    } // end if (WriteVDF_xyz){

  }
    
  
  
  //if (cycle%(col->getFieldOutputCycle())==0){
  //  for (int is = 0; is < ns; is++) {
  //    part[is].Add_vDist3D();
  //    part[is].Write_vDist3D(SaveDirName);
  //  }
  //}
}

void c_Solver::WriteChunks(int cycle) {

  if (outputChunk->PrintingChunks){
    // packs the field chunks                                                
  //outputChunk->PackChunks(EMf);                                 
  //cout << "End of PackChunks" << endl;             
  // and now the print 
    hdf5_agent.open_append(SaveDirName + "/Chunks_"+ num_proc.str() + ".hdf");
    output_mgr.output("Chunks", cycle);
    hdf5_agent.close();
  }
}


void c_Solver::WriteOutput(int cycle) {

  if (col->getWriteMethod() == "h5hut") {
    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */
    /* -------------------------------------------- */
    if (col->getFieldOutputCycle() != 0){
       if (cycle%(col->getFieldOutputCycle())==0)        WriteFieldsH5hut(ns, grid, EMf,  col, vct, cycle);}
    //if (cycle%(col->getParticlesOutputCycle())==0 &&
    //    cycle!=col->getLast_cycle() && cycle!=0)      WritePartclH5hut(ns, grid, part, col, vct, cycle);
   if (col->getParticlesOutputCycle() != 0)
       if (cycle%(col->getParticlesOutputCycle())==0 )      WritePartclH5hut(ns, grid, part, col, vct, cycle);    
    if (col->getTrackingOutputCycle() != 0) {
       for (int i=0; i< ns; i++){
           if (cycle%(col->getTrackingOutputCycle())==0)      WriteTestPartclH5hut(ns, grid, part, col, vct, cycle, "_track");
       }
    }
  }
  else
  {

    // OUTPUT to large file, called proc**
    if (cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
#ifdef EB
      output_mgr.output("Eall + Ball + rhos + Jsall + pressure + Bext + Flux", cycle);
#else
      output_mgr.output("Eall + Ball + rhos + Jsall + pressure + Flux", cycle);
#endif
      // Pressure tensor is available
      hdf5_agent.close();
    }
    if (cycle % (col->getParticlesOutputCycle()) == 0 && col->getParticlesOutputCycle() != 1) {
      //hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      hdf5_agent.open_append(SaveDirName + "/procPart" + num_proc.str() + ".hdf");  
      output_mgr.output("position + velocity + q ", cycle, 1);
      hdf5_agent.close();
    }
    // write the virtual satellite traces

    if (ns > 2) {
      ofstream my_file(cqsat.c_str(), fstream::app);
      for (int isat = 0; isat < nsat; isat++) {
        for (int jsat = 0; jsat < nsat; jsat++) {
          for (int ksat = 0; ksat < nsat; ksat++) {
            int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
            int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
            int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
            my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
            my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 0) + EMf->getJxs(index1, index2, index3, 2) << "\t" << EMf->getJys(index1, index2, index3, 0) + EMf->getJys(index1, index2, index3, 2) << "\t" << EMf->getJzs(index1, index2, index3, 0) + EMf->getJzs(index1, index2, index3, 2) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 1) + EMf->getJxs(index1, index2, index3, 3) << "\t" << EMf->getJys(index1, index2, index3, 1) + EMf->getJys(index1, index2, index3, 3) << "\t" << EMf->getJzs(index1, index2, index3, 1) + EMf->getJzs(index1, index2, index3, 3) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 0) + EMf->getRHOns(index1, index2, index3, 2) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 1) + EMf->getRHOns(index1, index2, index3, 3) << "\t";
          }}}
      my_file << endl;
      my_file.close();
    }
  }

//  for (int i=0; i< ns; i++){

//    if ((cycle % part[i].GetTrackingSpOutputCycle() == 0 || cycle == first_cycle) && part[i].GetTrackSpecies() ){
//      cout << "Proc " << vct->getCartesian_rank() << " sp " << i << " is writing tracking info" << endl;
//      part[i].WriteTracking(cycle, vct, col, EMf, grid);
      
//    }
//  }

}

void c_Solver::Finalize() {
  if (mem_avail == 0) {          // write the restart only if the simulation finished succesfully
    if (col->getWriteMethod() != "h5hut") {
      writeRESTART(RestartDirName, myrank, (col->getNcycles() + first_cycle) - 1, ns, mpi, vct, col, grid, EMf, part, 0);
    }
  }

  // stop profiling
  my_clock->stopTiming();

  // deallocate
  delete[]Ke;
  delete[]momentum;
  delete[]qx;
  delete[]qy;
  delete[]qz;
  // close MPI
  mpi->finalize_mpi();
}
