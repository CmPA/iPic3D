
#include "iPic3D.h"
#include "MyClock.h"

using namespace iPic3D;

extern MyClock *clocks;

int c_Solver::Init(int argc, char **argv) {
  // initialize MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;

  clocks = new MyClock(6);

  clocks->start(0);
  clocks->start(5);

  col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  verbose = col->getVerbose();
  restart_cycle = col->getRestartOutputCycle();
  SaveDirName = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart = col->getRestart_status();
  ns = col->getNs();            // get the number of particle species involved in simulation
  first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  qom  = new double[ns];
  for (int i=0; i < ns; i++)
      qom[i] = col->getQOM(i);

  x_center = col->getx_center();
  y_center = col->gety_center();
  z_center = col->getz_center();
  L_square = col->getL_square();
  L_outer = col->getL_outer();
  cylindrical = col->getcylindrical();

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
    if (col->getCase()=="KAWTurbulencePert") {
      double mime = fabs(col->getQOM(0)/col->getQOM(1));
      double TiTe = pow(col->getUth(1)/col->getUth(0), 2.0)*mime;
      EMf->initKAWTurbulencePert(vct, grid, col, mime, TiTe);
    }
    else if (col->getCase()=="DoubleHarrisRel_pairs") EMf->initDoubleHarrisRel_pairs(vct,grid,col);
    else if (col->getCase()=="DoubleHarrisRel_ionel") EMf->initDoubleHarrisRel_ionel(vct,grid,col);
    // OLD CASES FROM IPIC
//    else if (col->getCase()=="GEMnoPert") EMf->initGEMnoPert(vct,grid,col);
//    else if (col->getCase()=="ForceFree") EMf->initForceFree(vct,grid,col);
//    else if (col->getCase()=="ForceFreeHump") EMf->initForceFreeWithGaussianHumpPerturbation(vct,grid,col);
//    else if ((col->getCase()=="GEM") || (col->getCase()=="GEMRelativity"))  EMf->initGEM(vct, grid,col);
//    else if (col->getCase()=="KAWTurbulencePert") {
//      double mime = fabs(col->getQOM(0)/col->getQOM(1));
//      double TiTe = pow(col->getUth(1)/col->getUth(0), 2.0)*mime;
//      EMf->initKAWTurbulencePert(vct, grid, col, mime, TiTe);
//    }
//    else if (col->getCase()=="HarrisSteps")       EMf->initDoublePeriodicHarrisSteps(vct, grid,col);
//    else if (col->getCase()=="BATSRUS")   EMf->initBATSRUS(vct,grid,col);
//    else if (col->getCase()=="Dipole")    EMf->init(vct,grid,col);
//    else if (col->getCase()=="DoubleHarris")    EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
//    else if (col->getCase()=="Whistler")    EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
//    else if (col->getCase()=="WhistlerKappa")    EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
//    else if (col->getCase()=="Coils")  EMf->initWB8(vct,grid,col);
//    else if (col->getCase()=="CoilsMono")  EMf->initWB8(vct,grid,col);
//    else if (col->getCase()=="TwoCoils")  EMf->initTwoCoils(vct,grid,col);
//    else if (col->getCase()=="FluxRope")  EMf->initFluxRope(vct,grid,col);
//    else if (col->getCase()=="GEMNoVelShear")  EMf->initHarrisNoVelShear(vct, grid,col);
//    else if (col->getCase()=="Relativistic")  EMf->init(vct, grid, col);
    else {
      if (myrank==0) {
        cout << " =========================================================== " << endl;
        cout << " WARNING: The case '" << col->getCase() << "' was not recognized. " << endl;
        cout << "          Running simulation with the default initialization. " << endl;
        cout << " =========================================================== " << endl;
      }
      EMf->init(vct,grid,col);
    }
  }

  #ifdef __PETSC_SOLVER__
    // PETSc solver:
    petscSolver = new PetscSolver(EMf, grid, vct, col);
  #endif

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
    for (int i = 0; i < ns; i++) part[i].allocate(i, 0, col, vct, grid);

    // Initial Condition for PARTICLES if you are not starting from RESTART
    if (restart == 0) {
      // wave = new Planewave(col, EMf, grid, vct);
      // wave->Wave_Rotated(part); // Single Plane Wave

      if (myrank==0) cout << col->getCase() << endl;

      for (int i = 0; i < ns; i++) {
        if (col->getCase()=="KAWTurbulencePert") {
          double mime = fabs(col->getQOM(0)/col->getQOM(1));
          double TiTe = pow(col->getUth(1)/col->getUth(0), 2.0)*mime;
          part[i].KAWTurbulencePert(grid, EMf, vct, col->getB0x(), mime, TiTe, col->getPartSymmetric());
        }
        else if (col->getCase()=="DoubleHarrisRel_pairs") part[i].DoubleHarrisRel_pairs(grid, EMf, vct, col);
        else if (col->getCase()=="DoubleHarrisRel_ionel") part[i].DoubleHarrisRel_ionel(grid, EMf, vct, col);
//	// OLD CASES FROM IPIC
//        else if (col->getCase()=="ForceFree") part[i].force_free(grid,EMf,vct);
//        else if (col->getCase()=="BATSRUS")   part[i].MaxwellianFromFluid(grid,EMf,vct,col,i);
//        else if (col->getCase()=="DoubleHarris")    part[i].maxwellian_reversed(grid, EMf, vct);
//        else if (col->getCase()=="Whistler")    part[i].maxwellian_whistler(grid, EMf, vct);
//        else if (col->getCase()=="WhistlerKappa")    part[i].kappa(grid, EMf, vct);
//        else if (col->getCase()=="GEMRelativity")    part[i].relativistic_maxwellian(grid, EMf, vct);
//        else if (col->getCase()=="Relativistic")  part[i].twostream1D(grid, vct, 3);
//        else if (col->getCase()=="GEM" || col->getCase()=="GEMNoVelShear") {
//          if (i<2)
//            part[i].maxwellian(grid, EMf, vct);
//          else
//            if(col->getPartInit()=="Kappa")
//              part[i].kappa(grid, EMf, vct);
//            else
//              part[i].maxwellian(grid, EMf, vct);
//        }
//        else if (col->getCase()=="Coils") {
//          if (col->getRHOinit(i) > 0.0)
//            part[i].maxwell_box(grid,EMf,vct,L_square,x_center,y_center,z_center, 1.0); //generates maxwellian in a box
//          else
//               part[i].empty(grid, EMf, vct);
//        }
//        else if (col->getCase()=="TwoCoils"){
//          if (col->getRHOinit(i) > 0.0)
//            part[i].maxwell_box(grid,EMf,vct,L_square,x_center,y_center,z_center, 1.0); //generates maxwellian in a box
//          else
//            part[i].empty(grid, EMf, vct);
//        }
//        else if (col->getCase()=="CoilsMono"){
//          if (col->getRHOinit(i) > 0.0)
//            part[i].monoenergetic_box(grid,EMf,vct,L_square,x_center,y_center,z_center, 1.0); //generates maxwellian in a box
//          else
//            part[i].empty(grid, EMf, vct);
//        }
        else part[i].maxwellian(grid, EMf, vct);
      }
    }
  }

  num_proc << myrank;
  if (col->getWriteMethod() == "default") {
    // Initialize the output (simulation results and restart file)
    // PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    // myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
    hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
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

  Eenergy, Benergy, TOTenergy = 0.0, TOTmomentum = 0.0;
  Ke = new double[ns];
  BulkEnergy = new double[ns];
  momentum = new double[ns];
  Qtot = new double[ns]; Qtot[0] = Qtot[1] = 0;
  totParticles = new int[ns];
  globalTotParticles = new int[ns];
  speciesTemp = new double[ns];

  cq = SaveDirName + "/ConservedQuantities.txt";
  cq2 = SaveDirName + "/SummaryQuantities.txt";

  if (myrank == 0) {
    ofstream my_file(cq.c_str(), fstream::app);
    my_file.close();
    ofstream my_file2(cq2.c_str(), fstream::app);
    my_file2.close();

  }
  
  // // Distribution functions
  // nDistributionBins = 1000;
  // VelocityDist = new unsigned long[nDistributionBins];
  // ds = SaveDirName + "/DistributionFunctions.txt";
  // if (myrank == 0) {
  //   ofstream my_file(ds.c_str());
  //   my_file.close();
  // }

  //num_proc << myrank;
  cqsat = SaveDirName + "/VirtualSatelliteTraces" + num_proc.str() + ".txt";
  // if(myrank==0){
/*  ofstream my_file(cqsat.c_str(), fstream::binary);

  nsatx=3;
  nsaty=3;
  nsatz=3;
  if(vct->getXLEN() == 1){
      nsatx=1;
  }
  if(vct->getYLEN() == 1){
      nsaty=1;
  }
  if(vct->getZLEN() == 1){
      nsatz=1;
  }
  my_file << nsatx*nsaty*nsatz << "\t" <<nsatx << "\t" <<nsaty << "\t" <<nsatz << endl;
  for (int isat=0; isat < nsatx; isat++){
      for (int jsat=0; jsat < nsaty; jsat++){
          for (int ksat=0; ksat < nsatz; ksat++){
              int index1 = 1+isat*nx0/nsatx+nx0/nsatx/2;
              int index2 = 1+jsat*ny0/nsaty+ny0/nsaty/2;
              int index3 = 1+ksat*nz0/nsatz+nz0/nsatz/2;
            my_file <<  grid->getXC(index1,index2,index3) << "\t" << grid->getYC(index1,index2,index3) << "\t" << grid->getZC(index1,index2,index3) << endl;
      }}}
  my_file.close();
*/


  Qremoved = new double[ns];

  my_clock = new Timing(myrank);

  clocks->stop(0);
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
}

void c_Solver::UpdateCycleInfo(int cycle) {

  if (col->getCase()=="Dipole") EMf->UpdateFext(cycle);
  if (myrank == 0) cout << " Fext = " << EMf->getFext() << endl;
  if (cycle == first_cycle) {
    if (col->getCase()=="Dipole") {
      EMf->SetDipole_2Bext(vct,grid,col);
      EMf->SetLambda(grid);
    }
  }


}

void c_Solver::CalculateField() {

  // timeTasks.resetCycle();
  // interpolation
  // timeTasks.start(TimeTasks::MOMENTS);

  EMf->interpDensitiesN2C(vct, grid);       // calculate densities on centers from nodes
  EMf->calculateHatFunctions(grid, vct, col);    // calculate the hat quantities for the implicit method
  MPI_Barrier(MPI_COMM_WORLD);
  // timeTasks.end(TimeTasks::MOMENTS);

  // MAXWELL'S SOLVER
  // timeTasks.start(TimeTasks::FIELDS);
  #ifdef __PETSC_SOLVER__
    petscSolver->solveE();
  #else
    EMf->calculateE(grid, vct, col);               // calculate the E field
  #endif
  // timeTasks.end(TimeTasks::FIELDS);

}

void c_Solver::CalculateBField() {
  /* --------------------- */
  /* Calculate the B field */
  /* --------------------- */

  // timeTasks.start(TimeTasks::BFIELD);
  EMf->calculateB(grid, vct, col);   // calculate the B field
  // timeTasks.end(TimeTasks::BFIELD);

  // print out total time for all tasks
  // timeTasks.print_cycle_times();
}

bool c_Solver::PushParticles() {

  /*  ------------------------ */
  /*  Particle position pusher */
  /*  ------------------------ */

  // timeTasks.start(TimeTasks::PARTICLES);
  for (int i = 0; i < ns; i++) {  // move each species
    mem_avail = part[i].push_particles(grid, vct);
  }

  if (mem_avail < 0) { // not enough memory space allocated for particles: stop the simulation
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

bool c_Solver::ParticlesMover() {

  /*  -------------- */
  /*  Particle mover */
  /*  -------------- */

  // timeTasks.start(TimeTasks::PARTICLES);
  for (int i = 0; i < ns; i++) {  // move each species
    if (cylindrical) {
      // #pragma omp task inout(part[i]) in(grid) target_device(booster)
      mem_avail = part[i].mover_PC_sub_cyl(grid, vct, EMf); // use the Predictor Corrector scheme
    }
    else {
//      // #pragma omp task inout(part[i]) in(grid) target_device(booster)
//      //mem_avail = part[i].mover_PC_sub(grid, vct, EMf); // use the Predictor Corrector scheme
//      if(col->getCase()=="GEMRelativity" || col->getCase()=="Relativistic")
//        mem_avail = part[i].mover_relativistic(grid, vct, EMf);
//      else
//        mem_avail = part[i].mover_PC(grid, vct, EMf); // use the Predictor Corrector scheme
      mem_avail = part[i].mover_PC_rel(grid, vct, EMf); // use the Predictor Corrector scheme
    }
  }
  // timeTasks.end(TimeTasks::PARTICLES);

  if (mem_avail < 0) { // not enough memory space allocated for particles: stop the simulation
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

  if (col->getCase()=="Dipole") {
    for (int i=0; i < ns; i++){
      if (col->getRHOinject(i)>0.0)
           mem_avail = part[i].particle_repopulator(grid,vct,EMf,i);
        Qremoved[i] = part[i].deleteParticlesInsideSphere(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
    }
  }
  else if (col->getCase()=="Coils") {
    //Remove particles from outside the simulation box
    for (int i=0; i < ns; i++){
      //Qremoved[i] = part[i].deleteParticlesOutsideBox(col->getLx());
      /*
      Qremoved[i] = part[i].deleteParticlesOuterFrame(6.0,6.0,6.0);
      if (col->getRHOinject(i) > 0.0)
        mem_avail = part[i].injector_rand_box(grid,vct,EMf);
      */
      Qremoved[i] = part[i].ReturnToCenterCircle();
   }
  }
  else if (col->getCase()=="TwoCoils") {
    //Remove particles from outside the simulation box
    for (int i=0; i < ns; i++){
      //Qremoved[i] = part[i].deleteParticlesOutsideBox(col->getLx());
      // Qremoved[i] = part[i].deleteParticlesOuterFrame(6.0,6.0,6.0);
      /*
      Qremoved[i] =part[i].deleteParticlesOutsideSphere(L_outer, col->getx_center(), col->gety_center(), col->getz_center());
      if (col->getRHOinject(i) > 0.0){
        double x_center_inect = col->getx_center() ;
        double y_center_inect = col->gety_center() + col->getcoilSpacing()/2.0;
        double z_center_inect = col->getz_center() ;
        mem_avail = part[i].injector_rand_box(grid, vct, EMf, x_center_inect, y_center_inect, z_center_inect, L_square );
        x_center_inect = col->getx_center() ;
        y_center_inect = col->gety_center() - col->getcoilSpacing()/2.0;
        z_center_inect = col->getz_center() ;
        mem_avail = part[i].injector_rand_box(grid, vct, EMf, x_center_inect, y_center_inect, z_center_inect, L_square );
      }
      */
      Qremoved[i] = part[i].ReturnToCenterCircle();
    }
  }
  else if (col->getCase()=="CoilsMono") {
    // Remove particles from outside the simulation box
    for (int i=0; i < ns; i++){
      //Qremoved[i] = part[i].deleteParticlesOutsideBox(col->getLx());
      Qremoved[i] = part[i].deleteParticlesOuterFrame(6.0,6.0,6.0);
      if (col->getRHOinject(i) > 0.0)
        mem_avail = part[i].injector_rand_box_mono(grid,vct,EMf);
    }
  }
  else if (col->getCase()=="Relativistic");
  // Do nothing//
  else {

    // REPOPULATOR HERE: commented out for now

    // /* --------------------------------------- */
    // /* Remove particles from depopulation area */
    // /* --------------------------------------- */
    // for (int i=0; i < ns; i++){
    //   mem_avail = part[i].particle_repopulator(grid,vct,EMf,i);
    //   mem_avail = part[i].particle_reflector(grid,vct,EMf,i);
    // }
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
  if (cycle % col->getDiagnosticsOutputCycle() == 0) {
    // Old quantities (unused)
    //TOTenergy = 0.0;
    //TOTmomentum = 0.0;
    //totParticles[0] = totParticles[1] = 0;
    //globalTotParticles[0] = globalTotParticles[1] = 0;

    // Get EM energies
    Eenergy = EMf->getEenergy();
    Benergy = EMf->getBenergy();

    // Get particle energies
    for (int is = 0; is < ns; is++) {
      Ke[is] = part[is].getKe();

      // Old quantites (unused)
      //BulkEnergy[is] = EMf->getBulkEnergy(is);
      //momentum[is] = part[is].getP();
      //TOTmomentum += momentum[is];
      //totParticles[is] = part[is].getNOP();
      //MPI_Allreduce(&totParticles[is], &globalTotParticles[is], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      //Qtot[is] = part[is].getTotalQ();
      //speciesTemp[is] = (double)((double)Ke[is]/(Qtot[is]/qom[is]));
    }

    // Write to file
    // Structure: cycle - Bx By Bz Ex Ey Ez K0 K1 K2...
    if (myrank == 0) {
//      ofstream my_file(cq.c_str(), fstream::app);
//      my_file << cycle << " " << setprecision(15);
//      my_file << Benergy[0] << " " << Benergy[1] << " " << Benergy[2] << " "
//              << Eenergy[0] << " " << Eenergy[1] << " " << Eenergy[2];
//      for (int is=0; is<ns; is++) my_file << " " << Ke[is];
//      my_file << endl;
//      my_file.close();

      FILE* fd = fopen(cq.c_str(),"a");
      fprintf(fd,"%5d %.15e %.15e %.15e %.15e %.15e %.15e",
                cycle,
                Benergy[0], Benergy[1], Benergy[2], Eenergy[0], Eenergy[1], Eenergy[2]);
      for (int is = 0; is < ns; is++) fprintf(fd," %.15e", Ke[is]);
      fprintf(fd,"\n");
      fclose(fd);

      // Next is the diagnostics file (SummaryQuantities), to be made user-defined
      // Unused for now
      //ofstream my_file2(cq2.c_str(),fstream::app);
      //my_file2 << cycle << "\t" <<setprecision(15);
      //for (int i = 0; i < ns; i++)
      //    my_file2 << Ke[i] << "\t";
      //for (int i = 0; i < ns; i++)
      //    my_file2 << BulkEnergy[i] << "\t";
      //for (int i = 0; i < ns; i++)
      //    my_file2 << Qtot[i] << "\t";
      //for (int i = 0; i < ns; i++)
      //        my_file2 << globalTotParticles[i] << "\t";
      //for (int i = 0; i < ns; i++)
      //my_file2 << speciesTemp[i] << "\t";
      //my_file2 << endl;
      //my_file2.close();
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
  }
  
  //if (cycle%(col->getFieldOutputCycle())==0){
  //  for (int is = 0; is < ns; is++) {
  //    part[is].Add_vDist3D();
  //    part[is].Write_vDist3D(SaveDirName);
  //  }
  //}
}

void c_Solver::WriteOutput(int cycle) {

  if (col->getWriteMethod() == "h5hut") {

    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */

    if (cycle%(col->getFieldOutputCycle())==0 || (!(col->getSolInit()) && cycle==first_cycle)) {
      // set to zero the densities
      EMf->setZeroDensitiesOutput();
      // interpolate Particles to Grid(Nodes)
      for (int i = 0; i < ns; i++)
        part[i].interpP2GOutput(EMf, grid, vct); 
      // Write fields
      WriteFieldsH5hut(ns, grid, EMf,  col, vct, cycle);
    }
    //if (cycle%(col->getParticlesOutputCycle())==0 &&
    //    cycle!=col->getLast_cycle() && cycle!=0)      WritePartclH5hut(ns, grid, part, col, vct, cycle);
    if (cycle%(col->getParticlesOutputCycle())==0)    WritePartclH5hut(ns, grid, part, col, vct, cycle);

  }
  else
  {

    // OUTPUT to large file, called proc**
    if (cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("Eall + Ball + rhos + Jsall + pressure", cycle);
      // Pressure tensor is available
      hdf5_agent.close();
    }
    if (cycle % (col->getParticlesOutputCycle()) == 0 && col->getParticlesOutputCycle() != 1) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("position + velocity + q + ID", cycle, 1);
      hdf5_agent.close();
    }
  }
    // write the virtual satellite traces
/*
    bool binary_satellites = false;
    if(binary_satellites){
    float time_counter = cycle;
    float trace_counter = 0.0;
    float outta = 0.0;
    ofstream my_file(cqsat.c_str(),ofstream::app | ofstream::binary );
                     for (int isat=0; isat < nsatx; isat++)
                         for (int jsat=0; jsat < nsaty; jsat++)
                             for (int ksat=0; ksat < nsatz; ksat++){
                             int index1 = 1+isat*nx0/nsatx+nx0/nsatx/2;
                             int index2 = 1+jsat*ny0/nsaty+ny0/nsaty/2;
                              int index3 = 1+ksat*nz0/nsatz+nz0/nsatz/2;
                                   trace_counter = trace_counter +1.0;
                                 my_file.write((char*)&time_counter,sizeof(float));
                                 my_file.write((char*)&trace_counter,sizeof(float));
                                 //my_file << time_counter << "\t" <<trace_counter << "\t";

                                 outta = EMf->getBx(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";
                                 outta = EMf->getBy(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";
                                 outta = EMf->getBz(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";

                                 outta = EMf->getEx(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";
                                 outta = EMf->getEy(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";
                                 outta = EMf->getEz(index1,index2,index3);
                                 my_file.write((char*)&outta,sizeof(float));
                                 //my_file << outta << "\t";

                                 for (int is=0;is<2; is++){
                                     outta = EMf->getRHOns(index1,index2,index3,is);
                                     if(ns>3)outta += EMf->getRHOns(index1,index2,index3,is+2);
                                     my_file.write((char*)&outta,sizeof(float));
                                     //my_file << outta << "\t";
                                     outta = EMf->getJxs(index1,index2,index3,is);
                                     if(ns>3)outta += EMf->getJxs(index1,index2,index3,is+2);
                                     my_file.write((char*)&outta,sizeof(float));
                                     //my_file << outta << "\t";
                                     outta = EMf->getJys(index1,index2,index3,is);
                                     if(ns>3)outta += EMf->getJys(index1,index2,index3,is+2);
                                     my_file.write((char*)&outta,sizeof(float));
                                     //my_file << outta << "\t";
                                     outta = EMf->getJzs(index1,index2,index3,is);
                                     if(ns>3)outta += EMf->getJzs(index1,index2,index3,is+2);
                                     my_file.write((char*)&outta,sizeof(float));
                                     //my_file << outta << "\t"
                                 }
                                 //my_file  << endl;
                             }
                     my_file.close();
    }
    else {
    if (ns > 2) {
      ofstream my_file(cqsat.c_str(), fstream::app);
      for (int isat = 0; isat < nsatx; isat++)
        for (int jsat = 0; jsat < nsaty; jsat++)
          for (int ksat = 0; ksat < nsatz; ksat++) {
             int index1 = 1+isat*nx0/nsatx+nx0/nsatx/2;
             int index2 = 1+jsat*ny0/nsaty+ny0/nsaty/2;
              int index3 = 1+ksat*nz0/nsatz+nz0/nsatz/2;
            my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
            my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 0) + EMf->getJxs(index1, index2, index3, 2) << "\t" << EMf->getJys(index1, index2, index3, 0) + EMf->getJys(index1, index2, index3, 2) << "\t" << EMf->getJzs(index1, index2, index3, 0) + EMf->getJzs(index1, index2, index3, 2) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 1) + EMf->getJxs(index1, index2, index3, 3) << "\t" << EMf->getJys(index1, index2, index3, 1) + EMf->getJys(index1, index2, index3, 3) << "\t" << EMf->getJzs(index1, index2, index3, 1) + EMf->getJzs(index1, index2, index3, 3) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 0) + EMf->getRHOns(index1, index2, index3, 2) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 1) + EMf->getRHOns(index1, index2, index3, 3) << "\t";
          }
          my_file << endl;
          my_file.close();
        }
      if (ns == 2) {
              ofstream my_file(cqsat.c_str(), fstream::app);
              for (int isat = 0; isat < nsatx; isat++)
                for (int jsat = 0; jsat < nsaty; jsat++)
                  for (int ksat = 0; ksat < nsatz; ksat++) {
                       int index1 = 1+isat*nx0/nsatx+nx0/nsatx/2;
                       int index2 = 1+jsat*ny0/nsaty+ny0/nsaty/2;
                        int index3 = 1+ksat*nz0/nsatz+nz0/nsatz/2;
                    my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
                    my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
                    my_file << EMf->getJxs(index1, index2, index3, 0)  << "\t" << EMf->getJys(index1, index2, index3, 0)  << "\t" << EMf->getJzs(index1, index2, index3, 0)  << "\t";
                    my_file << EMf->getJxs(index1, index2, index3, 1)  << "\t" << EMf->getJys(index1, index2, index3, 1)  << "\t" << EMf->getJzs(index1, index2, index3, 1)  << "\t";
                    my_file << EMf->getRHOns(index1, index2, index3, 0)  << "\t";
                    my_file << EMf->getRHOns(index1, index2, index3, 1)  << "\t";
                  }

      my_file << endl;
      my_file.close();
    }
    }
*/
}

void c_Solver::Finalize() {
  if (mem_avail == 0) {          // write the restart only if the simulation finished succesfully
    if (col->getWriteMethod() != "h5hut") {
      writeRESTART(RestartDirName, myrank, (col->getNcycles() + first_cycle) - 1, ns, mpi, vct, col, grid, EMf, part, 0);
    }
  }

  // stop profiling
  clocks->stop(5);
  // stop profiling
  my_clock->stopTiming();

  //if(myrank == 0) cout << "Deallocating"<<endl;


  // deallocate
  delete[]Ke;
  delete[]BulkEnergy;
  delete[]momentum;
  delete[] Qtot;

  delete[] Qremoved;
  //delete[] VelocityDist;

  delete[] totParticles;
  delete[] globalTotParticles;
  delete[] speciesTemp;
  delete[] qom;


  //if(myrank == 0) cout << "Closing MPI"<<endl;

  // close MPI
  mpi->finalize_mpi();

  //if(myrank == 0) cout << "Closed MPI"<<endl;
}
