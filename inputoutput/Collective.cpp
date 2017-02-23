
#include "Collective.h"
#include "Alloc.h"   // mlmd: added for some mlmd vectors


/*! Read the input file from text file and put the data in a collective wrapper: if it's a restart read from input file basic sim data and load particles and EM field from restart file */
void Collective::ReadInput(string inputfile) {
  using namespace std;
  int test_verbose;
  // Loading the input file 
  ConfigFile config(inputfile);
  // the following variables are ALWAYS taken from inputfile, even if restarting 
  try {
  {

#ifdef BATSRUS
    if(RESTART1)
    {
      cout<<" The fluid interface can not handle RESTART yet, aborting!\n"<<flush;
      abort();
    }
#endif

    dt = config.read < double >("dt");
    ncycles = config.read < int >("ncycles");
    th = config.read < double >("th");
    config.readInto(Smooth, "Smooth");
    SaveDirName = config.read < string > ("SaveDirName");
    RestartDirName = config.read < string > ("RestartDirName");
    ns = config.read < int >("ns");
    NpMaxNpRatio = config.read < double >("NpMaxNpRatio");
    // GEM Challenge 
    B0x = config.read <double>("B0x");
    B0y = config.read <double>("B0y");
    B0z = config.read <double>("B0z");

    // Earth parameters
    B1x = 0.0;
    B1y = 0.0;
    B1z = 0.0;
    B1x = config.read <double>("B1x");
    B1y = config.read <double>("B1y");
    B1z = config.read <double>("B1z");

    delta = config.read < double >("delta");

    Case              = config.read<string>("Case");
    FieldsInit        = config.read<string>("FieldsInit");
    PartInit          = config.read<string>("PartInit");
    wmethod           = config.read<string>("WriteMethod");
    SimName           = config.read<string>("SimulationName");
    PoissonCorrection = config.read<string>("PoissonCorrection");

    rhoINIT = new double[ns];
    array_double rhoINIT0 = config.read < array_double > ("rhoINIT");
    rhoINIT[0] = rhoINIT0.a;
    if (ns > 1)
      rhoINIT[1] = rhoINIT0.b;
    if (ns > 2)
      rhoINIT[2] = rhoINIT0.c;
    if (ns > 3)
      rhoINIT[3] = rhoINIT0.d;
    if (ns > 4)
      rhoINIT[4] = rhoINIT0.e;
    if (ns > 5)
      rhoINIT[5] = rhoINIT0.f;

    rhoINJECT = new double[ns];
    array_double rhoINJECT0 = config.read<array_double>( "rhoINJECT" );
    rhoINJECT[0]=rhoINJECT0.a;
    if (ns > 1)
      rhoINJECT[1]=rhoINJECT0.b;
    if (ns > 2)
      rhoINJECT[2]=rhoINJECT0.c;
    if (ns > 3)
      rhoINJECT[3]=rhoINJECT0.d;
    if (ns > 4)
      rhoINJECT[4]=rhoINJECT0.e;
    if (ns > 5)
      rhoINJECT[5]=rhoINJECT0.f;

    // take the tolerance of the solvers
    CGtol = config.read < double >("CGtol");
    GMREStol = config.read < double >("GMREStol");
    NiterMover = config.read < int >("NiterMover");
    // take the injection of the particless
    Vinj = config.read < double >("Vinj");

    // take the output cycles
    FieldOutputCycle = config.read < int >("FieldOutputCycle");
    ParticlesOutputCycle = config.read < int >("ParticlesOutputCycle");
    RestartOutputCycle = config.read < int >("RestartOutputCycle");
    DiagnosticsOutputCycle = config.read < int >("DiagnosticsOutputCycle", FieldOutputCycle);
    
    // MPI topology and periodicity
    XLEN      = config.read < int > ("XLEN",1);
    YLEN      = config.read < int > ("YLEN",1);
    ZLEN      = config.read < int > ("ZLEN",1);
    PERIODICX = config.read < bool >("PERIODICX");
    PERIODICY = config.read < bool >("PERIODICY");
    PERIODICZ = config.read < bool >("PERIODICZ");

    // MLMD reads; always read from inputfile
    
    Ngrids = config.read < int >("Ngrids");
    array_int gridLevel0 = config.read < array_int > ("gridLevel");
    array_int RF0 = config.read < array_int > ("RF");
    array_int parentGrid0 = config.read < array_int > ("parentGrid");
    
    array_int nxc_mlmd0 = config.read < array_int > ("nxc_mlmd");
    array_int nyc_mlmd0 = config.read < array_int > ("nyc_mlmd");
    array_int nzc_mlmd0 = config.read < array_int > ("nzc_mlmd");
    
    array_double Ox0 = config.read < array_double > ("Ox_SW");
    array_double Oy0 = config.read < array_double > ("Oy_SW");
    array_double Oz0 = config.read < array_double > ("Oz_SW");

    gridLevel = new int[Ngrids];
    RF = new int[Ngrids];
    parentGrid = new int[Ngrids];

    nxc_mlmd = new int[Ngrids];
    nyc_mlmd = new int[Ngrids];
    nzc_mlmd = new int[Ngrids];

    Ox_SW = new double[Ngrids];
    Oy_SW = new double[Ngrids];
    Oz_SW = new double[Ngrids];

    gridLevel[0] = gridLevel0.a;
    RF[0] = RF0.a;
    parentGrid[0] = parentGrid0.a;
        
    nxc_mlmd[0] = nxc_mlmd0.a;
    nyc_mlmd[0] = nyc_mlmd0.a;
    nzc_mlmd[0] = nzc_mlmd0.a;
    
    Ox_SW[0] = Ox0.a;
    Oy_SW[0] = Oy0.a;
    Oz_SW[0] = Oz0.a;
    
    if (Ngrids >1) {
      gridLevel[1] = gridLevel0.b;
      RF[1]= RF0.b;
      parentGrid[1]= parentGrid0.b;
      
      nxc_mlmd[1] =nxc_mlmd0.b;
      nyc_mlmd[1] = nyc_mlmd0.b;
      nzc_mlmd[1] = nzc_mlmd0.b;
      
      Ox_SW[1]= Ox0.b;
      Oy_SW[1]= Oy0.b;
      Oz_SW[1]= Oz0.b;
    }
    
    if (Ngrids >2) {
      gridLevel[2] = gridLevel0.c;
      RF[2]= RF0.c;
      parentGrid[2]= parentGrid0.c;
      
      nxc_mlmd[2] =nxc_mlmd0.c;
      nyc_mlmd[2] = nyc_mlmd0.c;
      nzc_mlmd[2] = nzc_mlmd0.c;
      
      Ox_SW[2]= Ox0.c;
      Oy_SW[2]= Oy0.c;
      Oz_SW[2]= Oz0.c;
    }
    
    if (Ngrids >3) {
      gridLevel[3] = gridLevel0.d;
      RF[3]= RF0.d;
      parentGrid[3]= parentGrid0.d;
      
      nxc_mlmd[3] =nxc_mlmd0.d;
      nyc_mlmd[3] = nyc_mlmd0.d;
      nzc_mlmd[3] = nzc_mlmd0.d;
      
      Ox_SW[3]= Ox0.d;
      Oy_SW[3]= Oy0.d;
      Oz_SW[3]= Oz0.d;
    }
    
    if (Ngrids >4) {
      gridLevel[4] = gridLevel0.e;
      RF[4]= RF0.e;
      parentGrid[4]= parentGrid0.e;
      
      nxc_mlmd[4] =nxc_mlmd0.e;
      nyc_mlmd[4] = nyc_mlmd0.e;
      nzc_mlmd[4] = nzc_mlmd0.e;
      
      Ox_SW[4]= Ox0.e;
      Oy_SW[4]= Oy0.e;
      Oz_SW[4]= Oz0.e;
    }
    
    TopologyType = config.read < int > ("TopologyType");
    
    /* whether to perform mlmd operations */
    MLMD_BC = config.read < int > ("MLMD_BC");
    MLMD_PROJECTION = config.read < int > ("MLMD_PROJECTION");
    MLMD_ParticleREPOPULATION = config.read < int > ("MLMD_ParticleREPOPULATION");

    //cout << "MLMD_BC: " << MLMD_BC <<endl;
    
    // end MLMD reads
  }

  SolInit = false;
  if (RESTART1) {               // you are restarting
    cout << " WARNING: Performing a restart with the classic files. " << endl;

    RestartDirName = config.read < string > ("RestartDirName");
    ReadRestart(RestartDirName);
  }
  else if (SOLINIT1) {
#ifdef USEH5HUT
    SolInit  = true;
    initfile = FieldsInit;

    /* ------------------------------------ */
    /* Verify that the file name is correct */
    /* ------------------------------------ */
    stringstream ss;
    string       ff;
    unsigned first = initfile.find_last_of("-") + 1;
    unsigned last  = initfile.find_last_of("_");
    ss << initfile.substr(first, last-first);
    ss >> ff;

    if (ff=="Fields"){

      /* ---------------------------------- */
      /* Get the cycle number from the name */
      /* ---------------------------------- */
      first = initfile.find_last_of("_") + 1;
      last  = initfile.find_last_of(".h5");
      ss.clear();
      ss.str(string());
      ss << initfile.substr(first, last-first);
      ss >> last_cycle;

      /* --------------------------------- */
      /* Extract the base name of the file */
      /* --------------------------------- */
      last  = initfile.find_last_of("-");
      ss.clear();
      ss.str(string());
      ss << initfile.substr(0, last);
      ss >> ff;

      initfile = ff;

    }
    else {
      cout << " ERROR: The name of the initial file is incorrect, please verify that you are " << endl;
      cout << "        using an iPic3D file with the correct name: FileName-Fields_000000.h5" << endl;
    }

#else
    cout << " ERROR: You want to restart from an HDF5 initial solution using the H5hut library. " << endl;
    cout << "        However, the code was compiled without the H5hut options." << endl;
    abort();
#endif
  }

  if(!RESTART1) {

    if (!SOLINIT1) last_cycle = -1;
    restart_status = 0;
    c = config.read < double >("c");

#ifdef BATSRUS
    // set grid size and resolution based on the initial file from fluid code
    Lx =  getFluidLx();
    Ly =  getFluidLy();
    Lz =  getFluidLz();
    nxc = getFluidNxc();
    nyc = getFluidNyc();
    nzc = getFluidNzc();
#else
    Lx = config.read < double >("Lx");  // I need this read for the MLMD also
    Ly = config.read < double >("Ly");
    Lz = config.read < double >("Lz");
    nxc = config.read < int >("nxc");
    nyc = config.read < int >("nyc");
    nzc = config.read < int >("nzc");
#endif

    

    x_center = config.read < double >("x_center");
    y_center = config.read < double >("y_center");
    z_center = config.read < double >("z_center");
    L_square = config.read < double >("L_square");

    npcelx = new int[ns];
    npcely = new int[ns];
    npcelz = new int[ns];
    qom = new double[ns];
    uth = new double[ns];
    vth = new double[ns];
    wth = new double[ns];
    u0 = new double[ns];
    v0 = new double[ns];
    w0 = new double[ns];

    array_int npcelx0 = config.read < array_int > ("npcelx");
    array_int npcely0 = config.read < array_int > ("npcely");
    array_int npcelz0 = config.read < array_int > ("npcelz");
    array_double qom0 = config.read < array_double > ("qom");
    array_double uth0 = config.read < array_double > ("uth");
    array_double vth0 = config.read < array_double > ("vth");
    array_double wth0 = config.read < array_double > ("wth");
    array_double u00 = config.read < array_double > ("u0");
    array_double v00 = config.read < array_double > ("v0");
    array_double w00 = config.read < array_double > ("w0");

    npcelx[0] = npcelx0.a;
    npcely[0] = npcely0.a;
    npcelz[0] = npcelz0.a;
    qom[0] = qom0.a;
    uth[0] = uth0.a;
    vth[0] = vth0.a;
    wth[0] = wth0.a;
    u0[0] = u00.a;
    v0[0] = v00.a;
    w0[0] = w00.a;

    if (ns > 1) {
      npcelx[1] = npcelx0.b;
      npcely[1] = npcely0.b;
      npcelz[1] = npcelz0.b;
      qom[1] = qom0.b;
      uth[1] = uth0.b;
      vth[1] = vth0.b;
      wth[1] = wth0.b;
      u0[1] = u00.b;
      v0[1] = v00.b;
      w0[1] = w00.b;
    }
    if (ns > 2) {
      npcelx[2] = npcelx0.c;
      npcely[2] = npcely0.c;
      npcelz[2] = npcelz0.c;
      qom[2] = qom0.c;
      uth[2] = uth0.c;
      vth[2] = vth0.c;
      wth[2] = wth0.c;
      u0[2] = u00.c;
      v0[2] = v00.c;
      w0[2] = w00.c;
    }
    if (ns > 3) {
      npcelx[3] = npcelx0.d;
      npcely[3] = npcely0.d;
      npcelz[3] = npcelz0.d;
      qom[3] = qom0.d;
      uth[3] = uth0.d;
      vth[3] = vth0.d;
      wth[3] = wth0.d;
      u0[3] = u00.d;
      v0[3] = v00.d;
      w0[3] = w00.d;
    }
    if (ns > 4) {
      npcelx[4] = npcelx0.e;
      npcely[4] = npcely0.e;
      npcelz[4] = npcelz0.e;
      qom[4] = qom0.e;
      uth[4] = uth0.e;
      vth[4] = vth0.e;
      wth[4] = wth0.e;
      u0[4] = u00.e;
      v0[4] = v00.e;
      w0[4] = w00.e;
    }
    if (ns > 5) {
      npcelx[5] = npcelx0.f;
      npcely[5] = npcely0.f;
      npcelz[5] = npcelz0.f;
      qom[5] = qom0.f;
      uth[5] = uth0.f;
      vth[5] = vth0.f;
      wth[5] = wth0.f;
      u0[5] = u00.f;
      v0[5] = v00.f;
      w0[1] = w00.f;
    }

    verbose = config.read < bool > ("verbose");

    // PHI Electrostatic Potential 
    bcPHIfaceXright = config.read < int >("bcPHIfaceXright");
    bcPHIfaceXleft  = config.read < int >("bcPHIfaceXleft");
    bcPHIfaceYright = config.read < int >("bcPHIfaceYright");
    bcPHIfaceYleft  = config.read < int >("bcPHIfaceYleft");
    bcPHIfaceZright = config.read < int >("bcPHIfaceZright");
    bcPHIfaceZleft  = config.read < int >("bcPHIfaceZleft");

    // EM field boundary condition
    bcEMfaceXright = config.read < int >("bcEMfaceXright");
    bcEMfaceXleft  = config.read < int >("bcEMfaceXleft");
    bcEMfaceYright = config.read < int >("bcEMfaceYright");
    bcEMfaceYleft  = config.read < int >("bcEMfaceYleft");
    bcEMfaceZright = config.read < int >("bcEMfaceZright");
    bcEMfaceZleft  = config.read < int >("bcEMfaceZleft");

    /*  ---------------------------------------------------------- */
    /*  Electric and Magnetic field boundary conditions for BCface */
    /*  ---------------------------------------------------------- */
    // if bcEM* is 0: perfect conductor, if bcEM* is not 0: perfect mirror
    // perfect conductor: normal = free, perpendicular = 0   
    // perfect mirror   : normal = 0,    perpendicular = free
    /*  ---------------------------------------------------------- */

    /* X component in faces Xright, Xleft, Yright, Yleft, Zright and Zleft (0, 1, 2, 3, 4, 5) */
    bcEx[0] = bcEMfaceXright == 0 ? 2 : 1;   bcBx[0] = bcEMfaceXright == 0 ? 1 : 2;
    bcEx[1] = bcEMfaceXleft  == 0 ? 2 : 1;   bcBx[1] = bcEMfaceXleft  == 0 ? 1 : 2;
    bcEx[2] = bcEMfaceYright == 0 ? 1 : 2;   bcBx[2] = bcEMfaceYright == 0 ? 2 : 1;
    bcEx[3] = bcEMfaceYleft  == 0 ? 1 : 2;   bcBx[3] = bcEMfaceYleft  == 0 ? 2 : 1;
    bcEx[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBx[4] = bcEMfaceZright == 0 ? 2 : 1;
    bcEx[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBx[5] = bcEMfaceZleft  == 0 ? 2 : 1;
    /* Y component */
    bcEy[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBy[0] = bcEMfaceXright == 0 ? 2 : 1;
    bcEy[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBy[1] = bcEMfaceXleft  == 0 ? 2 : 1;
    bcEy[2] = bcEMfaceYright == 0 ? 2 : 1;   bcBy[2] = bcEMfaceYright == 0 ? 1 : 2;
    bcEy[3] = bcEMfaceYleft  == 0 ? 2 : 1;   bcBy[3] = bcEMfaceYleft  == 0 ? 1 : 2;
    bcEy[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBy[4] = bcEMfaceZright == 0 ? 2 : 1;
    bcEy[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBy[5] = bcEMfaceZleft  == 0 ? 2 : 1;
    /* Z component */
    bcEz[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBz[0] = bcEMfaceXright == 0 ? 2 : 1;
    bcEz[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBz[1] = bcEMfaceXleft  == 0 ? 2 : 1;
    bcEz[2] = bcEMfaceYright == 0 ? 1 : 1;   bcBz[2] = bcEMfaceYright == 0 ? 2 : 1;
    bcEz[3] = bcEMfaceYleft  == 0 ? 1 : 1;   bcBz[3] = bcEMfaceYleft  == 0 ? 2 : 1;
    bcEz[4] = bcEMfaceZright == 0 ? 2 : 1;   bcBz[4] = bcEMfaceZright == 0 ? 1 : 2;
    bcEz[5] = bcEMfaceZleft  == 0 ? 2 : 1;   bcBz[5] = bcEMfaceZleft  == 0 ? 1 : 2;

    // Particles Boundary condition
    bcPfaceXright = config.read < int >("bcPfaceXright");
    bcPfaceXleft  = config.read < int >("bcPfaceXleft");
    bcPfaceYright = config.read < int >("bcPfaceYright");
    bcPfaceYleft  = config.read < int >("bcPfaceYleft");
    bcPfaceZright = config.read < int >("bcPfaceZright");
    bcPfaceZleft  = config.read < int >("bcPfaceZleft");




  }
  TrackParticleID = new bool[ns];
  array_bool TrackParticleID0 = config.read < array_bool > ("TrackParticleID");
  TrackParticleID[0] = TrackParticleID0.a;
  if (ns > 1)
    TrackParticleID[1] = TrackParticleID0.b;
  if (ns > 2)
    TrackParticleID[2] = TrackParticleID0.c;
  if (ns > 3)
    TrackParticleID[3] = TrackParticleID0.d;
  if (ns > 4)
    TrackParticleID[4] = TrackParticleID0.e;
  if (ns > 5)
    TrackParticleID[5] = TrackParticleID0.f;

  } catch (ConfigFile::key_not_found k) {
    std::cout << " ERROR: Key not found = " << k.key.c_str() << std::endl;
    abort();
  }

  /*! this will be eventually in inputfile */
  MLMDVerbose = true; 

}
/*! Read the collective information from the RESTART file in HDF5 format There are three restart status: restart_status = 0 ---> new inputfile restart_status = 1 ---> RESTART and restart and result directories does not coincide restart_status = 2 ---> RESTART and restart and result directories coincide */



int Collective::ReadRestart(string inputfile) {

  if (Ngrids > 1) {
    cout<<" The MLMD version cannot handle RESTART yet, aborting...\n"<<flush;
    abort();
  }

  restart_status = 1;
  // hdf stuff 
  hid_t file_id;
  hid_t dataset_id;
  herr_t status;
  // Open the setting file for the restart.
  file_id = H5Fopen((inputfile + "/settings.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << inputfile << endl;
    return -1;
  }

  // read c
  dataset_id = H5Dopen2(file_id, "/collective/c", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c);
  status = H5Dclose(dataset_id);

  // read Lx 
  dataset_id = H5Dopen2(file_id, "/collective/Lx", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
  status = H5Dclose(dataset_id);
  // read Ly 
  dataset_id = H5Dopen2(file_id, "/collective/Ly", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
  status = H5Dclose(dataset_id);
  // read Lz 
  dataset_id = H5Dopen2(file_id, "/collective/Lz", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lz);
  status = H5Dclose(dataset_id);
  // read x_center
  dataset_id = H5Dopen2(file_id, "/collective/x_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x_center);
  status = H5Dclose(dataset_id);
  // read y_center
  dataset_id = H5Dopen2(file_id, "/collective/y_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y_center);
  status = H5Dclose(dataset_id);
  // read z_center
  dataset_id = H5Dopen2(file_id, "/collective/z_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z_center);
  status = H5Dclose(dataset_id);
  // read L_square
  dataset_id = H5Dopen2(file_id, "/collective/L_square", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L_square);
  status = H5Dclose(dataset_id);
  // read nxc
  dataset_id = H5Dopen2(file_id, "/collective/Nxc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nyc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nzc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzc);
  status = H5Dclose(dataset_id);
  // read ns
  dataset_id = H5Dopen2(file_id, "/collective/Ns", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ns);
  status = H5Dclose(dataset_id);


  /*! Boundary condition information */
  // read EMfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXleft);
  status = H5Dclose(dataset_id);
  // read EMfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXright);
  status = H5Dclose(dataset_id);
  // read EMfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYleft);
  status = H5Dclose(dataset_id);
  // read EMfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYright);
  status = H5Dclose(dataset_id);
  // read EMfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZleft);
  status = H5Dclose(dataset_id);
  // read EMfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZright);
  status = H5Dclose(dataset_id);

  // read PHIfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXright);
  status = H5Dclose(dataset_id);
  // read PHIfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYright);
  status = H5Dclose(dataset_id);
  // read PHIfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZright);
  status = H5Dclose(dataset_id);

  // read PfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXleft);
  status = H5Dclose(dataset_id);
  // read PfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXright);
  status = H5Dclose(dataset_id);
  // read PfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYleft);
  status = H5Dclose(dataset_id);
  // read PfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYright);
  status = H5Dclose(dataset_id);
  // read PfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZleft);
  status = H5Dclose(dataset_id);
  // read PfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZright);
  status = H5Dclose(dataset_id);
  // allocate fields depending on species
  npcelx = new int[ns];
  npcely = new int[ns];
  npcelz = new int[ns];
  qom = new double[ns];
  uth = new double[ns];
  vth = new double[ns];
  wth = new double[ns];
  u0 = new double[ns];
  v0 = new double[ns];
  w0 = new double[ns];
  // read data from species0, species 1, species2,...
  string *name_species = new string[ns];
  stringstream *ss = new stringstream[ns];

  for (int i = 0; i < ns; i++) {
    ss[i] << i;
    name_species[i] = "/collective/species_" + ss[i].str() + "/";
  }
  // npcelx for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelx").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  // npcely for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcely").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcely[i]);
    status = H5Dclose(dataset_id);
  }
  // npcelz for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelz").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelz[i]);
    status = H5Dclose(dataset_id);
  }
  // qom for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "qom").c_str(), H5P_DEFAULT);  // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &qom[i]);
    status = H5Dclose(dataset_id);
  }
  /*! not needed for restart * */
  for (int i = 0; i < ns; i++)
    uth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    vth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    wth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    u0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    v0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    w0[i] = 0.0;
  // verbose on
  verbose = 1;


  // if RestartDirName == SaveDirName overwrite dt,Th,Smooth (append to old hdf files)
  if (RestartDirName == SaveDirName) {
    restart_status = 2;
    // read dt
    dataset_id = H5Dopen2(file_id, "/collective/Dt", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dt);
    status = H5Dclose(dataset_id);
    // read th 
    dataset_id = H5Dopen2(file_id, "/collective/Th", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &th);
    status = H5Dclose(dataset_id);
    // read Smooth
    dataset_id = H5Dopen2(file_id, "/collective/Smooth", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Smooth);
    status = H5Dclose(dataset_id);
  }

  status = H5Fclose(file_id);


  // read last cycle (not from settings, but from restart0.hdf)

  file_id = H5Fopen((inputfile + "/restart0.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << inputfile << endl;
    return -1;
  }

  dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  // deallocate
  delete[]name_species;
  delete[]ss;

  return (0);
}
/*! constructor */
Collective::Collective(int argc, char **argv) {
  SOLINIT1 = false;
  RESTART1 = false;
  initfile  = "dummy-Fields_00000.h5";
  if (argc < 2) {
    inputfile = "inputfile";
  }
  else if (argc < 3) {
    inputfile = argv[1];
  }
  else if (argc < 4) {

    if (strcmp(argv[1], "restart") == 0) {
      inputfile = argv[2];
      RESTART1 = true;
    }
    else if (strcmp(argv[1], "solinit") == 0) {
      inputfile = argv[2];
      SOLINIT1 = true;
    }

  }
  else {
    cout << "Error: syntax error in mpirun arguments. Did you mean to 'restart' or 'soliniti' ?" << endl;
    return;
  }

  ReadInput(inputfile);
  /*! fourpi = 4 greek pi */
  fourpi = 16.0 * atan(1.0);
  /*! dx = space step - X direction */
  dx = Lx / (double) nxc;
  /*! dy = space step - Y direction */
  dy = Ly / (double) nyc;
  /*! dz = space step - Z direction */
  dz = Lz / (double) nzc;

  /*! MLMD:  dx_mlmd, dy_mlmd, dz_mlmd : resolution at grid level */

  dx_mlmd = new double[Ngrids];
  dy_mlmd = new double[Ngrids];
  dz_mlmd = new double[Ngrids];

  Lx_mlmd = new double[Ngrids];
  Ly_mlmd = new double[Ngrids];
  Lz_mlmd = new double[Ngrids];

  Ox_P = new double[Ngrids];
  Oy_P = new double[Ngrids];
  Oz_P = new double[Ngrids];
  
  Lx_mlmd[0] = Lx;
  Ly_mlmd[0] = Ly;
  Lz_mlmd[0] = Lz;

  for (int ng=1; ng < Ngrids; ng++) { // NB: RF is given with respect to the PARENT grid
    Lx_mlmd[ng]= Lx_mlmd[parentGrid[ng]]/RF[ng];
    Ly_mlmd[ng]= Ly_mlmd[parentGrid[ng]]/RF[ng];
    Lz_mlmd[ng]= Lz_mlmd[parentGrid[ng]]/RF[ng];
  }

  for (int ng=0; ng < Ngrids; ng++) { // NB: RF is given with respect to the PARENT grid 
    dx_mlmd[ng]= Lx_mlmd[ng]/ (double) nxc_mlmd[ng];
    dy_mlmd[ng]= Ly_mlmd[ng]/ (double) nyc_mlmd[ng];
    dz_mlmd[ng]= Lz_mlmd[ng]/ (double) nzc_mlmd[ng];
  }

  /*! Ox_P, Oy_P, Oz_P are in terms of the PARENT grid */
  Ox_P[0]= Ox_SW[0];
  Oy_P[0]= Oy_SW[0];
  Oz_P[0]= Oz_SW[0];
  for (int ng=1; ng< Ngrids; ng++) {
    Ox_P[ng]= Ox_SW[ng]- Ox_SW[parentGrid[ng]];
    Oy_P[ng]= Oy_SW[ng]- Oy_SW[parentGrid[ng]];
    Oz_P[ng]= Oz_SW[ng]- Oz_SW[parentGrid[ng]];
  }
  
  /* MLMD:  dx_mlmd, dy_mlmd, dz_mlmd : resolution at grid level */
  

  /*! npcel = number of particles per cell */
  npcel = new int[ns];
  /*! np = number of particles of different species */
  np = new long[ns];
  /*! npMax = maximum number of particles of different species */
  npMax = new long[ns];

  // MLMD: maybe np should be had per grid as well
  for (int i = 0; i < ns; i++) {
    npcel[i] = npcelx[i] * npcely[i] * npcelz[i];
    np[i] = npcel[i] * (nxc/XLEN) * (nyc/YLEN) * (nzc/ZLEN);
    npMax[i] = (long) (NpMaxNpRatio * np[i]);
  }

  /*! building the list of the children grids */
  childrenNum = new int[Ngrids];
  childrenGrids = newArr2(int, Ngrids, Ngrids);

  /*! this to trigger segFault in case of errors */
  for (int i=0; i< Ngrids; i++)
    for (int j=0; j< Ngrids; j++)
      childrenGrids[i][j]= -1;
  
  parentGrid[0] = -1;  
  /*! end this to trigger segFault in case of errors */

  for (int ng=0; ng< Ngrids; ng++){
    childrenNum[ng]=0;
  }
  for (int ng=1; ng< Ngrids; ng++){ // cycle on the children; we skip Grid 0
    int parent= parentGrid[ng];
    childrenGrids[parent][childrenNum[parent]]= ng;
    childrenNum[parent]++; 
    }

  /*! debug */
  int SpokePerson;
  MPI_Comm_rank(MPI_COMM_WORLD, &SpokePerson);
  MPI_Barrier(MPI_COMM_WORLD);
  if (0 && SpokePerson==0){
    for (int ng=0; ng<Ngrids; ng++) {
      cout << "Grid " << ng << " has " << childrenNum[ng] << " child(ren)" << endl;
      if (childrenNum[ng]>0){
	for (int c=0; c<childrenNum[ng]; c++)
	  cout << "Child grid " << childrenGrids[ng][c] << endl;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /*! end debug */

  /*! end building the list of the children grids */

  /*! a first sanity check on MLMD inputs, called at the end of the constructor; may abort internally */
  checkMLMDinputs();

}

/*! destructor */
Collective::~Collective() {
  delete[]np;
  delete[]npcel;
  delete[]npcelx;
  delete[]npcely;
  delete[]npcelz;
  delete[]npMax;
  delete[]qom;

  delete[]uth;
  delete[]vth;
  delete[]wth;

  delete[]u0;
  delete[]v0;
  delete[]w0;

  delete[]TrackParticleID;

  delete[]rhoINIT;
  delete[]rhoINJECT;

  // MLMD variables
  delete[]gridLevel;
  delete[]RF;
  delete[]parentGrid;
  delete[]nxc_mlmd;
  delete[]nyc_mlmd;
  delete[]nzc_mlmd;
  delete[]Ox_SW;
  delete[]Oy_SW;
  delete[]Oz_SW;
  delete[]Ox_P;
  delete[]Oy_P;
  delete[]Oz_P;
  delete[]dx_mlmd;
  delete[]dy_mlmd;
  delete[]dz_mlmd;
  // MLMD variables
}
/*! Print Simulation Parameters */
void Collective::Print() {
  cout << endl;
  cout << "Simulation Parameters" << endl;
  cout << "---------------------" << endl;
  cout << "Number of species    = " << ns << endl;
  for (int i = 0; i < ns; i++)
    cout << "Number of particles per proc of species " << i << " = " << np[i] << "\t (MAX = " << npMax[i] << ")" << "  QOM = " << qom[i] << endl;
  cout << "x-Length                 = " << Lx << endl;
  cout << "y-Length                 = " << Ly << endl;
  cout << "z-Length                 = " << Lz << endl;
  cout << "Number of cells (x)      = " << nxc << endl;
  cout << "Number of cells (y)      = " << nyc << endl;
  cout << "Number of cells (z)      = " << nzc << endl;
  cout << "Time step                = " << dt << endl;
  cout << "Number of cycles         = " << ncycles << endl;
  cout << "Results saved in  : " << SaveDirName << endl;
  cout << "Case type         : " << Case << endl;
  cout << "Simulation name   : " << SimName << endl;
  cout << "Poisson correction: " << PoissonCorrection << endl;
  cout << "---------------------" << endl;
  cout << "Check Simulation Constraints" << endl;
  cout << "---------------------" << endl;
  cout << "Accuracy Constraint:  " << endl;
  for (int i = 0; i < ns; i++) {
    cout << "u_th < dx/dt species " << i << ".....";
    if (uth[i] < (dx / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;

    cout << "v_th < dy/dt species " << i << "......";
    if (vth[i] < (dy / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;
  }
  cout << endl;
  cout << "Finite Grid Stability Constraint:  ";
  cout << endl;
  for (int is = 0; is < ns; is++) {
    if (uth[is] * dt / dx > .1)
      cout << "OK u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " > .1" << endl;
    else
      cout << "WARNING.  u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " < .1" << endl;

    if (vth[is] * dt / dy > .1)
      cout << "OK v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " > .1" << endl;
    else
      cout << "WARNING. v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " < .1"  << endl;

  }
  cout << "---------------------" << endl;
  cout << "MLMD variables" << endl;

  cout << "Number of MLMD grids          : " << Ngrids << endl;

  cout << "Grid levels:                  : " << endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout << gridLevel[ng] << "\t";
    }
  cout <<endl;

  cout << "Refinement factors            : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout << RF[ng] << "\t";
    }
  cout <<endl;
  
  cout << "Parent grid                   : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout << parentGrid[ng] << "\t";
    }
  cout <<endl;

  cout << "Grid origin, x y z direction (with respect to the COARSEST grid)      : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout <<"grid " << ng <<":\t"<< Ox_SW[ng] <<"\t" <<Oy_SW[ng] << "\t" <<Oz_SW[ng] <<endl;
    }
  cout <<endl;

  cout << "Spatial resolution, x y z direction      : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout <<"grid " << ng <<":\t"<< dx_mlmd[ng] <<"\t" <<dy_mlmd[ng] << "\t" <<dz_mlmd[ng] <<endl;
    }
  cout <<endl;
  
  cout << "Number of cells, x y z direction      : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      cout <<"grid " << ng <<":\t"<< nxc_mlmd[ng] <<"\t" <<nyc_mlmd[ng] << "\t" <<nzc_mlmd[ng] <<endl;
    }
  cout <<endl;

  cout << "---------------------" << endl;

}
/*! Print Simulation Parameters */
void Collective::save() {
  string temp;
  temp = SaveDirName + "/SimulationData.txt";
  ofstream my_file(temp.c_str());
  my_file << "---------------------------" << endl;
  my_file << "-  Simulation Parameters  -" << endl;
  my_file << "---------------------------" << endl;

  my_file << "Number of species    = " << ns << endl;
  for (int i = 0; i < ns; i++)
    my_file << "Number of particles per proc of species " << i << " = " << np[i] << "\t (MAX = " << npMax[i] << ")" << "  QOM = " << qom[i] << endl;
  my_file << "---------------------------" << endl;
  my_file << "x-Length                 = " << Lx << endl;
  my_file << "y-Length                 = " << Ly << endl;
  my_file << "z-Length                 = " << Lz << endl;
  my_file << "Number of cells (x)      = " << nxc << endl;
  my_file << "Number of cells (y)      = " << nyc << endl;
  my_file << "Number of cells (z)      = " << nzc << endl;
  my_file << "---------------------------" << endl;
  my_file << "Time step                = " << dt << endl;
  my_file << "Number of cycles         = " << ncycles << endl;
  my_file << "---------------------------" << endl;
  for (int is = 0; is < ns; is++){
    my_file << "rho init species   " << is << " = " << rhoINIT[is] << endl;
    my_file << "rho inject species " << is << " = " << rhoINJECT[is]  << endl;
  }
  my_file << "current sheet thickness  = " << delta << endl;
  my_file << "B0x                      = " << B0x << endl;
  my_file << "BOy                      = " << B0y << endl;
  my_file << "B0z                      = " << B0z << endl;
  my_file << "---------------------------" << endl;
  my_file << "Smooth                   = " << Smooth << endl;
  my_file << "GMRES error tolerance    = " << GMREStol << endl;
  my_file << "CG error tolerance       = " << CGtol << endl;
  my_file << "Mover error tolerance    = " << NiterMover << endl;
  my_file << "---------------------------" << endl;
  my_file << "Results saved in: " << SaveDirName << endl;
  my_file << "Restart saved in: " << RestartDirName << endl;
  my_file << "---------------------" << endl;

  my_file << "---------------------" << endl;
  my_file << "MLMD variables" << endl;

  my_file << "Number of MLMD grids          : " << Ngrids << endl;

  my_file << "Grid levels:                  : " << endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file << gridLevel[ng] << "\t";
    }
  my_file <<endl;

  my_file << "Refinement factors            : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file << RF[ng] << "\t";
    }
  my_file <<endl;
  
  my_file << "Parent grid                   : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file << parentGrid[ng] << "\t";
    }
  my_file <<endl;

  my_file << "Grid origin, x y z direction (with respect to the COARSEST grid)     : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file <<"grid " << ng <<":\t"<< Ox_SW[ng] <<"\t" <<Oy_SW[ng] << "\t" <<Oz_SW[ng] <<endl;
    }
  my_file <<endl;

  my_file << "Spatial resolution, x y z direction      : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file <<"grid " << ng <<":\t"<< dx_mlmd[ng] <<"\t" <<dy_mlmd[ng] << "\t" <<dz_mlmd[ng] <<endl;
    }
  my_file <<endl;
  
  my_file << "Number of cells, x y z direction      : " <<endl;
  for (int ng=0; ng < Ngrids; ng++)
    {
      my_file <<"grid " << ng <<":\t"<< nxc_mlmd[ng] <<"\t" <<nyc_mlmd[ng] << "\t" <<nzc_mlmd[ng] <<endl;
    }
  my_file <<endl;

  my_file << "---------------------" << endl;


  my_file.close();

}

/*! get the physical space dimensions */
int Collective::getDim() {
  return (dim);
}
/*! mlmd: use getLx_mlmd instead */
/*! get Lx */
/*double Collective::getLx() {
  return (Lx);
}*/
/*! get Ly */
/*double Collective::getLy() {
  return (Ly);
}*/
/*! get Lz */
/*double Collective::getLz() {
  return (Lz);
}*/
/*! end mlmd: use getLx_mlmd instead */ 
/*! get x_center */
double Collective::getx_center() {
  return (x_center);
}
/*! get y_center */
double Collective::gety_center() {
  return (y_center);
}
/*! get z_center */
double Collective::getz_center() {
  return (z_center);
}
/*! get L_square */
double Collective::getL_square() {
  return (L_square);
}
/*! mlmd: use getNxc_mlmd instead */
/*! get nxc */
/*int Collective::getNxc() {
  return (nxc);
}*/
/*! get nyx */
/*int Collective::getNyc() {
  return (nyc);
}*/
/*! get nzc */
/*!int Collective::getNzc() {
  return (nzc);
}*/
/*! get dx */
/*! mlmd: use getDx_mlmd instead*/
/*!double Collective::getDx() {
  return (dx);
}
/*! get dy */
/*double Collective::getDy() {
  return (dy);
}*/
/*! get dz */
/*double Collective::getDz() {
  return (dz);
}*/
/*! get the light speed */
double Collective::getC() {
  return (c);
}
/*! get the time step */
double Collective::getDt() {
  return (dt);
}
/*! get the decentering parameter */
double Collective::getTh() {
  return (th);
}
/*! get the smooth parameter */
double Collective::getSmooth() {
  return (Smooth);
}

/*! get the number of time cycles */
int Collective::getNcycles() {
  return (ncycles);
}
/*! get the number of species */
int Collective::getNs() {
  return (ns);
}
/*! get the number of particles per cell for species nspecies */
int Collective::getNpcel(int nspecies) {
  return (npcel[nspecies]);
}
/*! get the number of particles per cell for species nspecies - direction X */
int Collective::getNpcelx(int nspecies) {
  return (npcelx[nspecies]);
}
/*! get the number of particles per cell for species nspecies - direction Y */
int Collective::getNpcely(int nspecies) {
  return (npcely[nspecies]);
}
/*! get the number of particles per cell for species nspecies - direction Z */
int Collective::getNpcelz(int nspecies) {
  return (npcelz[nspecies]);
}
/*! get the number of particles for different species */
long Collective::getNp(int nspecies) {
  return (np[nspecies]);
}
/*! get maximum number of particles for different species */
long Collective::getNpMax(int nspecies) {
  return (npMax[nspecies]);
}
double Collective::getNpMaxNpRatio() {
  return (NpMaxNpRatio);
}
/*! get charge to mass ratio for different species */
double Collective::getQOM(int nspecies) {
  return (qom[nspecies]);
}
/*! get the background density for GEM challenge */
double Collective::getRHOinit(int nspecies) {
  return (rhoINIT[nspecies]);
}
/*! get the background density for GEM challenge */
double Collective::getRHOinject(int nspecies){
  return(rhoINJECT[nspecies]);
}
/*! get thermal velocity - Direction X */
double Collective::getUth(int nspecies) {
  return (uth[nspecies]);
}
/*! get thermal velocity - Direction Y */
double Collective::getVth(int nspecies) {
  return (vth[nspecies]);
}
/*! get thermal velocity - Direction Z */
double Collective::getWth(int nspecies) {
  return (wth[nspecies]);
}
/*! get beam velocity - Direction X */
double Collective::getU0(int nspecies) {
  return (u0[nspecies]);
}
/*! get beam velocity - Direction Y */
double Collective::getV0(int nspecies) {
  return (v0[nspecies]);
}
/*! get beam velocity - Direction Z */
double Collective::getW0(int nspecies) {
  return (w0[nspecies]);
}
/*! get Boundary Condition Particles: FaceXright */
int Collective::getBcPfaceXright() {
  return (bcPfaceXright);
}
/*! get Boundary Condition Particles: FaceXleft */
int Collective::getBcPfaceXleft() {
  return (bcPfaceXleft);
}
/*! get Boundary Condition Particles: FaceYright */
int Collective::getBcPfaceYright() {
  return (bcPfaceYright);
}
/*! get Boundary Condition Particles: FaceYleft */
int Collective::getBcPfaceYleft() {
  return (bcPfaceYleft);
}
/*! get Boundary Condition Particles: FaceZright */
int Collective::getBcPfaceZright() {
  return (bcPfaceZright);
}
/*! get Boundary Condition Particles: FaceZleft */
int Collective::getBcPfaceZleft() {
  return (bcPfaceZleft);
}
/*! get Boundary Condition Electrostatic Potential: FaceXright */
int Collective::getBcPHIfaceXright() {
  return (bcPHIfaceXright);
}
/*! get Boundary Condition Electrostatic Potential:FaceXleft */
int Collective::getBcPHIfaceXleft() {
  return (bcPHIfaceXleft);
}
/*! get Boundary Condition Electrostatic Potential:FaceYright */
int Collective::getBcPHIfaceYright() {
  return (bcPHIfaceYright);
}
/*! get Boundary Condition Electrostatic Potential:FaceYleft */
int Collective::getBcPHIfaceYleft() {
  return (bcPHIfaceYleft);
}
/*! get Boundary Condition Electrostatic Potential:FaceZright */
int Collective::getBcPHIfaceZright() {
  return (bcPHIfaceZright);
}
/*! get Boundary Condition Electrostatic Potential:FaceZleft */
int Collective::getBcPHIfaceZleft() {
  return (bcPHIfaceZleft);
}
/*! get Boundary Condition EM Field: FaceXright */
int Collective::getBcEMfaceXright() {
  return (bcEMfaceXright);
}
/*! get Boundary Condition EM Field: FaceXleft */
int Collective::getBcEMfaceXleft() {
  return (bcEMfaceXleft);
}
/*! get Boundary Condition EM Field: FaceYright */
int Collective::getBcEMfaceYright() {
  return (bcEMfaceYright);
}
/*! get Boundary Condition EM Field: FaceYleft */
int Collective::getBcEMfaceYleft() {
  return (bcEMfaceYleft);
}
/*! get Boundary Condition EM Field: FaceZright */
int Collective::getBcEMfaceZright() {
  return (bcEMfaceZright);
}
/*! get Boundary Condition EM Field: FaceZleft */
int Collective::getBcEMfaceZleft() {
  return (bcEMfaceZleft);
}
/*! Get GEM Challenge parameters */
double Collective::getDelta() {
  return (delta);
}
double Collective::getB0x() {
  return (B0x);
}
double Collective::getB0y() {
  return (B0y);
}
double Collective::getB0z() {
  return (B0z);
}
double Collective::getB1x(){
  return (B1x);
}
double Collective::getB1y(){
  return (B1y);
}
double Collective::getB1z(){
  return (B1z);
}
/*! get the boolean value for verbose results */
bool Collective::getVerbose() {
  return (verbose);
}
/*! get the boolean value for TrackParticleID */
bool Collective::getTrackParticleID(int nspecies) {
  return (TrackParticleID[nspecies]);
}
int Collective::getRestart_status() {
  return (restart_status);
}
/*! get SaveDirName */
string Collective::getSaveDirName() {
  return (SaveDirName);
}
/*! get RestartDirName */
string Collective::getRestartDirName() {
  return (RestartDirName);
}
/*! get inputfile */
string Collective::getinputfile() {
  return (inputfile);
}
/*! get initfile */
string Collective::getinitfile() {
  return (initfile);
}
/*! get Case type */
string Collective::getCase() {
  return (Case);
}
/*! get Particle initialization type */
string Collective::getPartInit() {
  return (PartInit);
}
/*! get simulation name */
string Collective::getSimName() {
  return (SimName);
}
/*! get output writing method */
string Collective::getWriteMethod() {
  return (wmethod);
}
/*! get Poisson correction flag */
string Collective::getPoissonCorrection() {
  return (PoissonCorrection);
}
/*! get initial solution flag */
bool Collective::getSolInit() {
  return (SolInit);
}
/*! get last_cycle */
int Collective::getLast_cycle() {
  return (last_cycle);
}
/*! get the velocity of injection of the plasma from the wall */
double Collective::getVinj() {
  return (Vinj);
}
/*! get the converging tolerance for CG solver */
double Collective::getCGtol() {
  return (CGtol);
}
/*! get the converging tolerance for GMRES solver */
double Collective::getGMREStol() {
  return (GMREStol);
}
/*! get the numbers of iteration for the PC mover */
int Collective::getNiterMover() {
  return (NiterMover);
}
/*! output of fields */
int Collective::getFieldOutputCycle() {
  return (FieldOutputCycle);
}
/*! output of particles */
int Collective::getParticlesOutputCycle() {
  return (ParticlesOutputCycle);
}
/*! restart cycle */
int Collective::getRestartOutputCycle() {
  return (RestartOutputCycle);
}
/*! output of fields */
int Collective::getDiagnosticsOutputCycle() {
  return (DiagnosticsOutputCycle);
}


/*! MLMD specific function */

/*! MLMD gets */
int Collective::getNgrids() {
  return Ngrids;
}
int Collective::getgridLevel(int numgrid) {
  return gridLevel[numgrid];
}
int Collective::getRF(int numgrid) {
  return RF[numgrid];
}
int Collective::getparentGrid(int numgrid) {
  return parentGrid[numgrid];
}
double Collective::getOx_SW(int numgrid) {
  return Ox_SW[numgrid];
}
double Collective::getOy_SW(int numgrid) {
  return Oy_SW[numgrid];
}
double Collective::getOz_SW(int numgrid) {
  return Oz_SW[numgrid];
}
double Collective::getOx_P(int numgrid) {
  return Ox_P[numgrid];
}
double Collective::getOy_P(int numgrid) {
  return Oy_P[numgrid];
}
double Collective::getOz_P(int numgrid) {
  return Oz_P[numgrid];
}
double Collective::getDx_mlmd(int numgrid) {
  return dx_mlmd[numgrid];
}
double Collective::getDy_mlmd(int numgrid) {
  return dy_mlmd[numgrid];
}
double Collective::getDz_mlmd(int numgrid) {
  return dz_mlmd[numgrid];
}
int Collective::getNxc_mlmd(int numgrid) {
  return nxc_mlmd[numgrid];
}
int Collective::getNyc_mlmd(int numgrid) {
  return nyc_mlmd[numgrid];
}
int Collective::getNzc_mlmd(int numgrid) {
  return nzc_mlmd[numgrid];
}
double Collective::getLx_mlmd(int numgrid) {
  return Lx_mlmd[numgrid];
}
double Collective::getLy_mlmd(int numgrid) {
  return Ly_mlmd[numgrid];
}
double Collective::getLz_mlmd(int numgrid) {
  return Lz_mlmd[numgrid];
}
int Collective::getTopologyType() {
  return TopologyType;
}
bool Collective::getMLMDVerbose() {
  return MLMDVerbose;
}
int Collective::getChildrenNum(int numgrid) {
  return childrenNum[numgrid];
}
int Collective::getChildrenGrids(int numgrid, int childnum) {
  if (childnum > childrenNum[numgrid] -1){
    cout << "WARNING!!! Collective::getChildrenGrids asked for a child grid which does not exists! WARNING!!" << endl;
    // end: later: manage this better //
  } else {
    return childrenGrids[numgrid][childnum];
  }
}
int Collective::getParentGrid(int numgrid) {
  if (numgrid==0){
    cout << "WARNING!!!!   Collective::getParentGrid(0)   WARNING!!!" <<endl;  
  }
  return parentGrid[numgrid];
}
int Collective::getMLMD_BC() {return MLMD_BC;}
int Collective::getMLMD_PROJECTION() {return MLMD_PROJECTION;}
int Collective::getMLMD_ParticleREPOPULATION() {return MLMD_ParticleREPOPULATION;}

/*! end MLMD gets */
/*! a first sanity check on MLMD inputs, called at the end of the constructor */
void Collective::checkMLMDinputs() {

  int rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  /*! check that the grids hierarchy is in increasing order in the inputfile; 
    it's just easier to process things that way  */
  for (int ng=1; ng< Ngrids; ng++) {
    if (gridLevel[ng] < gridLevel[ng-1]) { 
      if (rank ==0)
	cout<<" Please order the grid levels in increasing order in the inputfile, aborting...\n"<<flush;
      abort(); 
    } 
  } // end for grids

  /*! check that the refined grid is included in the parent grid  */
  bool ok, inBetween;
  double tol= 2*DBL_EPSILON;
  for (int ng=1; ng< Ngrids; ng++) {
    
    /*! x dir */
    // NB: Ox_P is in terms of the parent grid, not of the coarsest grid: hence we check it the Ox_P is between 0 and the length of the parent grid
    inBetween = isInBetween(Ox_P[ng], 0.0, Lx_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Ox_P[ng] < dx_mlmd[parentGrid[ng]])) {	
      if (rank==0)
	{
	  cout<<" Problem with x-position of grid " << ng <<":" <<endl <<"Lx parent= " << Lx_mlmd[parentGrid[ng]] << ", Ox_P=" << Ox_P[ng] <<", Ox_P+Lx=" << Ox_P[ng] + Lx_mlmd[ng] << ", keep Ox_P>" <<dx_mlmd[parentGrid[ng]]<< ", ok: " << ok << ",inBetween:" <<inBetween << ", aborting" << flush;
	}
      abort(); }
    
    inBetween = isInBetween(Ox_P[ng]+ Lx_mlmd[ng], 0.0, Lx_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Ox_P[ng] + Lx_mlmd[ng] > Lx_mlmd[parentGrid[ng]] - dx_mlmd[parentGrid[ng]]    )) {
      if (rank==0)
	{
	  cout<<" Problem with x-position of grid " << ng <<":" <<endl <<"Lx parent= " << Lx_mlmd[parentGrid[ng]] << ", Ox_P=" << Ox_P[ng] <<", Ox_P+Lx=" << Ox_P[ng] + Lx_mlmd[ng] << ", keep Ox_P + Lx<" <<Lx_mlmd[parentGrid[ng]]- dx_mlmd[parentGrid[ng]] << ", ok: " << ok << ", inBetween: " << inBetween << ", aborting" << flush;
	}
      abort(); }
    /*! y dir */
    inBetween = isInBetween(Oy_P[ng], 0.0, Ly_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Oy_P[ng] < dy_mlmd[parentGrid[ng]])) {
      if (rank==0)
	cout<<" Problem with y-position of grid " << ng <<":" <<endl <<"Ly parent= " << Ly_mlmd[parentGrid[ng]] << ", Oy_P=" << Oy_P[ng] <<", Oy_P+Ly=" << Oy_P[ng] + Ly_mlmd[ng] << ", keep Oy_P>" <<dy_mlmd[parentGrid[ng]]<< ", aborting" << flush;
      abort(); }

    inBetween = isInBetween(Oy_P[ng]+ Ly_mlmd[ng], 0.0, Ly_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Oy_P[ng] + Ly_mlmd[ng] > Ly_mlmd[parentGrid[ng]] - dy_mlmd[parentGrid[ng]]    )) {
      cout<<" Problem with y-position of grid " << ng <<":" <<endl <<"Ly parent= " << Ly_mlmd[parentGrid[ng]] << ", Oy_P=" << Oy_P[ng] <<", Oy_P+Ly=" << Oy_P[ng] + Ly_mlmd[ng] << ", keep Oy_P + Ly<" <<Ly_mlmd[parentGrid[ng]]- dy_mlmd[parentGrid[ng]]<< ", aborting" << flush;
      abort(); }
    
    /*! z dir */
    inBetween = isInBetween(Oz_P[ng], 0.0, Lz_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Oz_P[ng] < dz_mlmd[parentGrid[ng]])) {
      if (rank==0)
	cout<<" Problem with z-position of grid " << ng <<":" <<endl <<"Lz parent= " << Lz_mlmd[parentGrid[ng]] << ", Oz_P=" << Oz_P[ng] <<", Oz_P+Lz=" << Oz_P[ng] + Lz_mlmd[ng] << ", keep Oz_P>" <<dz_mlmd[parentGrid[ng]] <<", aborting" << flush;
      abort(); }

    inBetween = isInBetween(Oz_P[ng]+ Lz_mlmd[ng], 0.0, Lz_mlmd[parentGrid[ng]], tol, &ok);
    if ((! (ok && inBetween)) or (Oz_P[ng] + Lz_mlmd[ng] > Lz_mlmd[parentGrid[ng]] - dz_mlmd[parentGrid[ng]]    )) {
      if (rank==0)
	cout<<" Problem with z-position of grid " << ng <<":" <<endl <<"Lz parent= " << Lz_mlmd[parentGrid[ng]] << ", Oz_P=" << Oz_P[ng] <<", Oz_P+Lz=" << Oz_P[ng] + Lz_mlmd[ng] << ", keep Oz_P + Lz <" <<Lz_mlmd[parentGrid[ng]]- dz_mlmd[parentGrid[ng]]<< ", aborting" << flush;
      abort(); }
  } // end for grids

  /*! check that nxc/XLEN, nyc/YLEN, nzc/ZLEN are integer for all grids  */
  for (int ng=1; ng< Ngrids; ng++) {
    bool rem;
    if (nxc_mlmd[ng]%XLEN) 
      { if (rank==0)
	  cout << " Please check the # of cells in the x-dir of grid " <<ng <<", aborting ..." << flush; 
	abort(); }
    
    if (nyc_mlmd[ng]%YLEN) 
      { if (rank==0)
	  cout << " Please check the # of cells in the x-dir of grid " <<ng <<", aborting ..." << flush; 
	abort(); }
    
    if (nzc_mlmd[ng]%ZLEN) 
      { if (rank==0)
	  cout << " Please check the # of cells in the x-dir of grid " <<ng <<", aborting ..." << flush; 
	abort(); }
  } // end for grids
  
  /*! check that the requested topology type has been implemented */
  if (TopologyType)
    {
      if (rank==0){
	cout << "You requested TopologyType= " << TopologyType << endl;
	cout << "Only TopologyType = 0 is implemented, aborting ..." << flush;}
      abort();
    }
}
/*! end MLMD specific function */
