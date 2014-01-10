
#include "H5hut-io.h"

/* ==================== */
/* Auxiliary structures */
/* ==================== */

H5hutpart::H5hutpart(){
  allocated = false;
}

H5hutpart::~H5hutpart(){
  delete [] q;
  delete [] x;
  delete [] y;
  delete [] z;
  delete [] u;
  delete [] v;
  delete [] w;
}

void H5hutpart::memalloc(long long n){
  np = n;

  if (allocated == false) {
    q = new double[np];
    x = new double[np];
    y = new double[np];
    z = new double[np];
    u = new double[np];
    v = new double[np];
    w = new double[np];
    allocated = true;
  }
}

void H5hutpart::init(long long n, double chg){

  if (!allocated) memalloc(n);

  double LO = 0.1;
  double HI = 99.9;
  for (int i=0; i<np; i++){
    q[i] = chg;
    x[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
    y[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
    z[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
    u[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
    v[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
    w[i] = LO + double(rand()) / double(RAND_MAX/(HI-LO));
  }
}

void H5hutpart::info(){
  int i = int(np/2);
  std::cout << "H5hutpart with np=" << np << ", ex. position: " << x[i] << " " << y[i] << " " << z[i] << std::endl;
  std::cout << "H5hutpart with np=" << np << ", 1st position: " << x[0] << " " << y[0] << " " << z[0] << std::endl;
  std::cout << "H5hutpart with np=" << np << ", Lst position: " << x[np-1] << " " << y[np-1] << " " << z[np-1] << std::endl;
  std::cout << std::endl;
}

/* ====================== */
/*        OUTPUT          */
/* ====================== */

void H5output::SetBaseName(std::string name){
  basename = name;
}

void H5output::OpenFieldsCellFile(int nspec, int ntcx, int ntcy, int ntcz, int cycle, int *coord, int *dimns, MPI_Comm CART_COMM){

  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(5) << cycle;
  filename = basename + "-Fields" + "_" + filenmbr.str() + ".h5";

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  fldsfile = H5OpenFile(filename.c_str(), H5_O_WRONLY, CART_COMM);
  H5SetStep(fldsfile,0);
  H5WriteStepAttribInt32(fldsfile, "nspec", &nspec, 1);
  H5Block3dSetGrid(fldsfile, dimns[0], dimns[1], dimns[2]);
  H5Block3dSetDims(fldsfile, ntcx/dimns[0], ntcy/dimns[1], ntcz/dimns[2]);
  H5Block3dSetHalo(fldsfile, 1, 1, 1);

  int irange[2];
  int jrange[2];
  int krange[2];

  irange[0] =   coord[0]   * ntcx/dimns[0];
  irange[1] = ((coord[0]+1)* ntcx/dimns[0]) -1;
  jrange[0] =   coord[1]   * ntcy/dimns[1];
  jrange[1] = ((coord[1]+1)* ntcy/dimns[1]) -1;
  krange[0] =   coord[2]   * ntcz/dimns[2];
  krange[1] = ((coord[2]+1)* ntcz/dimns[2]) -1;

  /* --------------- */
  /* Write EM fields */
  /* --------------- */
  H5Block3dSetView(fldsfile, irange[0], irange[1],
                             jrange[0], jrange[1],
                             krange[0], krange[1]);

}

void H5output::CloseFieldsCellFile(){
  H5CloseFile(fldsfile);
}

void H5output::WriteFieldsCell(double ***field, std::string fname, int nxc, int nyc, int nzc){

  h5_float64_t*     buffer;

  buffer = new h5_float64_t[nzc*nyc*nxc];

  int n = 0;
  for (int k=0; k<nzc; k++) {
    for (int j=0; j<nyc; j++) {
      for (int i=0; i<nxc; i++) {
        buffer[n++] = field[i][j][k];
      }
    }
  }

  H5Block3dWriteScalarFieldFloat64(fldsfile, fname.c_str(), buffer);

  delete [] buffer;
}

void H5output::OpenPartclFile(int nspec, int cycle, MPI_Comm CART_COMM){
  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(5) << cycle;
  filename = basename + "-Partcl" + "_" + filenmbr.str() + ".h5";

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  partfile = H5OpenFile(filename.c_str(), H5_O_WRONLY, CART_COMM);
  H5SetStep(partfile,0);
  H5WriteStepAttribInt32(partfile, "nspec", &nspec, 1);

}

void H5output::ClosePartclFile(){
  H5CloseFile(partfile);
}

void H5output::WriteParticles(int ispec, long long np, double *q, double *x, double *y, double *z, double *u, double *v, double *w, MPI_Comm CART_COMM){

  /* --------------------------------------------------------------------- */
  /* Find out the total number of particles of species i in all the domain */
  /* --------------------------------------------------------------------- */

  long long ntpart;

  MPI_Allreduce(&np, &ntpart, 1, MPI_LONG_LONG, MPI_SUM, CART_COMM);
  const h5_int64_t ntp = ntpart;

  /* --------------------------------------------- */
  /* Write the number of particles as an attribute */
  /* --------------------------------------------- */

  std::stringstream sstm;

  sstm << "npart_" << ispec;
  std::string nparti = sstm.str();
  H5WriteStepAttribInt64(partfile, nparti.c_str(), &ntp, 1);
  sstm.str("");

  H5PartSetNumParticles(partfile,np);

  /* ------------------ */
  /* Write the datasets */
  /* ------------------ */

  sstm << "q_" << ispec;
  std::string dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),q);
  sstm.str("");

  sstm << "x_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),x);
  sstm.str("");

  sstm << "y_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),y);
  sstm.str("");

  sstm << "z_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),z);
  sstm.str("");

  sstm << "u_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),u);
  sstm.str("");

  sstm << "v_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),v);
  sstm.str("");

  sstm << "w_" << ispec;
  dtset = sstm.str();
  H5PartWriteDataFloat64(partfile,dtset.c_str(),w);
  sstm.str("");

}

/* ====================== */
/*         INPUT          */
/* ====================== */

void H5input::OpenPartclFile(int ns, std::string fn){
  nspec = ns;
  part = new H5hutpart[nspec];
  filename = fn;
}

void H5input::ReadParticles(int ndim, int rank, int nproc, int *dimns, double *L, MPI_Comm CART_COMM){
  LoadParticles(ndim, rank, nproc, dimns, L, CART_COMM);
  InitParticles(ndim, rank, CART_COMM);
}

void H5input::ClosePartclFile(){
  delete [] part;
}

void H5input::LoadLocalParticles(long long *np,
                        double *q_loc,
                        double *x_loc,
                        double *y_loc,
                        double *z_loc,
                        double *u_loc,
                        double *v_loc,
                        double *w_loc) {

  long long n;

  n = 0;
  for (int s=0; s<nspec; s++){

    part[s].memalloc(np[s]);

    std::copy(q_loc+n, q_loc+n+np[s], part[s].getQ());
    std::copy(x_loc+n, x_loc+n+np[s], part[s].getX());
    std::copy(y_loc+n, y_loc+n+np[s], part[s].getY());
    std::copy(z_loc+n, z_loc+n+np[s], part[s].getZ());
    std::copy(u_loc+n, u_loc+n+np[s], part[s].getU());
    std::copy(v_loc+n, v_loc+n+np[s], part[s].getV());
    std::copy(w_loc+n, w_loc+n+np[s], part[s].getW());

    n += np[s];
  }
}

void H5input::FindLocalParticles(int nproc, int ndim, h5_int64_t *h5npart, long long nop, int *dimns, double *L, MPI_Comm CART_COMM,
                        double *q,
                        double *x, double *y, double *z,
                        double *u, double *v, double *w){

  double     *q_loc;
  double     *x_loc;
  double     *y_loc;
  double     *z_loc;
  double     *u_loc;
  double     *v_loc;
  double     *w_loc;
  bool       *b;
  long long  locnp[nspec][nproc];
  long long  *locnop;
  long long  *np;
  long long  p;
  long long  nn;
  long long  MAX_NOP;
  int        irank;
  int        *pcoord;
  int        *inproc;

  pcoord = new int[ndim];
  inproc = new int[nop];
  np     = new long long[nspec+1];
  locnop = new long long[nproc];

  MAX_NOP = 1.1 * (nop/nproc);

  /* ------------------- */
  /* Initialize counters */
  /* ------------------- */

  for (int iproc=0; iproc<nproc; iproc++) {
    locnop[iproc] = 0;
    for (int i=0; i<nspec; i++) {
      locnp[i][iproc] = 0;
    }
  }

  /* ------------------------- */
  /* Tag particles with proc # */
  /* ------------------------- */

  nn = 0;
  for (int i=0; i<nspec; i++) {
    for (long long n=0; n < h5npart[i]; n++) {

      p = nn + n;

      /* --------------------------------------- */
      /* Get processor coordinate for particle p */
      /* --------------------------------------- */

      for (int d=0; d<ndim; d++) {
        double r;
        if (d==0) r = x[p];
        if (d==1) r = y[p];
        if (d==2) r = z[p];

        pcoord[d] = int(r/(L[d]/dimns[d]));
      }

      /* ----------------------------------------------------- */
      /* Detect processor rank and increment particle counters */
      /* ----------------------------------------------------- */

      MPI_Cart_rank(CART_COMM, pcoord, &irank);
      inproc[p] = irank;
      locnp[i][irank]++;
      locnop[irank]++;

    }

    nn += h5npart[i];
  }

  /* ----------------------- */
  /* Create the local arrays */
  /* ----------------------- */

  MAX_NOP = 0;
  for (int iproc=0; iproc<nproc; iproc++){
    if (locnop[iproc]>MAX_NOP) MAX_NOP = locnop[iproc];
  }

  q_loc = new double[MAX_NOP];
  x_loc = new double[MAX_NOP];
  y_loc = new double[MAX_NOP];
  z_loc = new double[MAX_NOP];
  u_loc = new double[MAX_NOP];
  v_loc = new double[MAX_NOP];
  w_loc = new double[MAX_NOP];

  /* ------------------------------- */
  /* Fill vectors for each processor */
  /* ------------------------------- */

  for (int iproc=0; iproc<nproc; iproc++) {

    /* ---------------- */
    /* Fill the vectors */
    /* ---------------- */

    p = 0;
    for (long long n=0; n<nop; n++) {
      if (inproc[n]==iproc) {
        q_loc[p] = q[n];
        x_loc[p] = x[n];
        y_loc[p] = y[n];
        z_loc[p] = z[n];
        u_loc[p] = u[n];
        v_loc[p] = v[n];
        w_loc[p] = w[n];
        p++;
      }
    }

    /* ------------------------------------------------------------------- */
    /* Copy particles in procesor 0 and send particles to other processors */
    /* ------------------------------------------------------------------- */

    for (int i=0; i<nspec; i++) np[i] = locnp[i][iproc];
    np[nspec] = locnop[iproc]; // The last element in np contains the total (all species combined) number of particles.

    if (iproc==0)  LoadLocalParticles(np, q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);
    else{
      MPI_Send(np    , nspec+1      , MPI_LONG_LONG, iproc, 0, CART_COMM);
      MPI_Send(q_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(x_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(y_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(z_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(u_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(v_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
      MPI_Send(w_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
    }
  }

  delete [] q_loc;
  delete [] x_loc;
  delete [] y_loc;
  delete [] z_loc;
  delete [] u_loc;
  delete [] v_loc;
  delete [] w_loc;

  delete [] pcoord;
  delete [] inproc;
  delete [] np;
  delete [] locnop;

}

void H5input::InitParticles(int ndim, int rank, MPI_Comm CART_COMM){

  long long  *np;
  long long  nop;
  double     *q_loc;
  double     *x_loc;
  double     *y_loc;
  double     *z_loc;
  double     *u_loc;
  double     *v_loc;
  double     *w_loc;

  if (rank!=0) {
    /* ---------------------------------------- */
    /* Receive particle information from proc 0 */
    /* ---------------------------------------- */

    np = new long long[nspec+1];

    MPI_Recv(np, nspec+1, MPI_LONG_LONG, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    nop = np[nspec];

    /* ----------------------- */
    /* Create the local arrays */
    /* ----------------------- */

    q_loc = new double[nop];
    x_loc = new double[nop];
    y_loc = new double[nop];
    z_loc = new double[nop];
    u_loc = new double[nop];
    v_loc = new double[nop];
    w_loc = new double[nop];

    MPI_Recv(q_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(x_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(y_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(z_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(u_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(v_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
    MPI_Recv(w_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);

    LoadLocalParticles(np, q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);

    delete [] q_loc;
    delete [] x_loc;
    delete [] y_loc;
    delete [] z_loc;
    delete [] u_loc;
    delete [] v_loc;
    delete [] w_loc;
  }
}

void H5input::LoadParticles(int ndim, int rank, int nproc, int *dimns, double *L, MPI_Comm CART_COMM){

  int        h5nspec;
  int        ranks_rdr[1] = {0};
  int        *ranks_rst;
  long long  offset;
  long long  nop;
  h5_int64_t *h5npart;
  double     *q;
  double     *x;
  double     *y;
  double     *z;
  double     *u;
  double     *v;
  double     *w;

  std::stringstream sstm;

  MPI_Comm  READER_COMM;
  MPI_Group org_grp;
  MPI_Group reader_grp;

  part = new H5hutpart[nspec];

  /* -------------------------------------------------- */
  /* Create a new communicator for the reader processor */
  /* This way only processor 0 reads the whole file.    */
  /* It would be better if each processor read one part */
  /* of the file...currently not possible with H5hut.   */
  /* -------------------------------------------------- */

  MPI_Comm_group (CART_COMM, &org_grp);
  ranks_rst = new int[nproc-1];
  for (int i=0; i<nproc; i++)  ranks_rst[i] = i;

  if (rank==0) MPI_Group_incl (org_grp, 1, ranks_rdr, &reader_grp);
  else         MPI_Group_incl (org_grp, nproc-1, ranks_rst, &reader_grp);

  delete [] ranks_rst;

  MPI_Comm_create(CART_COMM, reader_grp, &READER_COMM);

  if (rank==0) {

    /* -------------- */
    /* Open HDF5 file */
    /* -------------- */

    h5_file_t *file = H5OpenFile(filename.c_str(), H5_O_RDONLY, READER_COMM);
    H5SetStep(file, 0);
    H5ReadStepAttribInt32(file, "nspec", &h5nspec);

    if (nspec!=h5nspec) {
      std::cout << "ERROR in ReadParallelParticles: the number of species in the initial file" << std::endl;
      std::cout << "                                does not match the number of species requested." << std::endl;
      abort();
    }

    std::cout << "Reading "<< nspec << " species." << std::endl;

    h5npart = new h5_int64_t[nspec];

    nop = 0;
    for (int i=0; i<nspec; i++) {
      sstm << "npart_" << i;
      std::string nparti = sstm.str();
      H5ReadStepAttribInt64(file, nparti.c_str(), &h5npart[i]);
      sstm.str("");
      nop += h5npart[i];
    }

    q = new double[nop];
    x = new double[nop];
    y = new double[nop];
    z = new double[nop];
    u = new double[nop];
    v = new double[nop];
    w = new double[nop];

    offset = 0;
    for (int i=0; i<nspec; i++) {

      std::cout << " Reading a total of " << h5npart[i] << " particles of species " << i << std::endl;

      sstm << "q_" << i;
      std::string dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),q+offset);
      sstm.str("");

      sstm << "x_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),x+offset);
      sstm.str("");

      sstm << "y_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),y+offset);
      sstm.str("");

      sstm << "z_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),z+offset);
      sstm.str("");

      sstm << "u_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),u+offset);
      sstm.str("");

      sstm << "v_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),v+offset);
      sstm.str("");

      sstm << "w_" << i;
      dtset = sstm.str();
      H5PartReadDataFloat64(file,dtset.c_str(),w+offset);
      sstm.str("");

      offset += h5npart[i];

    }

    H5CloseFile(file);

    FindLocalParticles(nproc, ndim, h5npart, nop, dimns, L, CART_COMM, q, x, y, z, u, v, w);

    delete [] q;
    delete [] x;
    delete [] y;
    delete [] z;
    delete [] u;
    delete [] v;
    delete [] w;

  }

}
