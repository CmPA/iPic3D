
#include <fstream>
#include <cmath>

#include "H5hut-io.h"

/* ==================== */
/* Auxiliary structures */
/* ==================== */

H5hutpart::H5hutpart(){
  allocated = false;
}

H5hutpart::~H5hutpart(){
}

void H5hutpart::memalloc(long long n){

  SetNp(n);

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

void H5hutpart::init(long long n, double chg, int Lx, int Ly, int Lz){

  if (!allocated) memalloc(n);

  double LO = -50.0;
  double HI =  50.0;
  for (int i=0; i<np; i++){
    q[i] = chg;
    x[i] = 0.0 + double(rand()) / double(RAND_MAX/(Lx-0.0));
    y[i] = 0.0 + double(rand()) / double(RAND_MAX/(Ly-0.0));
    z[i] = 0.0 + double(rand()) / double(RAND_MAX/(Lz-0.0));
    u[i] = LO  + double(rand()) / double(RAND_MAX/(HI-LO));
    v[i] = LO  + double(rand()) / double(RAND_MAX/(HI-LO));
    w[i] = LO  + double(rand()) / double(RAND_MAX/(HI-LO));
  }
}

void H5hutpart::info(){
  int i = int(np/2);
  std::cout << "[PHDF5-io]" << " np=" << np << ", ex. position: " << x[i] << " " << y[i] << " " << z[i] << std::endl;
  std::cout << "[PHDF5-io]" << " np=" << np << ", 1st position: " << x[0] << " " << y[0] << " " << z[0] << std::endl;
  std::cout << "[PHDF5-io]" << " np=" << np << ", Lst position: " << x[np-1] << " " << y[np-1] << " " << z[np-1] << std::endl;
  std::cout << std::endl;
}

/* ====================== */
/*        OUTPUT          */
/* ====================== */


void H5output::SetNameCycle(std::string name, int c){
  basename = name;
  cycle    = c;
}

void H5output::OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *pdims, MPI_Comm CART_COMM){

  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(6) << cycle;
  filename = basename + "-Fields" + "_" + filenmbr.str() + ".h5";

  int d = -1;
  if (dtype=="Cell") d = 0;

  ntx = ntx + d;
  nty = nty + d;
  ntz = ntz + d;

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  fldsfile = H5OpenFile(filename.c_str(), H5_O_RDWR, CART_COMM);
  H5SetStep(fldsfile,0);
  H5WriteStepAttribInt32(fldsfile, "nspec", &nspec, 1);
  H5Block3dSetGrid(fldsfile, pdims[0], pdims[1], pdims[2]);
  H5Block3dSetDims(fldsfile, ntx/pdims[0], nty/pdims[1], ntz/pdims[2]);
  H5Block3dSetHalo(fldsfile, 1, 1, 1);

  int irange[2];
  int jrange[2];
  int krange[2];

  int nnx = ntx/pdims[0];
  int nny = nty/pdims[1];
  int nnz = ntz/pdims[2];

  d = 0;
  if (dtype=="Cell") d = -1;   // Yes, this line is the oposite of the previous 'if'

  irange[0] = coord[0]      * nnx;
  irange[1] =(coord[0] + 1) * nnx + d;
  jrange[0] = coord[1]      * nny;
  jrange[1] =(coord[1] + 1) * nny + d;
  krange[0] = coord[2]      * nnz;
  krange[1] =(coord[2] + 1) * nnz + d;

  /* -------------- */
  /* Set the "view" */
  /* -------------- */
  H5Block3dSetView(fldsfile, irange[0], irange[1],
                             jrange[0], jrange[1],
                             krange[0], krange[1]);

}

void H5output::CloseFieldsFile(){
  H5CloseFile(fldsfile);
}

void H5output::WriteFields(double ***field, std::string fname, int nx, int ny, int nz, int rank){

  h5_float64_t*     buffer;

  buffer = new h5_float64_t[nz*ny*nx];

  std::ofstream myfile;

  if (rank!=-1) {
    std::stringstream  ss;
    ss << "proc-" << rank << ".txt";
    std::string filename = ss.str();
    myfile.open(filename.c_str());
  }

  int n = 0;
  for (int k=1; k<nz-1; k++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=1; i<nx-1; i++) {
        buffer[n++] = field[i][j][k];
        if (rank!=-1) myfile << i << " " << j << " " << k << " : " << field[i][j][k] << std::endl;
      }
    }
  }
  if (rank!=-1) myfile.close();

  H5Block3dWriteScalarFieldFloat64(fldsfile, fname.c_str(), buffer);

  delete [] buffer;
}

void H5output::OpenPartclFile(int nspec, MPI_Comm CART_COMM){
  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(6) << cycle;
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

void H5input::SetNameCycle(std::string name, int rc){
  basename = name;
  recycle  = rc;
}

void H5input::OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *pdims, MPI_Comm CART_COMM){

  int file_nspec;

  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(6) << recycle;
  filename = basename + "-Fields" + "_" + filenmbr.str() + ".h5";

  int d = -1;
  if (dtype=="Cell") d = 0;

  ntx = ntx + d;
  nty = nty + d;
  ntz = ntz + d;

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  fldsfile = H5OpenFile(filename.c_str(), H5_O_RDONLY, CART_COMM);
  H5SetStep(fldsfile, 0);
  H5ReadStepAttribInt32(fldsfile, "nspec", &file_nspec);

  if (file_nspec!=nspec){
   std::cout << "[PHDF5-io]" << " ERROR: The number of species nspec=" << nspec << std::endl;
   std::cout << "[PHDF5-io]" << "        does not correspont to the number of species in the initial file." << std::endl;
   std::cout << "[PHDF5-io]" << "        file_nspec = " << file_nspec << std::endl; 
   abort();
  }

  /* ---------------------------------------------------------------------------- */
  /* Verify that the data in the file corresponds to the current case parametters */
  /* TODO: in a future version instead of stoping the code we should include the  */
  /*       input file into the initial file so the parameters are not separated.  */
  /* ---------------------------------------------------------------------------- */

  h5_size_t  f_rank, f_dims[3], e_rank;
  h5_int64_t f_type;

  H5BlockGetFieldInfoByName(fldsfile, "Ex", &f_rank, f_dims, &e_rank, &f_type);

  int ndim = 3; // By default 2D is considered as a 'flat' 3D

  if (f_rank!=ndim) {
   std::cout << "[PHDF5-io]" << " ERROR: The number of dimensions =" << ndim << std::endl;
   std::cout << "[PHDF5-io]" << "        does not correspont to the number of dimensions in the initial file." << std::endl;
   std::cout << "[PHDF5-io]" << "        f_ndim = " << f_rank << std::endl; 
   abort();
  }

  if (f_dims[0]!=ntx-d) {
   std::cout << "[PHDF5-io]" << " ERROR: The number of cells in X =" << ntx << std::endl;
   std::cout << "[PHDF5-io]" << "        does not correspont to the number of cells in the initial file." << std::endl;
   std::cout << "[PHDF5-io]" << "        ncell x = " << f_dims[0] << std::endl; 
   abort();
  }

  if (f_dims[1]!=nty-d) {
   std::cout << "[PHDF5-io]" << " ERROR: The number of cells in Y =" << nty << std::endl;
   std::cout << "[PHDF5-io]" << "        does not correspont to the number of cells in the initial file." << std::endl;
   std::cout << "[PHDF5-io]" << "        ncell y = " << f_dims[1] << std::endl; 
   abort();
  }

  if (f_dims[2]!=ntz-d) {
   std::cout << "[PHDF5-io]" << " ERROR: The number of cells in Z =" << ntz << std::endl;
   std::cout << "[PHDF5-io]" << "        does not correspont to the number of cells in the initial file." << std::endl;
   std::cout << "[PHDF5-io]" << "        ncell z = " << f_dims[2] << std::endl; 
   abort();
  }

  /* -------------------- */
  /* Set Block properties */
  /* -------------------- */

  H5Block3dSetGrid(fldsfile, pdims[0], pdims[1], pdims[2]);
  H5Block3dSetDims(fldsfile, ntx/pdims[0], nty/pdims[1], ntz/pdims[2]);
  H5Block3dSetHalo(fldsfile, 1, 1, 1);

  int irange[2];
  int jrange[2];
  int krange[2];

  d = 0;
  if (dtype=="Cell") d = -1;  // Yes, this is the oposite of the first 'if'

  irange[0] =   coord[0]   * ntx/pdims[0];
  irange[1] = ((coord[0]+1)* ntx/pdims[0]) + d;
  jrange[0] =   coord[1]   * nty/pdims[1];
  jrange[1] = ((coord[1]+1)* nty/pdims[1]) + d;
  krange[0] =   coord[2]   * ntz/pdims[2];
  krange[1] = ((coord[2]+1)* ntz/pdims[2]) + d;

  /* -------------- */
  /* Set the "view" */
  /* -------------- */
  H5Block3dSetView(fldsfile, irange[0], irange[1],
                             jrange[0], jrange[1],
                             krange[0], krange[1]);

}

void H5input::ReadFields(double ***field, std::string fname, int nx, int ny, int nz, int rank){

  h5_float64_t*     buffer;

  buffer = new h5_float64_t[nz*ny*nx];

  H5Block3dReadScalarFieldFloat64(fldsfile, fname.c_str(), buffer);

  std::ofstream myfile;

  if (rank!=-1) {
    std::stringstream  ss;
    ss << "proc-" << rank << ".txt";
    std::string filename = ss.str();
    myfile.open(filename.c_str());
  }

  int n = 0;
  for (int k=1; k<nz-1; k++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=1; i<nx-1; i++) {
        field[i][j][k] = buffer[n++];
        if (rank!=-1) myfile << i << " " << j << " " << k << " : " << field[i][j][k] << std::endl;
      }
    }
  }
  if (rank!=-1) myfile.close();

  delete [] buffer;
}

void H5input::CloseFieldsFile(){
  H5CloseFile(fldsfile);
}


void H5input::OpenPartclFile(int ns){
  nspec = ns;
  part = new H5hutpart[nspec];
}

// void H5input::ReadParticles(int rank, int nproc, int *pdims, double *L, MPI_Comm CART_COMM){
//   int ndim = 3; // By default 2D is considered as a 'flat' 3D
//   LoadParticles(ndim, rank, nproc, pdims, L, CART_COMM);
//   //LoadParticles(ndim, rank, nproc, pdims, L, CART_COMM);
//   //InitParticles(ndim, rank, CART_COMM);
// }

void H5input::ClosePartclFile(){
  delete [] part;
}

void H5input::LoadLocalParticles(int ispec, long long r_nop, long long r_beg,
                                 double *r_q,
                                 double *r_x,
                                 double *r_y,
                                 double *r_z,
                                 double *r_u,
                                 double *r_v,
                                 double *r_w) {

  std::copy(r_q, r_q+r_nop, part[ispec].getQ()+r_beg);
  std::copy(r_x, r_x+r_nop, part[ispec].getX()+r_beg);
  std::copy(r_y, r_y+r_nop, part[ispec].getY()+r_beg);
  std::copy(r_z, r_z+r_nop, part[ispec].getZ()+r_beg);
  std::copy(r_u, r_u+r_nop, part[ispec].getU()+r_beg);
  std::copy(r_v, r_v+r_nop, part[ispec].getV()+r_beg);
  std::copy(r_w, r_w+r_nop, part[ispec].getW()+r_beg);

}

// void H5input::LoadLocalParticles(long long *np,
//                         double *q_loc,
//                         double *x_loc,
//                         double *y_loc,
//                         double *z_loc,
//                         double *u_loc,
//                         double *v_loc,
//                         double *w_loc) {
// 
//   long long n;
// 
//   n = 0;
//   for (int s=0; s<nspec; s++){
// 
//     part[s].memalloc(np[s]);
// 
//     std::copy(q_loc+n, q_loc+n+np[s], part[s].getQ());
//     std::copy(x_loc+n, x_loc+n+np[s], part[s].getX());
//     std::copy(y_loc+n, y_loc+n+np[s], part[s].getY());
//     std::copy(z_loc+n, z_loc+n+np[s], part[s].getZ());
//     std::copy(u_loc+n, u_loc+n+np[s], part[s].getU());
//     std::copy(v_loc+n, v_loc+n+np[s], part[s].getV());
//     std::copy(w_loc+n, w_loc+n+np[s], part[s].getW());
// 
//     n += np[s];
//   }
// }

void H5input::SortParticles(int nproc, int rank, int ispec, int ndim, long long nop, int *pdims, double *L, MPI_Comm CART_COMM,
                        double *q,
                        double *x, double *y, double *z,
                        double *u, double *v, double *w){

  unsigned int *inproc;
  int          *pcoord;
  int          irank;
  long long    *s_nop;
  long long    *r_nop;
  long long    *r_beg;
  long long    f_nop;
  long long    S_MAX_NOP;
  long long    R_MAX_NOP;

  double       *s_buffer;
  double       *r_buffer;

  pcoord = new int[ndim];
  inproc = new unsigned int[nop];
  s_nop  = new long long[nproc];
  r_nop  = new long long[nproc];
  r_beg  = new long long[nproc];

  MPI_Request req;
  MPI_Status  status;

  /* ------------------------- */
  /* Tag particles with proc # */
  /* ------------------------- */

  for (int iproc=0; iproc<nproc; iproc++) s_nop[iproc] = 0;

  for (long long p=0; p < nop; p++) {

    /* --------------------------------------- */
    /* Get processor coordinate for particle p */
    /* --------------------------------------- */

    for (int d=0; d<ndim; d++) {
      double r;
      if (d==0) r = x[p];
      if (d==1) r = y[p];
      if (d==2) r = z[p];

      pcoord[d] = int(r/(L[d]/pdims[d]));
    }

    /* ----------------------------------------------------- */
    /* Detect processor rank and increment particle counters */
    /* ----------------------------------------------------- */

    MPI_Cart_rank(CART_COMM, pcoord, &irank);
    inproc[p] = irank;
    s_nop[irank]++;

  }

  delete [] pcoord;

  /* ----------------------------------- */
  /* Communicate sizes and define offset */
  /* ----------------------------------- */

  S_MAX_NOP = 0;
  for (int iproc=0; iproc<nproc; iproc++) {
    if (rank!=iproc) MPI_Send(&s_nop[iproc], 1, MPI_LONG_LONG, iproc, 0, CART_COMM);
    if (s_nop[iproc] > S_MAX_NOP) S_MAX_NOP = s_nop[iproc];
  }

  f_nop  = 0;
  R_MAX_NOP = 0;
  for (int iproc=0; iproc<nproc; iproc++) {

    if   (rank==iproc) r_nop[iproc] = s_nop[iproc];
    else MPI_Recv(&r_nop[iproc], 1, MPI_LONG_LONG, iproc, 0, CART_COMM, MPI_STATUS_IGNORE);

    if (iproc==0)      r_beg[iproc]   = 0;
    if (iproc<nproc-1) r_beg[iproc+1] = r_beg[iproc] + r_nop[iproc];

    if (r_nop[iproc] > R_MAX_NOP) R_MAX_NOP = r_nop[iproc];

    // if (rank==0) std::cout << rank << "++ proc= " << iproc << " :: beg = " << r_beg[iproc] << " + len = " << r_nop[iproc] << std::endl;

    f_nop += r_nop[iproc];

  }

  /* ------------------------------------------- */
  /* Allocate particle memory in local processor */
  /* ------------------------------------------- */

  if (rank==0) {
    double memsize = 7*f_nop*sizeof(double)*1e-6;
    std::cout << "[PHDF5-io] Species " << ispec << " -- Alloc final vectors memory: " << memsize << " MB " << std::endl;
  }

  part[ispec].memalloc(f_nop);

  /* --------------------------------- */
  /* Allocate send and receive buffers */
  /* --------------------------------- */

  if (rank==0) {
    double s_memsize = 7*S_MAX_NOP*sizeof(double)*1e-6;
    double r_memsize = 7*R_MAX_NOP*sizeof(double)*1e-6;
    std::cout << "[PHDF5-io] Species " << ispec << " -- Alloc send/recv buffer memory: " << s_memsize+r_memsize << " MB " << std::endl;
  }

  s_buffer = new double[7*S_MAX_NOP];
  r_buffer = new double[7*R_MAX_NOP];

  /* --------------------------------- */
  /* Fill send and receive the buffers */
  /* --------------------------------- */

  for (int iproc=0; iproc<nproc; iproc++) {

    /* -------------------- */
    /* Fill the send buffer */
    /* -------------------- */

    long long p = 0;

    for (long long n=0; n<nop; n++) {
      if (inproc[n]==iproc) {
        s_buffer[p]                    = q[n];
        s_buffer[p + (1*s_nop[iproc])] = x[n];
        s_buffer[p + (2*s_nop[iproc])] = y[n];
        s_buffer[p + (3*s_nop[iproc])] = z[n];
        s_buffer[p + (4*s_nop[iproc])] = u[n];
        s_buffer[p + (5*s_nop[iproc])] = v[n];
        s_buffer[p + (6*s_nop[iproc])] = w[n];
        
        p++;
      }
    }

    if (p != s_nop[iproc]) {
      std::cout << rank << " ERROR: The number of particles dont match: " << std::endl;
      std::cout << rank << "        p = " << p << " :: s_nop = " << s_nop[iproc] << std::endl;
    }
    if (p > S_MAX_NOP) {
      std::cout << rank << " ERROR: The number of particles in the send buffer is bigger than the buffer size: " << std::endl;
      std::cout << rank << "        p = " << p << " :: S_MAX_NOP = " << S_MAX_NOP << std::endl;
    }

    /* -------------------- */
    /* Send buffer to iproc */
    /* -------------------- */

    if (iproc!=rank) {
      //std::cout << rank << " ::  Send buffer size 7x" << s_nop[iproc] << " to " << iproc << std::endl;
      MPI_Send(s_buffer, 7*s_nop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM);
      //MPI_Isend(s_buffer, 7*s_nop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
      //MPI_Wait(&req, &status);
    }

    /* ---------------------------------- */
    /* Receive buffer in iproc from jproc */
    /* ---------------------------------- */

    else {

      for (int jproc=0; jproc<nproc; jproc++) {
        if (jproc==rank) { //COPY FROM s_variables
          std::copy(s_buffer, s_buffer+(7*r_nop[jproc]), r_buffer);
        }
        else {
          MPI_Recv(r_buffer, 7*r_nop[jproc], MPI_DOUBLE, jproc, 0, CART_COMM, MPI_STATUS_IGNORE);
        }

        /* ------------------------------------------- */
        /* Add received buffer to the particle vectors */
        /* ------------------------------------------- */

        LoadLocalParticles(ispec, r_nop[jproc], r_beg[jproc],
                                  r_buffer, 
                                  r_buffer+(1*r_nop[jproc]),
                                  r_buffer+(2*r_nop[jproc]),
                                  r_buffer+(3*r_nop[jproc]),
                                  r_buffer+(4*r_nop[jproc]),
                                  r_buffer+(5*r_nop[jproc]),
                                  r_buffer+(6*r_nop[jproc]));

      }

    }

  }

  delete [] inproc;
  delete [] s_nop;
  delete [] r_nop;
  delete [] r_beg;

  delete [] s_buffer;
  delete [] r_buffer;

  if (rank==0) std::cout << "[PHDF5-io] Species " << ispec << " -- Temporary memory deallocated" << std::endl;

}

// void H5input::FindLocalParticles_par(int nproc, int rank, int ispec, int ndim, long long nop, int *pdims, double *L, MPI_Comm CART_COMM,
//                         double *q,
//                         double *x, double *y, double *z,
//                         double *u, double *v, double *w){
// 
//   double     *q_loc;
//   double     *x_loc;
//   double     *y_loc;
//   double     *z_loc;
//   double     *u_loc;
//   double     *v_loc;
//   double     *w_loc;
//   bool       *b;
//   long long  *locnop;
//   long long  *glonop;
//   long long  ntp;
//   long long  MAX_NOP;
//   long long  *globeg;
//   long long  np;
//   long long  p;
//   int        irank;
//   int        *pcoord;
//   int        *inproc;
// 
//   MPI_Request req   = MPI_REQUEST_NULL;
//   MPI_Status  status;
// 
//   pcoord = new int[ndim];
//   inproc = new int[nop];
//   locnop = new long long[nproc];
//   glonop = new long long[nproc];
//   globeg = new long long[nproc];
// 
//   /* ------------------- */
//   /* Initialize counters */
//   /* ------------------- */
// 
//   for (int iproc=0; iproc<nproc; iproc++) locnop[iproc] = 0;
// 
//   /* ------------------------- */
//   /* Tag particles with proc # */
//   /* ------------------------- */
// 
//   for (long long p=0; p < nop; p++) {
// 
//     /* --------------------------------------- */
//     /* Get processor coordinate for particle p */
//     /* --------------------------------------- */
// 
//     for (int d=0; d<ndim; d++) {
//       double r;
//       if (d==0) r = x[p];
//       if (d==1) r = y[p];
//       if (d==2) r = z[p];
// 
//       pcoord[d] = int(r/(L[d]/pdims[d]));
//     }
// 
//     /* ----------------------------------------------------- */
//     /* Detect processor rank and increment particle counters */
//     /* ----------------------------------------------------- */
// 
//     MPI_Cart_rank(CART_COMM, pcoord, &irank);
//     inproc[p] = irank;
//     locnop[irank]++;
// 
//   }
// 
//   /* ----------------------------------- */
//   /* Communicate sizes and define offset */
//   /* ----------------------------------- */
// 
//   MAX_NOP = 0;
//   for (int iproc=0; iproc<nproc; iproc++) {
//     if (rank!=iproc) MPI_Send(&locnop[iproc], 1, MPI_LONG_LONG, iproc, 0, CART_COMM);
//     if (locnop[iproc] > MAX_NOP) MAX_NOP = locnop[iproc];
//   }
// 
//   ntp  = 0;
//   for (int iproc=0; iproc<nproc; iproc++) {
// 
//     if   (rank==iproc) glonop[iproc] = locnop[iproc];
//     else MPI_Recv(&glonop[iproc], 1, MPI_LONG_LONG, iproc, 0, CART_COMM, MPI_STATUS_IGNORE);
// 
//     if (iproc==0)      globeg[iproc]   = 0;
//     if (iproc<nproc-1) globeg[iproc+1] = globeg[iproc] + glonop[iproc];
// 
//     if (rank==0) std::cout << " proc= " << iproc << " :: beg = " << globeg[iproc] << " + len = " << glonop[iproc] << std::endl;
// 
//     ntp += glonop[iproc];
//   }
// 
//   /* ------------------------------------------- */
//   /* Allocate particle memory in local processor */
//   /* ------------------------------------------- */
// 
//   std::cout << rank << "  alloc part ntp = " << ntp << std::endl;
//   part[ispec].memalloc(ntp);
// 
//   /* ----------------------------- */
//   /* Create the send buffer arrays */
//   /* ----------------------------- */
// 
//   q_loc = new double[MAX_NOP];
//   x_loc = new double[MAX_NOP];
//   y_loc = new double[MAX_NOP];
//   z_loc = new double[MAX_NOP];
//   u_loc = new double[MAX_NOP];
//   v_loc = new double[MAX_NOP];
//   w_loc = new double[MAX_NOP];
// 
//   /* ------------------------------------------ */
//   /* Fill the sending buffer for each processor */
//   /* ------------------------------------------ */
// 
//   for (int iproc=0; iproc<nproc; iproc++) {
// 
//     /* ---------------- */
//     /* Fill the buffers */
//     /* ---------------- */
// 
//     p = 0;
//     for (long long n=0; n<nop; n++) {
//       if (inproc[n]==iproc) {
//         q_loc[p] = q[n];
//         x_loc[p] = x[n];
//         y_loc[p] = y[n];
//         z_loc[p] = z[n];
//         u_loc[p] = u[n];
//         v_loc[p] = v[n];
//         w_loc[p] = w[n];
//         
//         p++;
//       }
//     }
// 
//     if (p != locnop[iproc]) {
//       std::cout << rank << " ERROR: The number of particles dont match: " << std::endl;
//       std::cout << rank << "        p = " << p << " :: loncnop = " << locnop[iproc] << std::endl;
//     }
// 
//     /* ----------------------------------------------------------------------- */
//     /* Copy particles in local procesor and send particles to other processors */
//     /* ----------------------------------------------------------------------- */
// 
//     //if (iproc==rank)  LoadLocalParticles_par(ispec, glonop[iproc], globeg[iproc], q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);
// 
//     else{
//       MPI_Isend(q_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(x_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(y_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(z_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(u_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(v_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Isend(w_loc, p, MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//     }
// 
//   }
// 
//   delete [] q_loc;
//   delete [] x_loc;
//   delete [] y_loc;
//   delete [] z_loc;
//   delete [] u_loc;
//   delete [] v_loc;
//   delete [] w_loc;
// 
//   /* ------------------------------------- */
//   /* Receive buffers from other processors */
//   /* ------------------------------------- */
// 
//   for (int iproc=0; iproc<nproc; iproc++) {
// 
//     if (iproc!=rank) {
// 
//       q_loc = new double[glonop[iproc]];
//       x_loc = new double[glonop[iproc]];
//       y_loc = new double[glonop[iproc]];
//       z_loc = new double[glonop[iproc]];
//       u_loc = new double[glonop[iproc]];
//       v_loc = new double[glonop[iproc]];
//       w_loc = new double[glonop[iproc]];
// 
//       std::cout << rank << "  Wating arrays from " << iproc << std::endl;
//       MPI_Irecv(q_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(x_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(y_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(z_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(u_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(v_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
//       //MPI_Irecv(w_loc , glonop[iproc], MPI_DOUBLE, iproc, 0, CART_COMM, &req);
// 
//       std::cout << rank << " :: beg, len, last, ntp : " << globeg[iproc] << ", " << glonop[iproc] << ", "<< globeg[iproc] + glonop[iproc] << ", " << ntp << std::endl;
//       //LoadLocalParticles_par(ispec, glonop[iproc], globeg[iproc], q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);
//       std::cout << rank << " ++ proc " << iproc << " included! " << std::endl;
// 
//       delete [] q_loc;
//       delete [] x_loc;
//       delete [] y_loc;
//       delete [] z_loc;
//       delete [] u_loc;
//       delete [] v_loc;
//       delete [] w_loc;
// 
//     }
// 
//   }
// 
//   std::cout << rank << " Sync comms now! " << std::endl;
//   int done;
//   MPI_Test(&req,&done,&status);
//   if (!done) std::cout << rank << " MPI error done = " << done << std::endl;
//   std::cout << rank << " All buffers received and loaded " << done << std::endl;
// 
//   delete [] pcoord;
//   delete [] inproc;
//   delete [] locnop;
// 
// }
// 
// void H5input::FindLocalParticles(int nproc, int ndim, h5_int64_t *h5npart, long long nop, int *pdims, double *L, MPI_Comm CART_COMM,
//                         double *q,
//                         double *x, double *y, double *z,
//                         double *u, double *v, double *w){
// 
//   double     *q_loc;
//   double     *x_loc;
//   double     *y_loc;
//   double     *z_loc;
//   double     *u_loc;
//   double     *v_loc;
//   double     *w_loc;
//   bool       *b;
//   long long  locnp[nspec][nproc];
//   long long  *locnop;
//   long long  *np;
//   long long  p;
//   long long  nn;
//   long long  MAX_NOP;
//   int        irank;
//   int        *pcoord;
//   int        *inproc;
// 
//   pcoord = new int[ndim];
//   inproc = new int[nop];
//   np     = new long long[nspec+1];
//   locnop = new long long[nproc];
// 
//   MAX_NOP = 1.1 * (nop/nproc);
// 
//   /* ------------------- */
//   /* Initialize counters */
//   /* ------------------- */
// 
//   for (int iproc=0; iproc<nproc; iproc++) {
//     locnop[iproc] = 0;
//     for (int i=0; i<nspec; i++) {
//       locnp[i][iproc] = 0;
//     }
//   }
// 
//   /* ------------------------- */
//   /* Tag particles with proc # */
//   /* ------------------------- */
// 
//   nn = 0;
//   for (int i=0; i<nspec; i++) {
//     for (long long n=0; n < h5npart[i]; n++) {
// 
//       p = nn + n;
// 
//       /* --------------------------------------- */
//       /* Get processor coordinate for particle p */
//       /* --------------------------------------- */
// 
//       for (int d=0; d<ndim; d++) {
//         double r;
//         if (d==0) r = x[p];
//         if (d==1) r = y[p];
//         if (d==2) r = z[p];
// 
//         pcoord[d] = int(r/(L[d]/pdims[d]));
//       }
// 
//       /* ----------------------------------------------------- */
//       /* Detect processor rank and increment particle counters */
//       /* ----------------------------------------------------- */
// 
//       MPI_Cart_rank(CART_COMM, pcoord, &irank);
//       inproc[p] = irank;
//       locnp[i][irank]++;
//       locnop[irank]++;
// 
//     }
// 
//     nn += h5npart[i];
//   }
// 
//   /* ----------------------- */
//   /* Create the local arrays */
//   /* ----------------------- */
// 
//   MAX_NOP = 0;
//   for (int iproc=0; iproc<nproc; iproc++){
//     if (locnop[iproc]>MAX_NOP) MAX_NOP = locnop[iproc];
//   }
// 
//   q_loc = new double[MAX_NOP];
//   x_loc = new double[MAX_NOP];
//   y_loc = new double[MAX_NOP];
//   z_loc = new double[MAX_NOP];
//   u_loc = new double[MAX_NOP];
//   v_loc = new double[MAX_NOP];
//   w_loc = new double[MAX_NOP];
// 
//   /* ------------------------------- */
//   /* Fill vectors for each processor */
//   /* ------------------------------- */
// 
//   for (int iproc=0; iproc<nproc; iproc++) {
// 
//     /* ---------------- */
//     /* Fill the vectors */
//     /* ---------------- */
// 
//     p = 0;
//     for (long long n=0; n<nop; n++) {
//       if (inproc[n]==iproc) {
//         q_loc[p] = q[n];
//         x_loc[p] = x[n];
//         y_loc[p] = y[n];
//         z_loc[p] = z[n];
//         u_loc[p] = u[n];
//         v_loc[p] = v[n];
//         w_loc[p] = w[n];
//         p++;
//       }
//     }
// 
//     /* ------------------------------------------------------------------- */
//     /* Copy particles in procesor 0 and send particles to other processors */
//     /* ------------------------------------------------------------------- */
// 
//     for (int i=0; i<nspec; i++) np[i] = locnp[i][iproc];
//     np[nspec] = locnop[iproc]; // The last element in np contains the total (all species combined) number of particles.
// 
//     if (iproc==0)  LoadLocalParticles(np, q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);
//     else{
//       MPI_Send(np    , nspec+1      , MPI_LONG_LONG, iproc, 0, CART_COMM);
//       MPI_Send(q_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(x_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(y_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(z_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(u_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(v_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//       MPI_Send(w_loc , locnop[iproc], MPI_DOUBLE,    iproc, 0, CART_COMM);
//     }
//   }
// 
//   delete [] q_loc;
//   delete [] x_loc;
//   delete [] y_loc;
//   delete [] z_loc;
//   delete [] u_loc;
//   delete [] v_loc;
//   delete [] w_loc;
// 
//   delete [] pcoord;
//   delete [] inproc;
//   delete [] np;
//   delete [] locnop;
// 
// }
// 
// void H5input::InitParticles(int ndim, int rank, MPI_Comm CART_COMM){
// 
//   long long  *np;
//   long long  nop;
//   double     *q_loc;
//   double     *x_loc;
//   double     *y_loc;
//   double     *z_loc;
//   double     *u_loc;
//   double     *v_loc;
//   double     *w_loc;
// 
//   if (rank!=0) {
//     /* ---------------------------------------- */
//     /* Receive particle information from proc 0 */
//     /* ---------------------------------------- */
// 
//     np = new long long[nspec+1];
// 
//     MPI_Recv(np, nspec+1, MPI_LONG_LONG, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     nop = np[nspec];
// 
//     /* ----------------------- */
//     /* Create the local arrays */
//     /* ----------------------- */
// 
//     q_loc = new double[nop];
//     x_loc = new double[nop];
//     y_loc = new double[nop];
//     z_loc = new double[nop];
//     u_loc = new double[nop];
//     v_loc = new double[nop];
//     w_loc = new double[nop];
// 
//     MPI_Recv(q_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(x_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(y_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(z_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(u_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(v_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
//     MPI_Recv(w_loc, nop, MPI_DOUBLE, 0, 0, CART_COMM, MPI_STATUS_IGNORE);
// 
//     LoadLocalParticles(np, q_loc, x_loc, y_loc, z_loc, u_loc, v_loc, w_loc);
// 
//     delete [] q_loc;
//     delete [] x_loc;
//     delete [] y_loc;
//     delete [] z_loc;
//     delete [] u_loc;
//     delete [] v_loc;
//     delete [] w_loc;
//   }
// 
// }

void H5input::ReadParticles(int rank, int nproc, int *pdims, double *L, MPI_Comm CART_COMM){

  int        ndim = 3; // By default 2D is considered as a 'flat' 3D
  int        h5nspec;
  long long  nop;
  long long  MAX_NOPS;
  long long  *nops_beg;
  long long  *nops_end;
  long long  *nops;
  h5_int64_t *h5npart;
  double     *i_q;
  double     *i_x;
  double     *i_y;
  double     *i_z;
  double     *i_u;
  double     *i_v;
  double     *i_w;

  herr_t status;
  hid_t  file;
  hid_t  dataset;
  hid_t  attr;
  hid_t  group;
  hid_t  fapl;

  std::stringstream sstm;

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(6) << recycle;
  filename = basename + "-Partcl" + "_" + filenmbr.str() + ".h5";

  fapl   = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(fapl, CART_COMM, MPI_INFO_NULL);
  file   = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  status = H5Pclose(fapl);

  /* --------------------------- */
  /* Compare Attributes to input */
  /* --------------------------- */

  group    = H5Gopen2(file,"/Step#0",H5P_DEFAULT);
  attr     = H5Aopen_name(group, "nspec");
  status   = H5Aread(attr, H5T_NATIVE_INT, &h5nspec);
  status   = H5Aclose(attr);

  if (nspec!=h5nspec) {
    std::cout << "[PHDF5-io]" << "ERROR in ReadParallelParticles: the number of species in the initial file" << std::endl;
    std::cout << "[PHDF5-io]" << "                                does not match the number of species requested." << std::endl;
    abort();
  }

  /* ----------------------------------------- */
  /* Set the reading limits for each processor */
  /* ----------------------------------------- */

  h5npart  = new h5_int64_t[nspec];
  nops_beg = new long long [nspec];
  nops_end = new long long [nspec];
  nops     = new long long [nspec];

  nop = 0;
  MAX_NOPS = 0;
  for (int i=0; i<nspec; i++) {
    sstm << "npart_" << i;
    std::string nparti = sstm.str();

    attr     = H5Aopen_name(group, nparti.c_str());
    status   = H5Aread(attr, H5T_NATIVE_INT, &h5npart[i]);
    status   = H5Aclose(attr);
    sstm.str("");

    nops_beg[i] = rank     * ceil(h5npart[i]/nproc);
    nops_end[i] = (rank+1) * ceil(h5npart[i]/nproc) - 1;
    if (rank==nproc-1) nops_end[i] = h5npart[i]-1;
    nops[i] = nops_end[i] - nops_beg[i] + 1;

    if (nops[i]>MAX_NOPS) MAX_NOPS = nops[i];
    nop += nops[i];
    //std::cout << rank << " (()) nops " << i << " = " << nops[i] << " [[]] beg-end = " << nops_beg[i] << " " << nops_end[i] << std::endl;
  }

  delete [] nops_end;

  /* --------------------------- */
  /* Allocate the reading arrays */
  /* --------------------------- */

  if (rank==0) {
    double memsize = 7*MAX_NOPS*sizeof(double)*1e-9;
    if (memsize<1.0) std::cout << "[PHDF5-io]" << " Allocating " << memsize*1e3 << " MB per processor to read particles file" << std::endl;
    else             std::cout << "[PHDF5-io]" << " Allocating " << memsize << " GB per processor to read particles file" << std::endl;
  }

  i_q = new double[MAX_NOPS];
  i_x = new double[MAX_NOPS];
  i_y = new double[MAX_NOPS];
  i_z = new double[MAX_NOPS];
  i_u = new double[MAX_NOPS];
  i_v = new double[MAX_NOPS];
  i_w = new double[MAX_NOPS];

  /* ------------------------------------- */
  /* Read the hyperslabs from the datasets */
  /* ------------------------------------- */

  long long offset = 0;

  for (int i=0; i<nspec; i++) {

    if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Found   " << h5npart[i] << " particles in file" << std::endl;
    if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Reading " << nops[i] << " particles in each processor" << std::endl;

    sstm << "q_" << i;
    std::string dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_q);
    sstm.str("");

    sstm << "x_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_x);
    sstm.str("");

    sstm << "y_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_y);
    sstm.str("");

    sstm << "z_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_z);
    sstm.str("");

    sstm << "u_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_u);
    sstm.str("");

    sstm << "v_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_v);
    sstm.str("");

    sstm << "w_" << i;
    dtset = sstm.str();
    ReadPartDataset(group, dtset, nops[i], nops_beg[i], i_w);
    sstm.str("");

    if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Sorting particles and distributing among processors: " << std::endl;

    SortParticles(nproc, rank, i, ndim, nops[i], pdims, L, CART_COMM, i_q, i_x, i_y, i_z, i_u, i_v, i_w);

  }

  if (rank==0) std::cout << "[PHDF5-io]" << " Freeing memory " << std::endl;

  delete [] i_q;
  delete [] i_x;
  delete [] i_y;
  delete [] i_z;
  delete [] i_u;
  delete [] i_v;
  delete [] i_w;

  delete [] h5npart;
  delete [] nops_beg;
  delete [] nops;
  
}

void H5input::ReadPartDataset(hid_t group, std::string dtset, long long nops, long long nops_beg, double *arr) {

  herr_t  status;
  hid_t   dsetid;
  hid_t   ds_mem;
  hid_t   ds_file;
  hid_t   dxpl;

  hsize_t* count  = new hsize_t[1];
  hsize_t* start  = new hsize_t[1];
 
  count [0] = nops;
  start [0] = nops_beg;

  dsetid  = H5Dopen(group, dtset.c_str(), H5P_DEFAULT);
  ds_file = H5Dget_space(dsetid);
  status  = H5Sselect_hyperslab(ds_file, H5S_SELECT_SET, start, NULL, count, NULL);
  ds_mem  = H5Screate_simple(1, count, NULL);
  dxpl    = H5Pcreate(H5P_DATASET_XFER);
  status  = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  status  = H5Dread(dsetid, H5T_NATIVE_DOUBLE, ds_mem, ds_file, dxpl, arr);

  status  = H5Pclose(dxpl);
  status  = H5Sclose(ds_file);
  status  = H5Sclose(ds_mem);
  status  = H5Dclose(dsetid);

  delete [] count;
  delete [] start;
}

// void H5input::LoadParticles(int ndim, int rank, int nproc, int *pdims, double *L, MPI_Comm CART_COMM){
// 
//   int        h5nspec;
//   int        ranks_rdr[1] = {0};
//   int        *ranks_rst;
//   long long  offset;
//   long long  nop;
//   h5_int64_t *h5npart;
//   double     *q;
//   double     *x;
//   double     *y;
//   double     *z;
//   double     *u;
//   double     *v;
//   double     *w;
// 
//   std::stringstream sstm;
// 
//   MPI_Comm  READER_COMM;
//   MPI_Group org_grp;
//   MPI_Group reader_grp;
// 
//   /* -------------------------------------------------- */
//   /* Create a new communicator for the reader processor */
//   /* This way only processor 0 reads the whole file.    */
//   /* It would be better if each processor read one part */
//   /* of the file...currently not possible with H5hut.   */
//   /* -------------------------------------------------- */
// 
//   MPI_Comm_group (CART_COMM, &org_grp);
//   ranks_rst = new int[nproc-1];
//   for (int i=0; i<nproc-1; i++)  ranks_rst[i] = i;
// 
//   if (rank==0) MPI_Group_incl (org_grp, 1, ranks_rdr, &reader_grp);
//   else         MPI_Group_incl (org_grp, nproc-1, ranks_rst, &reader_grp);
// 
//   delete [] ranks_rst;
// 
//   MPI_Comm_create(CART_COMM, reader_grp, &READER_COMM);
// 
//   if (rank==0) {
// 
//     /* -------------- */
//     /* Open HDF5 file */
//     /* -------------- */
// 
//     std::stringstream filenmbr;
//     std::string       filename;
// 
//     filenmbr << std::setfill('0') << std::setw(6) << recycle;
//     filename = basename + "-Partcl" + "_" + filenmbr.str() + ".h5";
// 
//     partfile = H5OpenFile(filename.c_str(), H5_O_RDONLY, READER_COMM);
//     H5SetStep(partfile, 0);
//     H5ReadStepAttribInt32(partfile, "nspec", &h5nspec);
// 
//     if (nspec!=h5nspec) {
//       std::cout << "[PHDF5-io]" << "ERROR in ReadParallelParticles: the number of species in the initial file" << std::endl;
//       std::cout << "[PHDF5-io]" << "                                does not match the number of species requested." << std::endl;
//       abort();
//     }
// 
//     h5npart = new h5_int64_t[nspec];
// 
//     nop = 0;
//     for (int i=0; i<nspec; i++) {
//       sstm << "npart_" << i;
//       std::string nparti = sstm.str();
//       H5ReadStepAttribInt64(partfile, nparti.c_str(), &h5npart[i]);
//       sstm.str("");
//       nop += h5npart[i];
//     }
// 
//     q = new double[nop];
//     x = new double[nop];
//     y = new double[nop];
//     z = new double[nop];
//     u = new double[nop];
//     v = new double[nop];
//     w = new double[nop];
// 
//     offset = 0;
//     for (int i=0; i<nspec; i++) {
// 
//       std::cout << "[PHDF5-io]" << " Reading a total of " << h5npart[i] << " particles of species " << i << std::endl;
// 
//       sstm << "q_" << i;
//       std::string dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),q+offset);
//       sstm.str("");
// 
//       sstm << "x_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),x+offset);
//       sstm.str("");
// 
//       sstm << "y_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),y+offset);
//       sstm.str("");
// 
//       sstm << "z_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),z+offset);
//       sstm.str("");
// 
//       sstm << "u_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),u+offset);
//       sstm.str("");
// 
//       sstm << "v_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),v+offset);
//       sstm.str("");
// 
//       sstm << "w_" << i;
//       dtset = sstm.str();
//       H5PartReadDataFloat64(partfile,dtset.c_str(),w+offset);
//       sstm.str("");
// 
//       offset += h5npart[i];
// 
//     }
// 
//     H5CloseFile(partfile);
// 
//     FindLocalParticles(nproc, ndim, h5npart, nop, pdims, L, CART_COMM, q, x, y, z, u, v, w);
// 
//     delete [] q;
//     delete [] x;
//     delete [] y;
//     delete [] z;
//     delete [] u;
//     delete [] v;
//     delete [] w;
// 
//     delete [] h5npart;
// 
//   }
// 
// }

void H5input::DumpPartclX(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getX(i);
  part[s].freeX();
}

void H5input::DumpPartclY(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getY(i);
  part[s].freeY();
}

void H5input::DumpPartclZ(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getZ(i);
  part[s].freeZ();
}

void H5input::DumpPartclU(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getU(i);
  part[s].freeU();
}

void H5input::DumpPartclV(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getV(i);
  part[s].freeV();
}

void H5input::DumpPartclW(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getW(i);
  part[s].freeW();
}

void H5input::DumpPartclQ(double*& tgt, int s){
  for (int i=0; i<part[s].getNP(); i++) tgt[i] = part[s].getQ(i);
  part[s].freeQ();
}
