
#include <fstream>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <new>

#include "H5hut-io.h"

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

  int n = 0;
  for (int k=1; k<nz-1; k++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=1; i<nx-1; i++) {
        buffer[n++] = field[i][j][k];
      }
    }
  }

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

  int n = 0;
  for (int k=1; k<nz-1; k++) {
    for (int j=1; j<ny-1; j++) {
      for (int i=1; i<nx-1; i++) {
        field[i][j][k] = buffer[n++];
      }
    }
  }

  delete [] buffer;
}

void H5input::CloseFieldsFile(){
  H5CloseFile(fldsfile);
}

void H5input::ClosePartclFile(){

  H5Gclose(pclgroup);
  H5Fclose(partfile);

  delete [] h5npart;
  delete [] nops_beg;
  delete [] nops;
  
}

void H5input::FillPartVectors(long long sizevec, int rank, int jproc, int ispec, long long r_nop, long long r_beg, double* r_buffer,
                              double *q, double *x, double *y, double *z, double *u, double *v, double *w) {

  long long l = 0;
  for (long long k=r_beg; k<r_beg+r_nop; k++){
    try {
      q[k] = r_buffer[l + 0*r_nop];
      x[k] = r_buffer[l + 1*r_nop];
      y[k] = r_buffer[l + 2*r_nop];
      z[k] = r_buffer[l + 3*r_nop];
      u[k] = r_buffer[l + 4*r_nop];
      v[k] = r_buffer[l + 5*r_nop];
      w[k] = r_buffer[l + 6*r_nop];
      l++;
    }
    catch (const std::out_of_range& e) {
      std::cout << rank << ": Out of Range error at element " << k << ": " << e.what() << std::endl;
    } 
  }

}

void H5input::SortParticles(int nproc, int rank, int ispec, int ndim, long long i_nop, int *pdims, double *L, MPI_Comm CART_COMM){

  int          *pcoord;
  int          irank;

  try {
    pcoord = new int[ndim];
    inproc = new unsigned int[i_nop];
    s_nop  = new long long[nproc];
    r_nop  = new long long[nproc];
    r_beg  = new long long[nproc];
  }
  catch(std::bad_alloc& exc)
  {
    std::cout << rank << " : Bad allocation of auxiliary vectors in SortParticles: " << exc.what() << std::endl;
  }

  MPI_Request req;
  MPI_Status  status;

  /* ------------------------- */
  /* Tag particles with proc # */
  /* ------------------------- */

  for (int iproc=0; iproc<nproc; iproc++) s_nop[iproc] = 0;

  for (long long p=0; p < i_nop; p++) {

    /* --------------------------------------- */
    /* Get processor coordinate for particle p */
    /* --------------------------------------- */

    for (int d=0; d<ndim; d++) {
      double r;
      if (d==0) r = i_x[p];
      if (d==1) r = i_y[p];
      if (d==2) r = i_z[p];

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

  long long ntp = 0;
  for (int iproc=0; iproc<nproc; iproc++) {
    if (s_nop[iproc]>0) {
      ntp+= s_nop[iproc];
    }
  }

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

    f_nop += r_nop[iproc];

  }

}

void H5input::ExchangeParticles(long long sizevec, int nproc, int rank, int ispec, long long i_nop, int *pdims, double *L, MPI_Comm CART_COMM,
                        double *q,
                        double *x, double *y, double *z,
                        double *u, double *v, double *w){

  double       *s_buffer;
  double       *r_buffer;

  /* --------------------------------- */
  /* Allocate send and receive buffers */
  /* --------------------------------- */

  if (rank==0) {
    double s_memsize = 7*S_MAX_NOP*sizeof(double)*1e-6;
    double r_memsize = 7*R_MAX_NOP*sizeof(double)*1e-6;
    std::cout << "[PHDF5-io] Species " << ispec << " -- Alloc send/recv buffer memory: " << s_memsize+r_memsize << " MB " << std::endl;
  }

  try{
    s_buffer = new double[7*S_MAX_NOP];
    r_buffer = new double[7*R_MAX_NOP];
  }
  catch(std::bad_alloc& exc)
  {
    std::cout << rank << " : Bad allocation of send/receive buffers: " << exc.what() << std::endl;
  }

  /* --------------------------------- */
  /* Fill send and receive the buffers */
  /* --------------------------------- */

  for (int iproc=0; iproc<nproc; iproc++) {

    /* -------------------- */
    /* Fill the send buffer */
    /* -------------------- */

    long long p = 0;

    for (long long n=0; n<i_nop; n++) {
      if (inproc[n]==iproc) {
        s_buffer[p]                    = i_q[n];
        s_buffer[p + (1*s_nop[iproc])] = i_x[n];
        s_buffer[p + (2*s_nop[iproc])] = i_y[n];
        s_buffer[p + (3*s_nop[iproc])] = i_z[n];
        s_buffer[p + (4*s_nop[iproc])] = i_u[n];
        s_buffer[p + (5*s_nop[iproc])] = i_v[n];
        s_buffer[p + (6*s_nop[iproc])] = i_w[n];
        
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

    /* ---------------------------------- */
    /* Receive buffer in iproc from jproc */
    /* ---------------------------------- */

    if (iproc==rank) {

      for (int jproc=0; jproc<nproc; jproc++) {

        if (r_nop[jproc]>0) {

          if (jproc==rank) {
            for (int k=0; k<7*r_nop[jproc]; k++) r_buffer[k] = s_buffer[k];
          }
          else {
            MPI_Recv(r_buffer, 7*r_nop[jproc], MPI_DOUBLE, jproc, iproc, CART_COMM, MPI_STATUS_IGNORE);
          }

          /* ------------------------------------------- */
          /* Add received buffer to the particle vectors */
          /* ------------------------------------------- */

          FillPartVectors(sizevec, rank, jproc, ispec, r_nop[jproc], r_beg[jproc], r_buffer, q, x, y, z, u, v, w);

        }
        // else{} // Nothing to comm

      }

    }

    /* -------------------- */
    /* Send buffer to iproc */
    /* -------------------- */

    else {
      if (s_nop[iproc]>0) {
        MPI_Send(s_buffer, 7*s_nop[iproc], MPI_DOUBLE, iproc, iproc, CART_COMM);
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

void H5input::OpenPartclFile(int ns, int rank, int nproc, MPI_Comm CART_COMM){

  long long         nops_end;
  int               h5nspec;
  herr_t            status;
  hid_t             attr;
  hid_t             fapl;
  std::stringstream sstm;

  nspec = ns;

  /* -------------- */
  /* Open HDF5 file */
  /* -------------- */

  std::stringstream filenmbr;
  std::string       filename;

  filenmbr << std::setfill('0') << std::setw(6) << recycle;
  filename = basename + "-Partcl" + "_" + filenmbr.str() + ".h5";

  fapl     = H5Pcreate(H5P_FILE_ACCESS);
  status   = H5Pset_fapl_mpio(fapl, CART_COMM, MPI_INFO_NULL);
  partfile = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  status = H5Pclose(fapl);

  /* --------------------------- */
  /* Compare Attributes to input */
  /* --------------------------- */

  pclgroup = H5Gopen2(partfile,"/Step#0",H5P_DEFAULT);
  attr     = H5Aopen_name(pclgroup, "nspec");
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

  h5npart  = new long long [nspec];
  nops_beg = new long long [nspec];
  nops     = new long long [nspec];

  for (int i=0; i<nspec; i++) {

    sstm << "npart_" << i;
    std::string nparti = sstm.str();

    attr     = H5Aopen_name(pclgroup, nparti.c_str());
    status   = H5Aread(attr, H5T_NATIVE_LLONG, &h5npart[i]);
    status   = H5Aclose(attr);
    sstm.str("");

    nops_beg[i] = rank * ceil(h5npart[i]/nproc);
    nops_end = (rank==nproc-1) ? h5npart[i] : (rank+1) * ceil(h5npart[i]/nproc);

    nops[i] = nops_end - nops_beg[i];

  }

}

void H5input::ReadParticles(int rank, int nproc, int i, int *pdims, double *L, MPI_Comm CART_COMM){

  int               ndim = 3; // By default 2D is considered as a 'flat' 3D
  long long         i_nop;
  long long         i_beg;
  std::stringstream sstm;

  /* --------------------------- */
  /* Allocate the reading arrays */
  /* --------------------------- */

  i_nop = nops    [i];
  i_beg = nops_beg[i];

  if (rank==0) {
    double memsize = 7*i_nop*sizeof(double)*1e-9;
    if (memsize<1.0) std::cout << "[PHDF5-io]" << " Allocating " << memsize*1e3 << " MB per processor to read " << i_nop << " particles from file" << std::endl;
    else             std::cout << "[PHDF5-io]" << " Allocating " << memsize << " GB per processor to read " << i_nop << " particles from file" << std::endl;
  }

  try{
    i_q = new double[i_nop];
    i_x = new double[i_nop];
    i_y = new double[i_nop];
    i_z = new double[i_nop];
    i_u = new double[i_nop];
    i_v = new double[i_nop];
    i_w = new double[i_nop];
  }
  catch(std::bad_alloc& exc)
  {
    std::cout << rank << " : Bad allocation of input vectors: " << exc.what() << std::endl;
  }

  /* ------------------------------------- */
  /* Read the hyperslabs from the datasets */
  /* ------------------------------------- */

  if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Found   " << h5npart[i] << " particles in file" << std::endl;
  if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Reading " << i_nop << " particles in each processor" << std::endl;

  sstm << "q_" << i;
  std::string dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_q);
  sstm.str("");

  sstm << "x_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_x);
  sstm.str("");

  sstm << "y_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_y);
  sstm.str("");

  sstm << "z_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_z);
  sstm.str("");

  sstm << "u_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_u);
  sstm.str("");

  sstm << "v_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_v);
  sstm.str("");

  sstm << "w_" << i;
  dtset = sstm.str();
  ReadPartDataset(pclgroup, dtset, i_nop, i_beg, i_w);
  sstm.str("");

  if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Sorting particles and distributing among processors: " << std::endl;

  SortParticles(nproc, rank, i, ndim, i_nop, pdims, L, CART_COMM);

}


void H5input::LoadParticles(long long sizevec, int rank, int nproc, int i, int *pdims, double *L, MPI_Comm CART_COMM,
                            double *q, double *x, double *y, double *z, 
                            double *u, double *v, double *w){

  ExchangeParticles(sizevec, nproc, rank, i, nops[i], pdims, L, CART_COMM, q, x, y, z, u, v, w);

  if (rank==0) std::cout << "[PHDF5-io] Species " << i << " Freeing memory " << std::endl;

  delete [] i_q;
  delete [] i_x;
  delete [] i_y;
  delete [] i_z;
  delete [] i_u;
  delete [] i_v;
  delete [] i_w;
  
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
