
#include <mpi.h>
#include "EMfields3D.h"
#include "Particles3Dcomm.h"
#include "TimeTasks.h"
#include "Moments.h"
#include "ompdefs.h"

/*! constructor */
EMfields3D::EMfields3D(Collective * col, Grid * grid) {
  nxc = grid->getNXC();
  nxn = grid->getNXN();
  nyc = grid->getNYC();
  nyn = grid->getNYN();
  nzc = grid->getNZC();
  nzn = grid->getNZN();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  invVOL = grid->getInvVOL();
  xStart = grid->getXstart();
  xEnd = grid->getXend();
  yStart = grid->getYstart();
  yEnd = grid->getYend();
  zStart = grid->getZstart();
  zEnd = grid->getZend();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  ns = col->getNs();
  c = col->getC();
  dt = col->getDt();
  th = col->getTh();
  ue0 = col->getU0(0);
  ve0 = col->getV0(0);
  we0 = col->getW0(0);
  x_center = col->getx_center();
  y_center = col->gety_center();
  z_center = col->getz_center();
  L_square = col->getL_square();

  delt = c * th * dt;
  PoissonCorrection = false;
  if (col->getPoissonCorrection()=="yes") PoissonCorrection = true;
  CGtol = col->getCGtol();
  GMREStol = col->getGMREStol();
  qom = new double[ns];
  for (int i = 0; i < ns; i++)
    qom[i] = col->getQOM(i);
  // boundary conditions: PHI and EM fields
  bcPHIfaceXright = col->getBcPHIfaceXright();
  bcPHIfaceXleft = col->getBcPHIfaceXleft();
  bcPHIfaceYright = col->getBcPHIfaceYright();
  bcPHIfaceYleft = col->getBcPHIfaceYleft();
  bcPHIfaceZright = col->getBcPHIfaceZright();
  bcPHIfaceZleft = col->getBcPHIfaceZleft();

  bcEMfaceXright = col->getBcEMfaceXright();
  bcEMfaceXleft = col->getBcEMfaceXleft();
  bcEMfaceYright = col->getBcEMfaceYright();
  bcEMfaceYleft = col->getBcEMfaceYleft();
  bcEMfaceZright = col->getBcEMfaceZright();
  bcEMfaceZleft = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  // Earth Simulation
  B1x = col->getB1x();
  B1y = col->getB1y();
  B1z = col->getB1z();
  delta = col->getDelta();
  Smooth = col->getSmooth();
  // get the density background for the gem Challange
  rhoINIT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i = 0; i < ns; i++) {
    rhoINIT[i] = col->getRHOinit(i);
    if ((fabs(col->getW0(i)) != 0) || (fabs(col->getU0(i)) != 0)) // GEM and LHDI
      DriftSpecies[i] = true;
    else
      DriftSpecies[i] = false;
  }
  /*! parameters for GEM challenge */
  FourPI = 16 * atan(1.0);
  /*! Restart */
  restart1 = col->getRestart_status();
  RestartDirName = col->getRestartDirName();
  Case = col->getCase();

  // OpenBC
  injFieldsLeft   = new injInfoFields(nxn, nyn, nzn);
  injFieldsRight  = new injInfoFields(nxn, nyn, nzn);
  injFieldsTop    = new injInfoFields(nxn, nyn, nzn);
  injFieldsBottom = new injInfoFields(nxn, nyn, nzn);
  injFieldsFront  = new injInfoFields(nxn, nyn, nzn);
  injFieldsRear   = new injInfoFields(nxn, nyn, nzn);

  // arrays allocation: nodes
  Ex = newArr3(double, nxn, nyn, nzn);
  Ey = newArr3(double, nxn, nyn, nzn);
  Ez = newArr3(double, nxn, nyn, nzn);
  Exth = newArr3(double, nxn, nyn, nzn);
  Eyth = newArr3(double, nxn, nyn, nzn);
  Ezth = newArr3(double, nxn, nyn, nzn);
  Bxn = newArr3(double, nxn, nyn, nzn);
  Byn = newArr3(double, nxn, nyn, nzn);
  Bzn = newArr3(double, nxn, nyn, nzn);
  rhon = newArr3(double, nxn, nyn, nzn);
  Jx = newArr3(double, nxn, nyn, nzn);
  Jy = newArr3(double, nxn, nyn, nzn);
  Jz = newArr3(double, nxn, nyn, nzn);
  Jxh = newArr3(double, nxn, nyn, nzn);
  Jyh = newArr3(double, nxn, nyn, nzn);
  Jzh = newArr3(double, nxn, nyn, nzn);
  // External imposed fields
  Bx_ext = newArr3(double,nxn,nyn,nzn);
  By_ext = newArr3(double,nxn,nyn,nzn);
  Bz_ext = newArr3(double,nxn,nyn,nzn);
  Jx_ext = newArr3(double,nxn,nyn,nzn);
  Jy_ext = newArr3(double,nxn,nyn,nzn);
  Jz_ext = newArr3(double,nxn,nyn,nzn);
  // involving species
  rhons = newArr4(double, ns, nxn, nyn, nzn);
  rhocs = newArr4(double, ns, nxc, nyc, nzc);
  Jxs = newArr4(double, ns, nxn, nyn, nzn);
  Jys = newArr4(double, ns, nxn, nyn, nzn);
  Jzs = newArr4(double, ns, nxn, nyn, nzn);
  pXXsn = newArr4(double, ns, nxn, nyn, nzn);
  pXYsn = newArr4(double, ns, nxn, nyn, nzn);
  pXZsn = newArr4(double, ns, nxn, nyn, nzn);
  pYYsn = newArr4(double, ns, nxn, nyn, nzn);
  pYZsn = newArr4(double, ns, nxn, nyn, nzn);
  pZZsn = newArr4(double, ns, nxn, nyn, nzn);
  // arrays allocation: central points 
  PHI = newArr3(double, nxc, nyc, nzc);
  Bxc = newArr3(double, nxc, nyc, nzc);
  Byc = newArr3(double, nxc, nyc, nzc);
  Bzc = newArr3(double, nxc, nyc, nzc);
  rhoc = newArr3(double, nxc, nyc, nzc);
  rhoh = newArr3(double, nxc, nyc, nzc);

  // temporary arrays
  tempXC = newArr3(double, nxc, nyc, nzc);
  tempYC = newArr3(double, nxc, nyc, nzc);
  tempZC = newArr3(double, nxc, nyc, nzc);

  tempXN = newArr3(double, nxn, nyn, nzn);
  tempYN = newArr3(double, nxn, nyn, nzn);
  tempZN = newArr3(double, nxn, nyn, nzn);
  tempC = newArr3(double, nxc, nyc, nzc);
  tempX = newArr3(double, nxn, nyn, nzn);
  tempY = newArr3(double, nxn, nyn, nzn);
  tempZ = newArr3(double, nxn, nyn, nzn);
  temp2X = newArr3(double, nxn, nyn, nzn);
  temp2Y = newArr3(double, nxn, nyn, nzn);
  temp2Z = newArr3(double, nxn, nyn, nzn);
  imageX = newArr3(double, nxn, nyn, nzn);
  imageY = newArr3(double, nxn, nyn, nzn);
  imageZ = newArr3(double, nxn, nyn, nzn);
  Dx = newArr3(double, nxn, nyn, nzn);
  Dy = newArr3(double, nxn, nyn, nzn);
  Dz = newArr3(double, nxn, nyn, nzn);
  vectX = newArr3(double, nxn, nyn, nzn);
  vectY = newArr3(double, nxn, nyn, nzn);
  vectZ = newArr3(double, nxn, nyn, nzn);
  divC = newArr3(double, nxc, nyc, nzc);
  sizeMomentsArray = omp_thread_count();
  momentsArray = new Moments*[sizeMomentsArray];
  for(int i=0;i<sizeMomentsArray;i++)
  {
    momentsArray[i] = new Moments;
    momentsArray[i]->init(nxn,nyn,nzn);
  }
}

// This was Particles3Dcomm::interpP2G()
void EMfields3D::sumMoments(const Particles3Dcomm& pcls, Grid * grid, VirtualTopology3D * vct)
{
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  double const*const x = pcls.getXall();
  double const*const y = pcls.getYall();
  double const*const z = pcls.getZall();
  double const*const u = pcls.getUall();
  double const*const v = pcls.getVall();
  double const*const w = pcls.getWall();
  double const*const q = pcls.getQall();
  //
  const int is = pcls.get_ns();
  double* rhons1d = rhons[is][0][0];
  double* Jxs1d   = Jxs  [is][0][0];
  double* Jys1d   = Jys  [is][0][0];
  double* Jzs1d   = Jzs  [is][0][0];
  double* pXXsn1d = pXXsn[is][0][0];
  double* pXYsn1d = pXYsn[is][0][0];
  double* pXZsn1d = pXZsn[is][0][0];
  double* pYYsn1d = pYYsn[is][0][0];
  double* pYZsn1d = pYZsn[is][0][0];
  double* pZZsn1d = pZZsn[is][0][0];
  //
  const long long nop_ll = pcls.getNOP();
  const int nop = pcls.getNOP();
  assert_le(nop_ll,INT_MAX); // else would need to use long long
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP
  #pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    Moments& speciesMoments = fetch_momentsArray(thread_num);
    speciesMoments.set_to_zero();
    //
    double*** rho = speciesMoments.fetch_rho();
    double*** Jx  = speciesMoments.fetch_Jx();
    double*** Jy  = speciesMoments.fetch_Jy();
    double*** Jz  = speciesMoments.fetch_Jz();
    double*** Pxx = speciesMoments.fetch_Pxx();
    double*** Pxy = speciesMoments.fetch_Pxy();
    double*** Pxz = speciesMoments.fetch_Pxz();
    double*** Pyy = speciesMoments.fetch_Pyy();
    double*** Pyz = speciesMoments.fetch_Pyz();
    double*** Pzz = speciesMoments.fetch_Pzz();
    // The following loop is expensive, so it is wise to assume that the
    // compiler is stupid.  Therefore we should on the one hand
    // expand things out and on the other hand avoid repeating computations.
    #pragma omp for
    for (int i = 0; i < nop; i++)
    {
      // compute the quadratic moments of velocity
      //
      const double ui=u[i];
      const double vi=v[i];
      const double wi=w[i];
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      //const double xi0   = x[i] - grid->getXN(ix - 1, iy, iz);
      //const double eta0  = y[i] - grid->getYN(ix, iy - 1, iz);
      //const double zeta0 = z[i] - grid->getZN(ix, iy, iz - 1);
      //const double xi1   = grid->getXN(ix, iy, iz) - x[i];
      //const double eta1  = grid->getYN(ix, iy, iz) - y[i];
      //const double zeta1 = grid->getZN(ix, iy, iz) - z[i];
      const double xi0   = x[i] - grid->getXN(ix-1);
      const double eta0  = y[i] - grid->getYN(iy-1);
      const double zeta0 = z[i] - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - x[i];
      const double eta1  = grid->getYN(iy) - y[i];
      const double zeta1 = grid->getZN(iz) - z[i];
      const double qi = q[i];
      const double weight000 = qi * xi0 * eta0 * zeta0 * invVOL;
      const double weight001 = qi * xi0 * eta0 * zeta1 * invVOL;
      const double weight010 = qi * xi0 * eta1 * zeta0 * invVOL;
      const double weight011 = qi * xi0 * eta1 * zeta1 * invVOL;
      const double weight100 = qi * xi1 * eta0 * zeta0 * invVOL;
      const double weight101 = qi * xi1 * eta0 * zeta1 * invVOL;
      const double weight110 = qi * xi1 * eta1 * zeta0 * invVOL;
      const double weight111 = qi * xi1 * eta1 * zeta1 * invVOL;
      //
      // use the weight to distribute the moments
      //
      // add charge density
      //speciesMoments.addRho(weight, ix, iy, iz);
      rho[ix  ][iy  ][iz  ] += weight000;
      rho[ix  ][iy  ][iz-1] += weight001;
      rho[ix  ][iy-1][iz  ] += weight010;
      rho[ix  ][iy-1][iz-1] += weight011;
      rho[ix-1][iy  ][iz  ] += weight100;
      rho[ix-1][iy  ][iz-1] += weight101;
      rho[ix-1][iy-1][iz  ] += weight110;
      rho[ix-1][iy-1][iz-1] += weight111;
      // add current density - X
      //speciesMoments.addJx(temp, ix, iy, iz);
      Jx[ix  ][iy  ][iz  ] += ui*weight000;
      Jx[ix  ][iy  ][iz-1] += ui*weight001;
      Jx[ix  ][iy-1][iz  ] += ui*weight010;
      Jx[ix  ][iy-1][iz-1] += ui*weight011;
      Jx[ix-1][iy  ][iz  ] += ui*weight100;
      Jx[ix-1][iy  ][iz-1] += ui*weight101;
      Jx[ix-1][iy-1][iz  ] += ui*weight110;
      Jx[ix-1][iy-1][iz-1] += ui*weight111;
      // add current density - Y
      //speciesMoments.addJy(temp, ix, iy, iz);
      Jy[ix  ][iy  ][iz  ] += vi*weight000;
      Jy[ix  ][iy  ][iz-1] += vi*weight001;
      Jy[ix  ][iy-1][iz  ] += vi*weight010;
      Jy[ix  ][iy-1][iz-1] += vi*weight011;
      Jy[ix-1][iy  ][iz  ] += vi*weight100;
      Jy[ix-1][iy  ][iz-1] += vi*weight101;
      Jy[ix-1][iy-1][iz  ] += vi*weight110;
      Jy[ix-1][iy-1][iz-1] += vi*weight111;
      // add current density - Z
      //speciesMoments.addJz(temp, ix, iy, iz);
      Jz[ix  ][iy  ][iz  ] += wi*weight000;
      Jz[ix  ][iy  ][iz-1] += wi*weight001;
      Jz[ix  ][iy-1][iz  ] += wi*weight010;
      Jz[ix  ][iy-1][iz-1] += wi*weight011;
      Jz[ix-1][iy  ][iz  ] += wi*weight100;
      Jz[ix-1][iy  ][iz-1] += wi*weight101;
      Jz[ix-1][iy-1][iz  ] += wi*weight110;
      Jz[ix-1][iy-1][iz-1] += wi*weight111;
      // Pxx - add pressure tensor
      //speciesMoments.addPxx(temp, ix, iy, iz);
      Pxx[ix  ][iy  ][iz  ] += uui*weight000;
      Pxx[ix  ][iy  ][iz-1] += uui*weight001;
      Pxx[ix  ][iy-1][iz  ] += uui*weight010;
      Pxx[ix  ][iy-1][iz-1] += uui*weight011;
      Pxx[ix-1][iy  ][iz  ] += uui*weight100;
      Pxx[ix-1][iy  ][iz-1] += uui*weight101;
      Pxx[ix-1][iy-1][iz  ] += uui*weight110;
      Pxx[ix-1][iy-1][iz-1] += uui*weight111;
      // Pxy - add pressure tensor
      //speciesMoments.addPxy(temp, ix, iy, iz);
      Pxy[ix  ][iy  ][iz  ] += uvi*weight000;
      Pxy[ix  ][iy  ][iz-1] += uvi*weight001;
      Pxy[ix  ][iy-1][iz  ] += uvi*weight010;
      Pxy[ix  ][iy-1][iz-1] += uvi*weight011;
      Pxy[ix-1][iy  ][iz  ] += uvi*weight100;
      Pxy[ix-1][iy  ][iz-1] += uvi*weight101;
      Pxy[ix-1][iy-1][iz  ] += uvi*weight110;
      Pxy[ix-1][iy-1][iz-1] += uvi*weight111;
      // Pxz - add pressure tensor
      //speciesMoments.addPxz(temp, ix, iy, iz);
      Pxz[ix  ][iy  ][iz  ] += uwi*weight000;
      Pxz[ix  ][iy  ][iz-1] += uwi*weight001;
      Pxz[ix  ][iy-1][iz  ] += uwi*weight010;
      Pxz[ix  ][iy-1][iz-1] += uwi*weight011;
      Pxz[ix-1][iy  ][iz  ] += uwi*weight100;
      Pxz[ix-1][iy  ][iz-1] += uwi*weight101;
      Pxz[ix-1][iy-1][iz  ] += uwi*weight110;
      Pxz[ix-1][iy-1][iz-1] += uwi*weight111;
      // Pyy - add pressure tensor
      //speciesMoments.addPyy(temp, ix, iy, iz);
      Pyy[ix  ][iy  ][iz  ] += vvi*weight000;
      Pyy[ix  ][iy  ][iz-1] += vvi*weight001;
      Pyy[ix  ][iy-1][iz  ] += vvi*weight010;
      Pyy[ix  ][iy-1][iz-1] += vvi*weight011;
      Pyy[ix-1][iy  ][iz  ] += vvi*weight100;
      Pyy[ix-1][iy  ][iz-1] += vvi*weight101;
      Pyy[ix-1][iy-1][iz  ] += vvi*weight110;
      Pyy[ix-1][iy-1][iz-1] += vvi*weight111;
      // Pyz - add pressure tensor
      //speciesMoments.addPyz(temp, ix, iy, iz);
      Pyz[ix  ][iy  ][iz  ] += vwi*weight000;
      Pyz[ix  ][iy  ][iz-1] += vwi*weight001;
      Pyz[ix  ][iy-1][iz  ] += vwi*weight010;
      Pyz[ix  ][iy-1][iz-1] += vwi*weight011;
      Pyz[ix-1][iy  ][iz  ] += vwi*weight100;
      Pyz[ix-1][iy  ][iz-1] += vwi*weight101;
      Pyz[ix-1][iy-1][iz  ] += vwi*weight110;
      Pyz[ix-1][iy-1][iz-1] += vwi*weight111;
      // Pzz - add pressure tensor
      //speciesMoments.addPzz(temp, ix, iy, iz);
      Pzz[ix  ][iy  ][iz  ] += wwi*weight000;
      Pzz[ix  ][iy  ][iz-1] += wwi*weight001;
      Pzz[ix  ][iy-1][iz  ] += wwi*weight010;
      Pzz[ix  ][iy-1][iz-1] += wwi*weight011;
      Pzz[ix-1][iy  ][iz  ] += wwi*weight100;
      Pzz[ix-1][iy  ][iz-1] += wwi*weight101;
      Pzz[ix-1][iy-1][iz  ] += wwi*weight110;
      Pzz[ix-1][iy-1][iz-1] += wwi*weight111;
    }
    // The following way is an equivalent reduction but less
    // efficient for a large number of threads.
    //
    //#pragma omp critical
    //addToSpeciesMoments(speciesMoments,is);
    //
    // Instead we split up the reduction tasks.
    //
    // One-dimensional array access is presumably
    // more efficient on poor compilers.
    //
    const double*const rho1d = rho[0][0];
    const double*const Jx1d  = Jx [0][0];
    const double*const Jy1d  = Jy [0][0];
    const double*const Jz1d  = Jz [0][0];
    const double*const Pxx1d = Pxx[0][0];
    const double*const Pxy1d = Pxy[0][0];
    const double*const Pxz1d = Pxz[0][0];
    const double*const Pyy1d = Pyy[0][0];
    const double*const Pyz1d = Pyz[0][0];
    const double*const Pzz1d = Pzz[0][0];
    //
    assert_eq(speciesMoments.get_nx(), nxn);
    assert_eq(speciesMoments.get_ny(), nyn);
    assert_eq(speciesMoments.get_nz(), nzn);
    const int numel = nxn*nyn*nzn;
    #pragma omp critical
    for(int i=0;i<numel;i++) rhons1d[i] += invVOL*rho1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) Jxs1d  [i] += invVOL*Jx1d [i];
    #pragma omp critical
    for(int i=0;i<numel;i++) Jys1d  [i] += invVOL*Jy1d [i];
    #pragma omp critical
    for(int i=0;i<numel;i++) Jzs1d  [i] += invVOL*Jz1d [i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pXXsn1d[i] += invVOL*Pxx1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pXYsn1d[i] += invVOL*Pxy1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pXZsn1d[i] += invVOL*Pxz1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pYYsn1d[i] += invVOL*Pyy1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pYZsn1d[i] += invVOL*Pyz1d[i];
    #pragma omp critical
    for(int i=0;i<numel;i++) pZZsn1d[i] += invVOL*Pzz1d[i];
  }
  communicateGhostP2G(is, 0, 0, 0, 0, vct);
}

/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;
  double ***divE = newArr3(double, nxc, nyc, nzc);
  double ***gradPHIX = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIY = newArr3(double, nxn, nyn, nzn);
  double ***gradPHIZ = newArr3(double, nxn, nyn, nzn);

  double *xkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 E components
  double *bkrylov = new double[3 * (nxn - 2) * (nyn - 2) * (nzn - 2)];  // 3 components

  double *xkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  double *bkrylovPoisson = new double[(nxc - 2) * (nyc - 2) * (nzc - 2)];
  // set to zero all the stuff 
  eqValue(0.0, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, bkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, divE, nxc, nyc, nzc);
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, gradPHIX, nxn, nyn, nzn);
  eqValue(0.0, gradPHIY, nxn, nyn, nzn);
  eqValue(0.0, gradPHIZ, nxn, nyn, nzn);
  // Adjust E calculating laplacian(PHI) = div(E) -4*PI*rho DIVERGENCE CLEANING
  if (PoissonCorrection) {
    if (vct->getCartesian_rank() == 0)
      cout << "*** DIVERGENCE CLEANING ***" << endl;
    grid->divN2C(divE, Ex, Ey, Ez);
    scale(tempC, rhoc, -FourPI, nxc, nyc, nzc);
    sum(divE, tempC, nxc, nyc, nzc);
    // move to krylov space 
    phys2solver(bkrylovPoisson, divE, nxc, nyc, nzc);
    // use conjugate gradient first
    if (!CG(xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 3000, CGtol, &Field::PoissonImage, grid, vct, this)) {
      if (vct->getCartesian_rank() == 0)
        cout << "CG not Converged. Trying with GMRes. Consider to increase the number of the CG iterations" << endl;
      eqValue(0.0, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2));
      GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc - 2) * (nyc - 2) * (nzc - 2), bkrylovPoisson, 20, 200, GMREStol, grid, vct, this);
    }
    solver2phys(PHI, xkrylovPoisson, nxc, nyc, nzc);
    communicateCenterBC(nxc, nyc, nzc, PHI, 2, 2, 2, 2, 2, 2, vct);
    // calculate the gradient
    grid->gradC2N(gradPHIX, gradPHIY, gradPHIZ, PHI);
    // sub
    sub(Ex, gradPHIX, nxn, nyn, nzn);
    sub(Ey, gradPHIY, nxn, nyn, nzn);
    sub(Ez, gradPHIZ, nxn, nyn, nzn);

  }                             // end of divergence cleaning 
  if (vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER ***" << endl;
  // prepare the source 
  MaxwellSource(bkrylov, grid, vct, col);
  phys2solver(xkrylov, Ex, Ey, Ez, nxn, nyn, nzn);
  // solver
  GMRES(&Field::MaxwellImage, xkrylov, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2), bkrylov, 20, 200, GMREStol, grid, vct, this);
  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, xkrylov, nxn, nyn, nzn);

  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);

  // apply to smooth to electric field 3 times
  smoothE(Smooth, vct, col);
  smoothE(Smooth, vct, col);
  smoothE(Smooth, vct, col);

  // communicate so the interpolation can have good values
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

  // OpenBC
  BoundaryConditionsE(Exth, Eyth, Ezth, nxn, nyn, nzn, grid, vct);
  BoundaryConditionsE(Ex, Ey, Ez, nxn, nyn, nzn, grid, vct);

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;
  delete[]xkrylovPoisson;
  delete[]bkrylovPoisson;
  delArr3(divE, nxc, nyc);
  delArr3(gradPHIX, nxn, nyn);
  delArr3(gradPHIY, nxn, nyn);
  delArr3(gradPHIZ, nxn, nyn);

}

/*! Calculate sorgent for Maxwell solver */
void EMfields3D::MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, tempXN, nxn, nyn, nzn);
  eqValue(0.0, tempYN, nxn, nyn, nzn);
  eqValue(0.0, tempZN, nxn, nyn, nzn);
  eqValue(0.0, temp2X, nxn, nyn, nzn);
  eqValue(0.0, temp2Y, nxn, nyn, nzn);
  eqValue(0.0, temp2Z, nxn, nyn, nzn);
  // communicate
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  if (Case=="ForceFree") fixBforcefree(grid,vct);
  if (Case=="GEM")       fixBgem(grid, vct);
  if (Case=="GEMnoPert") fixBgem(grid, vct);

  // OpenBC:
  BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc,grid,vct);

  // prepare curl of B for known term of Maxwell solver: for the source term
  grid->curlC2N(tempXN, tempYN, tempZN, Bxc, Byc, Bzc);
  scale(temp2X, Jxh, -FourPI / c, nxn, nyn, nzn);
  scale(temp2Y, Jyh, -FourPI / c, nxn, nyn, nzn);
  scale(temp2Z, Jzh, -FourPI / c, nxn, nyn, nzn);

  // -- dipole SOURCE version using J_ext
  addscale(-FourPI/c,temp2X,Jx_ext,nxn,nyn,nzn);
  addscale(-FourPI/c,temp2Y,Jy_ext,nxn,nyn,nzn);
  addscale(-FourPI/c,temp2Z,Jz_ext,nxn,nyn,nzn);
  // -- end of dipole SOURCE version using J_ext

  sum(temp2X, tempXN, nxn, nyn, nzn);
  sum(temp2Y, tempYN, nxn, nyn, nzn);
  sum(temp2Z, tempZN, nxn, nyn, nzn);
  scale(temp2X, delt, nxn, nyn, nzn);
  scale(temp2Y, delt, nxn, nyn, nzn);
  scale(temp2Z, delt, nxn, nyn, nzn);

  communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct);
  grid->gradC2N(tempX, tempY, tempZ, rhoh);

  scale(tempX, -delt * delt * FourPI, nxn, nyn, nzn);
  scale(tempY, -delt * delt * FourPI, nxn, nyn, nzn);
  scale(tempZ, -delt * delt * FourPI, nxn, nyn, nzn);
  // sum E, past values
  sum(tempX, Ex, nxn, nyn, nzn);
  sum(tempY, Ey, nxn, nyn, nzn);
  sum(tempZ, Ez, nxn, nyn, nzn);
  // sum curl(B) + jhat part
  sum(tempX, temp2X, nxn, nyn, nzn);
  sum(tempY, temp2Y, nxn, nyn, nzn);
  sum(tempZ, temp2Z, nxn, nyn, nzn);

  // Boundary condition in the known term
  // boundary condition: Xleft
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 0);
  // boundary condition: Xright
  if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 0);
  // boundary condition: Yleft
  if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 1);
  // boundary condition: Yright
  if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 1);
  // boundary condition: Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 2);
  // boundary condition: Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 2);

  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, nxn, nyn, nzn);

}
/*! Mapping of Maxwell image to give to solver */
void EMfields3D::MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct) {
  eqValue(0.0, im, 3 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  eqValue(0.0, imageX, nxn, nyn, nzn);
  eqValue(0.0, imageY, nxn, nyn, nzn);
  eqValue(0.0, imageZ, nxn, nyn, nzn);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, Dx, nxn, nyn, nzn);
  eqValue(0.0, Dy, nxn, nyn, nzn);
  eqValue(0.0, Dz, nxn, nyn, nzn);
  // move from krylov space to physical space
  solver2phys(vectX, vectY, vectZ, vector, nxn, nyn, nzn);
  grid->lapN2N(imageX, vectX, vct);
  grid->lapN2N(imageY, vectY, vct);
  grid->lapN2N(imageZ, vectZ, vct);
  neg(imageX, nxn, nyn, nzn);
  neg(imageY, nxn, nyn, nzn);
  neg(imageZ, nxn, nyn, nzn);
  // grad(div(mu dot E(n + theta)) mu dot E(n + theta) = D
  MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, grid);
  grid->divN2C(divC, Dx, Dy, Dz);
  // communicate you should put BC 
  // think about the Physics 
  // communicateCenterBC(nxc,nyc,nzc,divC,1,1,1,1,1,1,vct);
  communicateCenterBC(nxc, nyc, nzc, divC, 2, 2, 2, 2, 2, 2, vct);  // GO with Neumann, now then go with rho

  grid->gradC2N(tempX, tempY, tempZ, divC);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta))
  sub(imageX, tempX, nxn, nyn, nzn);
  sub(imageY, tempY, nxn, nyn, nzn);
  sub(imageZ, tempZ, nxn, nyn, nzn);

  // scale delt*delt
  scale(imageX, delt * delt, nxn, nyn, nzn);
  scale(imageY, delt * delt, nxn, nyn, nzn);
  scale(imageZ, delt * delt, nxn, nyn, nzn);

  // -lap(E(n +theta)) - grad(div(mu dot E(n + theta)) + eps dot E(n + theta)
  sum(imageX, Dx, nxn, nyn, nzn);
  sum(imageY, Dy, nxn, nyn, nzn);
  sum(imageZ, Dz, nxn, nyn, nzn);
  sum(imageX, vectX, nxn, nyn, nzn);
  sum(imageY, vectY, nxn, nyn, nzn);
  sum(imageZ, vectZ, nxn, nyn, nzn);

  // boundary condition: Xleft
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 0, grid);
  // boundary condition: Xright
  if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 0, grid);
  // boundary condition: Yleft
  if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 1, grid);
  // boundary condition: Yright
  if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 1, grid);
  // boundary condition: Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
    perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 2, grid);
  // boundary condition: Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
    perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 2, grid);

  // OpenBC
  BoundaryConditionsEImage(imageX, imageY, imageZ, vectX, vectY, vectZ, nxn, nyn, nzn, vct, grid);

  // move from physical space to krylov space
  phys2solver(im, imageX, imageY, imageZ, nxn, nyn, nzn);


}

/*! Calculate PI dot (vectX, vectY, vectZ) */
void EMfields3D::PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid * grid) {
  double beta, edotb, omcx, omcy, omcz, denom;
  beta = .5 * qom[ns] * dt / c;
  for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++) {
        omcx = beta * (Bxn[i][j][k] + Bx_ext[i][j][k]);
        omcy = beta * (Byn[i][j][k] + By_ext[i][j][k]);
        omcz = beta * (Bzn[i][j][k] + Bz_ext[i][j][k]);
        edotb = vectX[i][j][k] * omcx + vectY[i][j][k] * omcy + vectZ[i][j][k] * omcz;
        denom = 1 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        PIdotX[i][j][k] += (vectX[i][j][k] + (vectY[i][j][k] * omcz - vectZ[i][j][k] * omcy + edotb * omcx)) * denom;
        PIdotY[i][j][k] += (vectY[i][j][k] + (vectZ[i][j][k] * omcx - vectX[i][j][k] * omcz + edotb * omcy)) * denom;
        PIdotZ[i][j][k] += (vectZ[i][j][k] + (vectX[i][j][k] * omcy - vectY[i][j][k] * omcx + edotb * omcz)) * denom;
      }


}
/*! Calculate MU dot (vectX, vectY, vectZ) */
void EMfields3D::MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid * grid) {
  double beta, edotb, omcx, omcy, omcz, denom;
  for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++) {
        MUdotX[i][j][k] = 0.0;
        MUdotY[i][j][k] = 0.0;
        MUdotZ[i][j][k] = 0.0;
      }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int k = 1; k < nzn - 1; k++) {
          omcx = beta * (Bxn[i][j][k] + Bx_ext[i][j][k]);
          omcy = beta * (Byn[i][j][k] + By_ext[i][j][k]);
          omcz = beta * (Bzn[i][j][k] + Bz_ext[i][j][k]);
          edotb = vectX[i][j][k] * omcx + vectY[i][j][k] * omcy + vectZ[i][j][k] * omcz;
          denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
          MUdotX[i][j][k] += (vectX[i][j][k] + (vectY[i][j][k] * omcz - vectZ[i][j][k] * omcy + edotb * omcx)) * denom;
          MUdotY[i][j][k] += (vectY[i][j][k] + (vectZ[i][j][k] * omcx - vectX[i][j][k] * omcz + edotb * omcy)) * denom;
          MUdotZ[i][j][k] += (vectZ[i][j][k] + (vectX[i][j][k] * omcy - vectY[i][j][k] * omcx + edotb * omcz)) * denom;
        }

  }

}
/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smooth(double value, double ***vector, int type, Grid * grid, VirtualTopology3D * vct) {

  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {

    if (value != 1.0) {
      double alpha;
      int nx, ny, nz;
      switch (type) {
        case (0):
          nx = grid->getNXC();
          ny = grid->getNYC();
          nz = grid->getNZC();
          communicateCenterBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);

          break;
        case (1):
          nx = grid->getNXN();
          ny = grid->getNYN();
          nz = grid->getNZN();
          communicateNodeBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct);
          break;
      }
      double ***temp = newArr3(double, nx, ny, nz);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++)
          for (int k = 1; k < nz - 1; k++)
            temp[i][j][k] = value * vector[i][j][k] + alpha * (vector[i - 1][j][k] + vector[i + 1][j][k] + vector[i][j - 1][k] + vector[i][j + 1][k] + vector[i][j][k - 1] + vector[i][j][k + 1]);
      for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++)
          for (int k = 1; k < nz - 1; k++)
            vector[i][j][k] = temp[i][j][k];
      delArr3(temp, nx, ny);
    }
  }
}
/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smoothE(double value, VirtualTopology3D * vct, Collective *col) {

  int nvolte = 6;
  for (int icount = 1; icount < nvolte + 1; icount++) {
    if (value != 1.0) {
      double alpha;
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Exth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Ex[i][j][k] + alpha * (Ex[i - 1][j][k] + Ex[i + 1][j][k] + Ex[i][j - 1][k] + Ex[i][j + 1][k] + Ex[i][j][k - 1] + Ex[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Ex[i][j][k] = temp[i][j][k];
      // Eyth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Ey[i][j][k] + alpha * (Ey[i - 1][j][k] + Ey[i + 1][j][k] + Ey[i][j - 1][k] + Ey[i][j + 1][k] + Ey[i][j][k - 1] + Ey[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Ey[i][j][k] = temp[i][j][k];
      // Ezth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Ez[i][j][k] + alpha * (Ez[i - 1][j][k] + Ez[i + 1][j][k] + Ez[i][j - 1][k] + Ez[i][j + 1][k] + Ez[i][j][k - 1] + Ez[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Ez[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
}

/* SPECIES: Interpolation smoothing TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector type = 1 --> node based vector */
void EMfields3D::smooth(double value, double ****vector, int is, int type, Grid * grid, VirtualTopology3D * vct) {
  cout << "Smoothing for Species not implemented in 3D" << endl;
}

/*! fix the B boundary when running gem */
void EMfields3D::fixBgem(Grid * grid, VirtualTopology3D * vct) {
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Bxc[i][nyc - 2][k] = Bxc[i][nyc - 1][k];
        Bxc[i][nyc - 3][k] = Bxc[i][nyc - 1][k];
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z;
        Bzc[i][nyc - 2][k] = B0z;
        Bzc[i][nyc - 3][k] = B0z;
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bxc[i][1][k] = Bxc[i][0][k];
        Bxc[i][2][k] = Bxc[i][0][k];
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z;
        Bzc[i][1][k] = B0z;
        Bzc[i][2][k] = B0z;
      }
  }

}

/*! fix the B boundary when running forcefree */
void EMfields3D::fixBforcefree(Grid * grid, VirtualTopology3D * vct) {
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z / cosh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 2][k] = B0z / cosh((grid->getYC(i, nyc - 2, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 3][k] = B0z / cosh((grid->getYC(i, nyc - 3, k) - Ly / 2) / delta);
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z / cosh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bzc[i][1][k] = B0z / cosh((grid->getYC(i, 1, k) - Ly / 2) / delta);
        Bzc[i][2][k] = B0z / cosh((grid->getYC(i, 2, k) - Ly / 2) / delta);
      }
  }

}


/*! adjust densities on boundaries that are not periodic */
void EMfields3D::adjustNonPeriodicDensities(int is, VirtualTopology3D * vct) {
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
      for (int k = 1; k < nzn - 1; k++) {
        rhons[is][1][i][k] += rhons[is][1][i][k];
        Jxs  [is][1][i][k] += Jxs  [is][1][i][k];
        Jys  [is][1][i][k] += Jys  [is][1][i][k];
        Jzs  [is][1][i][k] += Jzs  [is][1][i][k];
        pXXsn[is][1][i][k] += pXXsn[is][1][i][k];
        pXYsn[is][1][i][k] += pXYsn[is][1][i][k];
        pXZsn[is][1][i][k] += pXZsn[is][1][i][k];
        pYYsn[is][1][i][k] += pYYsn[is][1][i][k];
        pYZsn[is][1][i][k] += pYZsn[is][1][i][k];
        pZZsn[is][1][i][k] += pZZsn[is][1][i][k];
      }
  }
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int k = 1; k < nzn - 1; k++) {
        rhons[is][i][1][k] += rhons[is][i][1][k];
        Jxs  [is][i][1][k] += Jxs  [is][i][1][k];
        Jys  [is][i][1][k] += Jys  [is][i][1][k];
        Jzs  [is][i][1][k] += Jzs  [is][i][1][k];
        pXXsn[is][i][1][k] += pXXsn[is][i][1][k];
        pXYsn[is][i][1][k] += pXYsn[is][i][1][k];
        pXZsn[is][i][1][k] += pXZsn[is][i][1][k];
        pYYsn[is][i][1][k] += pYYsn[is][i][1][k];
        pYZsn[is][i][1][k] += pYZsn[is][i][1][k];
        pZZsn[is][i][1][k] += pZZsn[is][i][1][k];
      }
  }
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++) {
        rhons[is][i][j][1] += rhons[is][i][j][1];
        Jxs  [is][i][j][1] += Jxs  [is][i][j][1];
        Jys  [is][i][j][1] += Jys  [is][i][j][1];
        Jzs  [is][i][j][1] += Jzs  [is][i][j][1];
        pXXsn[is][i][j][1] += pXXsn[is][i][j][1];
        pXYsn[is][i][j][1] += pXYsn[is][i][j][1];
        pXZsn[is][i][j][1] += pXZsn[is][i][j][1];
        pYYsn[is][i][j][1] += pYYsn[is][i][j][1];
        pYZsn[is][i][j][1] += pYZsn[is][i][j][1];
        pZZsn[is][i][j][1] += pZZsn[is][i][j][1];
      }
  }
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
      for (int k = 1; k < nzn - 1; k++) {
        rhons[is][nxn - 2][i][k] += rhons[is][nxn - 2][i][k];
        Jxs  [is][nxn - 2][i][k] += Jxs  [is][nxn - 2][i][k];
        Jys  [is][nxn - 2][i][k] += Jys  [is][nxn - 2][i][k];
        Jzs  [is][nxn - 2][i][k] += Jzs  [is][nxn - 2][i][k];
        pXXsn[is][nxn - 2][i][k] += pXXsn[is][nxn - 2][i][k];
        pXYsn[is][nxn - 2][i][k] += pXYsn[is][nxn - 2][i][k];
        pXZsn[is][nxn - 2][i][k] += pXZsn[is][nxn - 2][i][k];
        pYYsn[is][nxn - 2][i][k] += pYYsn[is][nxn - 2][i][k];
        pYZsn[is][nxn - 2][i][k] += pYZsn[is][nxn - 2][i][k];
        pZZsn[is][nxn - 2][i][k] += pZZsn[is][nxn - 2][i][k];
      }
  }
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int k = 1; k < nzn - 1; k++) {
        rhons[is][i][nyn - 2][k] += rhons[is][i][nyn - 2][k];
        Jxs  [is][i][nyn - 2][k] += Jxs  [is][i][nyn - 2][k];
        Jys  [is][i][nyn - 2][k] += Jys  [is][i][nyn - 2][k];
        Jzs  [is][i][nyn - 2][k] += Jzs  [is][i][nyn - 2][k];
        pXXsn[is][i][nyn - 2][k] += pXXsn[is][i][nyn - 2][k];
        pXYsn[is][i][nyn - 2][k] += pXYsn[is][i][nyn - 2][k];
        pXZsn[is][i][nyn - 2][k] += pXZsn[is][i][nyn - 2][k];
        pYYsn[is][i][nyn - 2][k] += pYYsn[is][i][nyn - 2][k];
        pYZsn[is][i][nyn - 2][k] += pYZsn[is][i][nyn - 2][k];
        pZZsn[is][i][nyn - 2][k] += pZZsn[is][i][nyn - 2][k];
      }
  }
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++) {
        rhons[is][i][j][nzn - 2] += rhons[is][i][j][nzn - 2];
        Jxs  [is][i][j][nzn - 2] += Jxs  [is][i][j][nzn - 2];
        Jys  [is][i][j][nzn - 2] += Jys  [is][i][j][nzn - 2];
        Jzs  [is][i][j][nzn - 2] += Jzs  [is][i][j][nzn - 2];
        pXXsn[is][i][j][nzn - 2] += pXXsn[is][i][j][nzn - 2];
        pXYsn[is][i][j][nzn - 2] += pXYsn[is][i][j][nzn - 2];
        pXZsn[is][i][j][nzn - 2] += pXZsn[is][i][j][nzn - 2];
        pYYsn[is][i][j][nzn - 2] += pYYsn[is][i][j][nzn - 2];
        pYZsn[is][i][j][nzn - 2] += pYZsn[is][i][j][nzn - 2];
        pZZsn[is][i][j][nzn - 2] += pZZsn[is][i][j][nzn - 2];
      }
  }
}

void EMfields3D::ConstantChargeOpenBCv2(Grid * grid, VirtualTopology3D * vct) {

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

    ff = 1.0;
    if (is == 0) ff = -1.0;

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = rhons[is][4][j][k];
          rhons[is][1][j][k] = rhons[is][4][j][k];
          rhons[is][2][j][k] = rhons[is][4][j][k];
          rhons[is][3][j][k] = rhons[is][4][j][k];
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-3][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-2][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-1][j][k] = rhons[is][nx-5][j][k];
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = rhons[is][i][4][k];
          rhons[is][i][1][k] = rhons[is][i][4][k];
          rhons[is][i][2][k] = rhons[is][i][4][k];
          rhons[is][i][3][k] = rhons[is][i][4][k];
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-3][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-2][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-1][k] = rhons[is][i][ny-5][k];
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = rhons[is][i][j][4];
          rhons[is][i][j][1] = rhons[is][i][j][4];
          rhons[is][i][j][2] = rhons[is][i][j][4];
          rhons[is][i][j][3] = rhons[is][i][j][4];
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-3] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-2] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-1] = rhons[is][i][j][nz-5];
        }
    }
  }

}

void EMfields3D::ConstantChargeOpenBC(Grid * grid, VirtualTopology3D * vct) {

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

    ff = 1.0;
    if (is == 0) ff = -1.0;

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][1][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][3][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-3][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-1][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][1][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][3][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-3][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-1][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][1] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][3] = ff * rhoINIT[is] / FourPI;
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-3] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-1] = ff * rhoINIT[is] / FourPI;
        }
    }
  }

}

void EMfields3D::ConstantChargePlanet(Grid * grid, VirtualTopology3D * vct, double R, double x_center, double y_center, double z_center) {

  double xd;
  double yd;
  double zd;
  double ff;

  for (int is = 0; is < ns; is++) {
    ff = 1.0;
    if (is == 0) ff = -1.0;
    for (int i = 1; i < nxn; i++) {
      for (int j = 1; j < nyn; j++) {
        for (int k = 1; k < nzn; k++) {

          xd = grid->getXN(i,j,k) - x_center;
          yd = grid->getYN(i,j,k) - y_center;
          zd = grid->getZN(i,j,k) - z_center;

          if ((xd*xd+yd*yd+zd*zd) <= R*R) {
            rhons[is][i][j][k] = ff * rhoINIT[is] / FourPI;
          }

        }
      }
    }
  }

}

/*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
void EMfields3D::calculateB(Grid * grid, VirtualTopology3D * vct, Collective *col) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** B CALCULATION ***" << endl;

  // calculate the curl of Eth
  grid->curlN2C(tempXC, tempYC, tempZC, Exth, Eyth, Ezth);
  // update the magnetic field
  addscale(-c * dt, 1, Bxc, tempXC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Byc, tempYC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Bzc, tempZC, nxc, nyc, nzc);
  // communicate ghost 
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  if (Case=="ForceFree") fixBforcefree(grid,vct);
  if (Case=="GEM")       fixBgem(grid, vct);
  if (Case=="GEMnoPert") fixBgem(grid, vct);

  // OpenBC:
  BoundaryConditionsB(Bxc,Byc,Bzc,nxc,nyc,nzc,grid,vct);

  // interpolate C2N
  grid->interpC2N(Bxn, Bxc);
  grid->interpC2N(Byn, Byc);
  grid->interpC2N(Bzn, Bzc);

  communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);


}
/*! initialize EM field with transverse electric waves 1D and rotate anticlockwise (theta degrees) */
void EMfields3D::initEM_rotate(VirtualTopology3D * vct, Grid * grid, Collective *col, double B, double theta) {
  // initialize E and rhos on nodes
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      Ex[i][j][0] = 0.0;
      Ey[i][j][0] = 0.0;
      Ez[i][j][0] = 0.0;
      Bxn[i][j][0] = B * cos(theta * M_PI / 180);
      Byn[i][j][0] = B * sin(theta * M_PI / 180);
      Bzn[i][j][0] = 0.0;
      rhons[0][i][j][0] = 0.07957747154595; // electrons: species is now first index
      rhons[1][i][j][0] = 0.07957747154595; // protons: species is now first index
    }
  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);


  for (int is = 0; is < ns; is++)
    grid->interpN2C(rhocs, is, rhons);

}
/*!Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid) {

  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  ne_mod *= alpha;
  ni_mod *= alpha;
  // cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      rhons[0][i][j][0] += ne_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ne_phase);
      rhons[1][i][j][0] += ni_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ni_phase);
    }

  for (int is = 0; is < ns; is++)
    grid->interpN2C(rhocs, is, rhons);
}


/*!Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid) {

  double alpha;

  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;

  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      Ex[i][j][0] += Ex_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ex_phase);
      Ey[i][j][0] += Ey_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ey_phase);
      Ez[i][j][0] += Ez_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ez_phase);
      Bxn[i][j][0] += Bx_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bx_phase);
      Byn[i][j][0] += By_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + By_phase);
      Bzn[i][j][0] += Bz_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bz_phase);

    }

  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);


}


/*! Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
void EMfields3D::calculateHatFunctions(Grid * grid, VirtualTopology3D * vct) {
  // smoothing
  smooth(Smooth, rhoc, 0, grid, vct);
  // calculate j hat

  for (int is = 0; is < ns; is++) {
    grid->divSymmTensorN2C(tempXC, tempYC, tempZC, pXXsn, pXYsn, pXZsn, pYYsn, pYZsn, pZZsn, is);

    scale(tempXC, -dt / 2.0, nxc, nyc, nzc);
    scale(tempYC, -dt / 2.0, nxc, nyc, nzc);
    scale(tempZC, -dt / 2.0, nxc, nyc, nzc);
    // communicate before interpolating
    communicateCenterBC_P(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
    communicateCenterBC_P(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
    communicateCenterBC_P(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);

    grid->interpC2N(tempXN, tempXC);
    grid->interpC2N(tempYN, tempYC);
    grid->interpC2N(tempZN, tempZC);
    sum(tempXN, Jxs, nxn, nyn, nzn, is);
    sum(tempYN, Jys, nxn, nyn, nzn, is);
    sum(tempZN, Jzs, nxn, nyn, nzn, is);
    // PIDOT
    PIdot(Jxh, Jyh, Jzh, tempXN, tempYN, tempZN, is, grid);

  }
  // smooth j
  smooth(Smooth, Jxh, 1, grid, vct);
  smooth(Smooth, Jyh, 1, grid, vct);
  smooth(Smooth, Jzh, 1, grid, vct);

  // calculate rho hat = rho - (dt*theta)div(jhat)
  grid->divN2C(tempXC, Jxh, Jyh, Jzh);
  scale(tempXC, -dt * th, nxc, nyc, nzc);
  sum(tempXC, rhoc, nxc, nyc, nzc);
  eq(rhoh, tempXC, nxc, nyc, nzc);
  // communicate rhoh
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct);
}
/*! Image of Poisson Solver */
void EMfields3D::PoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct) {
  // allocate 2 three dimensional service vectors
  double ***temp = newArr3(double, nxc, nyc, nzc);
  double ***im = newArr3(double, nxc, nyc, nzc);
  eqValue(0.0, image, (nxc - 2) * (nyc - 2) * (nzc - 2));
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, nxc, nyc, nzc);
  // calculate the laplacian
  grid->lapC2Cpoisson(im, temp, vct);
  // move from physical space to krylov space
  phys2solver(image, im, nxc, nyc, nzc);
  // deallocate temporary array and objects
  delArr3(temp, nxc, nyc);
  delArr3(im, nxc, nyc);
}
/*! interpolate charge density and pressure density from node to center */
void EMfields3D::interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid) {
  // do we need communication or not really?
  grid->interpN2C(rhoc, rhon);
}
/*! communicate ghost for grid -> Particles interpolation */
void EMfields3D::communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct) {
  // interpolate adding common nodes among processors
  timeTasks.start_communicate();

  communicateInterp(nxn, nyn, nzn, ns, rhons, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jxs, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jys, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, Jzs, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXXsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXYsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pXZsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pYYsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pYZsn, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, pZZsn, 0, 0, 0, 0, 0, 0, vct);
  // calculate the correct densities on the boundaries
  adjustNonPeriodicDensities(ns, vct);
  // put the correct values on ghost cells
  timeTasks.addto_communicate();

  communicateNode_P(nxn, nyn, nzn, rhons, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jxs, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jys, ns, vct);
  communicateNode_P(nxn, nyn, nzn, Jzs, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXXsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXYsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pXZsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pYYsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pYZsn, ns, vct);
  communicateNode_P(nxn, nyn, nzn, pZZsn, ns, vct);
}

/* add moments (e.g. from an OpenMP thread) to the accumulated moments */
void EMfields3D::addToSpeciesMoments(const Moments & in, int is) {
  assert_eq(in.get_nx(), nxn);
  assert_eq(in.get_ny(), nyn);
  assert_eq(in.get_nz(), nzn);
  for (register int i = 0; i < nxn; i++) {
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
        rhons[is][i][j][k] += invVOL*in.get_rho(i, j, k);
        Jxs  [is][i][j][k] += invVOL*in.get_Jx(i, j, k);
        Jys  [is][i][j][k] += invVOL*in.get_Jy(i, j, k);
        Jzs  [is][i][j][k] += invVOL*in.get_Jz(i, j, k);
        pXXsn[is][i][j][k] += invVOL*in.get_pXX(i, j, k);
        pXYsn[is][i][j][k] += invVOL*in.get_pXY(i, j, k);
        pXZsn[is][i][j][k] += invVOL*in.get_pXZ(i, j, k);
        pYYsn[is][i][j][k] += invVOL*in.get_pYY(i, j, k);
        pYZsn[is][i][j][k] += invVOL*in.get_pYZ(i, j, k);
        pZZsn[is][i][j][k] += invVOL*in.get_pZZ(i, j, k);
      }
  }
}

/*! set to 0 all the densities fields */
void EMfields3D::setZeroDensities() {
  for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
        Jx  [i][j][k] = 0.0;
        Jxh [i][j][k] = 0.0;
        Jy  [i][j][k] = 0.0;
        Jyh [i][j][k] = 0.0;
        Jz  [i][j][k] = 0.0;
        Jzh [i][j][k] = 0.0;
        rhon[i][j][k] = 0.0;
      }
  for (register int i = 0; i < nxc; i++)
    for (register int j = 0; j < nyc; j++)
      for (register int k = 0; k < nzc; k++) {
        rhoc[i][j][k] = 0.0;
        rhoh[i][j][k] = 0.0;
      }
  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          rhons[kk][i][j][k] = 0.0;
          Jxs  [kk][i][j][k] = 0.0;
          Jys  [kk][i][j][k] = 0.0;
          Jzs  [kk][i][j][k] = 0.0;
          pXXsn[kk][i][j][k] = 0.0;
          pXYsn[kk][i][j][k] = 0.0;
          pXZsn[kk][i][j][k] = 0.0;
          pYYsn[kk][i][j][k] = 0.0;
          pYZsn[kk][i][j][k] = 0.0;
          pZZsn[kk][i][j][k] = 0.0;
        }

}
/*!SPECIES: Sum the charge density of different species on NODES */
void EMfields3D::sumOverSpecies(VirtualTopology3D * vct) {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++)
          rhon[i][j][k] += rhons[is][i][j][k];
}

/*!SPECIES: Sum current density for different species */
void EMfields3D::sumOverSpeciesJ() {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          Jx[i][j][k] += Jxs[is][i][j][k];
          Jy[i][j][k] += Jys[is][i][j][k];
          Jz[i][j][k] += Jzs[is][i][j][k];
        }
}



/*! initialize Magnetic and Electric Field with initial configuration */
void EMfields3D::init(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++) {
      for (int j = 0; j < nyn; j++) {
        for (int k = 0; k < nzn; k++) {
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = B0x;
          Byn[i][j][k] = B0y;
          Bzn[i][j][k] = B0z;
        }
      }
    }

    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {                        // READING FROM RESTART
    if (vct->getCartesian_rank() == 0)
      cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << name_file << endl;
      cout << "RESTART NOT POSSIBLE" << endl;
    }

    dataset_id = H5Dopen2(file_id, "/fields/Bx/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);



    // Bxn
    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    dataset_id = H5Dopen2(file_id, "/fields/By/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    dataset_id = H5Dopen2(file_id, "/fields/Bz/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    dataset_id = H5Dopen2(file_id, "/fields/Ex/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey 
    dataset_id = H5Dopen2(file_id, "/fields/Ey/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez 
    dataset_id = H5Dopen2(file_id, "/fields/Ez/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species

    stringstream *species_name = new stringstream[ns];
    for (int is = 0; is < ns; is++) {
      species_name[is] << is;
      string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
      dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int jj = 1; jj < nzn - 1; jj++)
            rhons[is][i][j][jj] = temp_storage[k++];
      communicateNode_P(nxn, nyn, nzn, rhons, is, vct);
      status = H5Dclose(dataset_id);

    }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // communicate E
    communicateNodeBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
    delete[]species_name;
  }
}

#ifdef BATSRUS
/*! initiliaze EM for GEM challange */
void EMfields3D::initBATSRUS(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  cout << "------------------------------------------" << endl;
  cout << "         Initialize from BATSRUS          " << endl;
  cout << "------------------------------------------" << endl;

  // loop over species and cell centers: fill in charge density
  for (int is=0; is < ns; is++)
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
        for (int k=0; k < nzc; k++)
        {
          // WARNING getFluidRhoCenter contains "case" statment
          rhocs[is][i][j][k] = col->getFluidRhoCenter(i,j,k,is);
        }

  // loop over cell centers and fill in magnetic and electric fields
  for (int i=0; i < nxc; i++)
    for (int j=0; j < nyc; j++)
      for (int k=0; k < nzc; k++)
      {
        // WARNING getFluidRhoCenter contains "case" statment
        col->setFluidFieldsCenter(&Ex[i][j][k],&Ey[i][j][k],&Ez[i][j][k],
            &Bxc[i][j][k],&Byc[i][j][k],&Bzc[i][j][k],i,j,k);
      }

  // interpolate from cell centers to nodes (corners of cells)
  for (int is=0 ; is<ns; is++)
    grid->interpC2N(rhons[is],rhocs[is]);
  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);
}
#endif

/*! initiliaze EM for GEM challange */
void EMfields3D::initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
          Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXN(i, j, k) - Lx / 2;
          ypert = grid->getYN(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly );
          Byc[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXC(i, j, k) - Lx / 2;
          ypert = grid->getYC(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxc[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byc[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initOriginalGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          const double yM = grid->getYN(i, j, k) - .5 * Ly;
          Bxn[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          Bxn[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byn[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          const double yM = grid->getYC(i, j, k) - .5 * Ly;
          Bxc[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXC(i, j, k) - .5 * Lx;
          Bxc[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byc[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzc[i][j][k] = B0z;
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double deltax = 8. * delta;
  const double deltay = 4. * delta;
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              const double sech_yBd = 1. / cosh(yBd);
              const double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxn[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byn[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxn[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byn[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);

          // guide field
          Bzn[i][j][k] = B0z;
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          Bxc[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxc[i][j][k] += 0.;
          Byc[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxc[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byc[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxc[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byc[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);
          // guide field
          Bzc[i][j][k] = B0z;
        }
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


/*! initialize GEM challenge with no Perturbation with dipole-like tail topology */
void EMfields3D::initGEMDipoleLikeTailNoPert(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  // parameters controling the field topology
  // e.g., x1=Lx/5,x2=Lx/4 give 'separated' fields, x1=Lx/4,x2=Lx/3 give 'reconnected' topology

  double x1 = Lx / 6.0;         // minimal position of the gaussian peak 
  double x2 = Lx / 4.0;         // maximal position of the gaussian peak (the one closer to the center)
  double sigma = Lx / 15;       // base sigma of the gaussian - later it changes with the grid
  double stretch_curve = 2.0;   // stretch the sin^2 function over the x dimension - also can regulate the number of 'knots/reconnecitons points' if less than 1
  double skew_parameter = 0.50; // skew of the shape of the gaussian
  double pi = 3.1415927;
  double r1, r2, delta_x1x2;

  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYN(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYN(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYC(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYC(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }

}

/*! initialize GEM challenge with no Perturbation */
void EMfields3D::initGEMnoPert(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initRandomField(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  double **modes_seed = newArr2(double, 7, 7);
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    double phixy;
    double phix;
    double phiy;
    double phiz;
    double kx;
    double ky;
    phixy = rand() / (double) RAND_MAX;
    phiz = rand() / (double) RAND_MAX;
    phix = rand() / (double) RAND_MAX;
    phiy = rand() / (double) RAND_MAX;
    for (int m = -3; m < 4; m++)
      for (int n = -3; n < 4; n++) {
        modes_seed[m + 3][n + 3] = rand() / (double) RAND_MAX;
      }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = 0.0;
          Byn[i][j][k] = 0.0;
          Bzn[i][j][k] = 0.0;
          for (int m = -3; m < 4; m++)
            for (int n = -3; n < 4; n++) {

              kx = 2.0 * M_PI * m / Lx;
              ky = 2.0 * M_PI * n / Ly;
              Bxn[i][j][k] += -B0x * ky * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
              Byn[i][j][k] += B0x * kx * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
              Bzn[i][j][k] += B0x * cos(grid->getXN(i, j, k) * kx + grid->getYN(i, j, k) * ky + 2.0 * M_PI * modes_seed[m + 3][n + 3]);
            }

          /* for (int m=1; m < 4; m++) for (int n=1; n < 4; n++){ kx=2.0*M_PI*m/Lx; ky=2.0*M_PI*n/Ly; Bxn[i][j][k] += B0x/kx*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phixy); Byn[i][j][k] += B0x/ky*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phixy); Bzn[i][j][k] += B0x/(kx+ky)*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*phiz); } for(int n=1; n < 4; n++){ ky=2.0*M_PI*n/Ly; Bxn[i][j][k] += B0x/(2.0*M_PI/Lx)*cos(grid->getYN(i,j,k)*ky+2.0*M_PI*phix); } for(int m=1; m < 4; m++){ kx=2.0*M_PI*m/Lx; Byn[i][j][k] += B0x/(2.0*M_PI/Ly)*cos(grid->getXN(i,j,k)*kx+2.0*M_PI*phiy); } */
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
  delArr2(modes_seed, 7);
}

/*! Init Force Free (JxB=0) */
void EMfields3D::initForceFree(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  if (restart1 == 0) {

    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "----------------------------------------" << endl;
      cout << "Initialize Force Free with Perturbation" << endl;
      cout << "----------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
      }
      cout << "Smoothing Factor = " << Smooth << endl;
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          Bxn[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXN(i, j, k) / Lx) * sin(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          Byn[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXN(i, j, k) / Lx) * cos(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzn[i][j][k] = B0z / cosh((grid->getYN(i, j, k) - Ly / 2) / delta);
        }
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the perturbation
          Bxc[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXC(i, j, k) / Lx) * sin(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          Byc[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXC(i, j, k) / Lx) * cos(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzc[i][j][k] = B0z / cosh((grid->getYC(i, j, k) - Ly / 2) / delta);
        }

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}
/*! Initialize the EM field with constants values or from restart */
void EMfields3D::initBEAM(VirtualTopology3D * vct, Grid * grid, Collective *col, double x_center, double y_center, double z_center, double radius) {
  double distance;
  // initialize E and rhos on nodes
  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = 0.0;
          Byn[i][j][k] = 0.0;
          Bzn[i][j][k] = 0.0;
          distance = (grid->getXN(i, j, k) - x_center) * (grid->getXN(i, j, k) - x_center) / (radius * radius) + (grid->getYN(i, j, k) - y_center) * (grid->getYN(i, j, k) - y_center) / (radius * radius) + (grid->getZN(i, j, k) - z_center) * (grid->getZN(i, j, k) - z_center) / (4 * radius * radius);
          // plasma
          rhons[0][i][j][k] = rhoINIT[0] / FourPI;  // initialize with constant density
          // electrons
          rhons[1][i][j][k] = rhoINIT[1] / FourPI;
          // beam
          if (distance < 1.0)
            rhons[2][i][j][k] = rhoINIT[2] / FourPI;
          else
            rhons[2][i][j][k] = 0.0;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = 0.0;
          Byc[i][j][k] = 0.0;
          Bzc[i][j][k] = 0.0;


        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {                        // EM initialization from RESTART
    init(vct, grid, col);            // use the fields from restart file
  }

}

/*! Initialise a combination of magnetic dipoles */
void EMfields3D::initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col){

  double distance;

  if (restart1 == 0){

    double ebc[3];
    cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
    scale(ebc,-1.0,3);

    for (int i=0; i < nxn; i++){
      for (int j=0; j < nyn; j++){
        for (int k=0; k < nzn; k++){
          for (int is=0; is < ns; is++){
            rhons[is][i][j][k] = rhoINIT[is]/FourPI;
          }
          Ex[i][j][k] = ebc[0];
          Ey[i][j][k] = ebc[1];
          Ez[i][j][k] = ebc[2];

          double blp[3];
          // Set coil diameter
          double a=delta;

          double xc=x_center;
          double yc=y_center;
          double zc=z_center;

          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
          double z = grid->getZN(i,j,k);

          // Compute dipolar field B_ext
          loopZ(blp, x, y, z, a, xc, yc, zc, B1z);
          Bx_ext[i][j][k]  = blp[0];
          By_ext[i][j][k]  = blp[1];
          Bz_ext[i][j][k]  = blp[2];
          loopX(blp, x, y, z, a, xc, yc, zc, B1x);
          Bx_ext[i][j][k] += blp[0];
          By_ext[i][j][k] += blp[1];
          Bz_ext[i][j][k] += blp[2];
          loopY(blp, x, y, z, a, xc, yc, zc, B1y);
          Bx_ext[i][j][k] += blp[0];
          By_ext[i][j][k] += blp[1];
          Bz_ext[i][j][k] += blp[2];

          Bxn[i][j][k] = B0x + Bx_ext[i][j][k];
          Byn[i][j][k] = B0y + By_ext[i][j][k];
          Bzn[i][j][k] = B0z + Bz_ext[i][j][k];

          // -- Uncomment if using the J_ext method:
          Bx_ext[i][j][k]  = 0.0;
          By_ext[i][j][k]  = 0.0;
          Bz_ext[i][j][k]  = 0.0;
          // -- end Uncomment
        }
      }
    }

    grid->interpN2C(Bxc,Bxn);
    grid->interpN2C(Byc,Byn);
    grid->interpN2C(Bzc,Bzn);

    communicateCenterBC_P(nxc,nyc,nzc,Bxc,col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5],vct);
    communicateCenterBC_P(nxc,nyc,nzc,Byc,col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5],vct);
    communicateCenterBC_P(nxc,nyc,nzc,Bzc,col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5],vct);

    // initialize J_ext =c/4*pi curl(B) on nodes (current due to the dipole)
    grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc);
    scale(Jx_ext,tempXN,c/FourPI,nxn,nyn,nzn);
    scale(Jy_ext,tempYN,c/FourPI,nxn,nyn,nzn);
    scale(Jz_ext,tempZN,c/FourPI,nxn,nyn,nzn);

    for (int is=0 ; is<ns; is++)
      grid->interpN2C(rhocs,is,rhons);

  }
  else { // EM initialization from RESTART
    init(vct,grid,col);  // use the fields from restart file
  }

}


/*! Calculate the susceptibility on the boundary leftX */
void EMfields3D::sustensorLeftX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[1][j][k];
        omcy = beta * Byn[1][j][k];
        omcz = beta * Bzn[1][j][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][1][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary rightX */
void EMfields3D::sustensorRightX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[nxn - 2][j][k];
        omcy = beta * Byn[nxn - 2][j][k];
        omcz = beta * Bzn[nxn - 2][j][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][nxn - 2][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[i][1][k];
        omcy = beta * Byn[i][1][k];
        omcz = beta * Bzn[i][1][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][1][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * Bxn[i][nyn - 2][k];
        omcy = beta * Byn[i][nyn - 2][k];
        omcz = beta * Bzn[i][nyn - 2][k];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][nyn - 2][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * Bxn[i][j][1];
        omcy = beta * Byn[i][j][1];
        omcz = beta * Bzn[i][j][1];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][1] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * Bxn[i][j][nzn - 2];
        omcy = beta * Byn[i][j][nzn - 2];
        omcz = beta * Bzn[i][j][nzn - 2];
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][nyn - 2] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }
}

/*! Perfect conductor boundary conditions: LEFT wall */
void EMfields3D::perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0:  // boundary condition on X-DIRECTION 
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorLeftX(susxx, susyx, suszx);
      for (int i=1; i <  nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[1][i][j] = vectorX[1][i][j] - (Ex[1][i][j] - susyx[i][j]*vectorY[1][i][j] - suszx[i][j]*vectorZ[1][i][j] - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[1][i][j] = vectorY[1][i][j] - 0.0*vectorY[2][i][j];
          imageZ[1][i][j] = vectorZ[1][i][j] - 0.0*vectorZ[2][i][j];
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorLeftY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[i][1][j] = vectorX[i][1][j] - 0.0*vectorX[i][2][j];
          imageY[i][1][j] = vectorY[i][1][j] - (Ey[i][1][j] - susxy[i][j]*vectorX[i][1][j] - suszy[i][j]*vectorZ[i][1][j] - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][1][j] = vectorZ[i][1][j] - 0.0*vectorZ[i][2][j];
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorLeftZ(susxy, susyy, suszy);
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          imageX[i][j][1] = vectorX[i][j][1];
          imageY[i][j][1] = vectorX[i][j][1];
          imageZ[i][j][1] = vectorZ[i][j][1] - (Ez[i][j][1] - susxz[i][j]*vectorX[i][j][1] - susyz[i][j]*vectorY[i][j][1] - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions: RIGHT wall */
void EMfields3D::perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  double beta, omcx, omcy, omcz, denom;
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorRightX(susxx, susyx, suszx);
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[nxn-2][i][j] = vectorX[nxn-2][i][j] - (Ex[nxn-2][i][j] - susyx[i][j]*vectorY[nxn-2][i][j] - suszx[i][j]*vectorZ[nxn-2][i][j] - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[nxn-2][i][j] = vectorY[nxn-2][i][j] - 0.0 * vectorY[nxn-3][i][j];
          imageZ[nxn-2][i][j] = vectorZ[nxn-2][i][j] - 0.0 * vectorZ[nxn-3][i][j];
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);       
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorRightY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          imageX[i][nyn-2][j] = vectorX[i][nyn-2][j] - 0.0*vectorX[i][nyn-3][j];
          imageY[i][nyn-2][j] = vectorY[i][nyn-2][j] - (Ey[i][nyn-2][j] - susxy[i][j]*vectorX[i][nyn-2][j] - suszy[i][j]*vectorZ[i][nyn-2][j] - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][nyn-2][j] = vectorZ[i][nyn-2][j] - 0.0*vectorZ[i][nyn-3][j];
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION RIGHT
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorRightZ(susxz, susyz, suszz);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nyn-1;j++){
          imageX[i][j][nzn-2] = vectorX[i][j][nzn-2];
          imageY[i][j][nzn-2] = vectorY[i][j][nzn-2];
          imageZ[i][j][nzn-2] = vectorZ[i][j][nzn-2] - (Ez[i][j][nzn-2] - susxz[i][j]*vectorX[i][j][nzn-2] - susyz[i][j]*vectorY[i][j][nzn-2] - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);       
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION LEFT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[1][i][j] = 0.0;
          vectorY[1][i][j] = ebc[1];
          vectorZ[1][i][j] = ebc[2];
          //+//          vectorX[1][i][j] = 0.0;
          //+//          vectorY[1][i][j] = 0.0;
          //+//          vectorZ[1][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][1][j] = ebc[0];
          vectorY[i][1][j] = 0.0;
          vectorZ[i][1][j] = ebc[2];
          //+//          vectorX[i][1][j] = 0.0;
          //+//          vectorY[i][1][j] = 0.0;
          //+//          vectorZ[i][1][j] = 0.0;
        }
      break;
    case 2: // boundary condition on Z-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][1] = ebc[0];
          vectorY[i][j][1] = ebc[1];
          vectorZ[i][j][1] = 0.0;
          //+//          vectorX[i][j][1] = 0.0;
          //+//          vectorY[i][j][1] = 0.0;
          //+//          vectorZ[i][j][1] = 0.0;
        }
      break;
  }
}

/*! Perfect conductor boundary conditions for source: RIGHT WALL */
void EMfields3D::perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[nxn-2][i][j] = 0.0;
          vectorY[nxn-2][i][j] = ebc[1];
          vectorZ[nxn-2][i][j] = ebc[2];
          //+//          vectorX[nxn-2][i][j] = 0.0;
          //+//          vectorY[nxn-2][i][j] = 0.0;
          //+//          vectorZ[nxn-2][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][nyn-2][j] = ebc[0];
          vectorY[i][nyn-2][j] = 0.0;
          vectorZ[i][nyn-2][j] = ebc[2];
          //+//          vectorX[i][nyn-2][j] = 0.0;
          //+//          vectorY[i][nyn-2][j] = 0.0;
          //+//          vectorZ[i][nyn-2][j] = 0.0;
        }
      break;
    case 2:
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][nzn-2] = ebc[0];
          vectorY[i][j][nzn-2] = ebc[1];
          vectorZ[i][j][nzn-2] = 0.0;
          //+//          vectorX[i][j][nzn-2] = 0.0;
          //+//          vectorY[i][j][nzn-2] = 0.0;
          //+//          vectorZ[i][j][nzn-2] = 0.0;
        }
      break;
  }
}


// OpenBCs

injInfoFields* EMfields3D::get_InfoFieldsTop() {return injFieldsTop;}
injInfoFields* EMfields3D::get_InfoFieldsBottom() {return injFieldsBottom;}
injInfoFields* EMfields3D::get_InfoFieldsLeft() {return injFieldsLeft;}
injInfoFields* EMfields3D::get_InfoFieldsRight() {return injFieldsRight;}
injInfoFields* EMfields3D::get_InfoFieldsFront() {return injFieldsFront;}
injInfoFields* EMfields3D::get_InfoFieldsRear() {return injFieldsRear;}

// Open Boundary conditions implementation

void EMfields3D::updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col){

  double u_0, v_0, w_0;
  u_0=col->getU0(0);
  v_0=col->getV0(0);
  w_0=col->getW0(0);

  if (vct->getXleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< 3;i++)
      for (int j=0; j<nyn;j++)
        for (int k=0; k<nzn;k++){

          injFieldsLeft->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsLeft->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsLeft->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsLeft->BxITemp[i][j][k]=B0x;
          injFieldsLeft->ByITemp[i][j][k]=B0y;
          injFieldsLeft->BzITemp[i][j][k]=B0z;
        }
  }

  if (vct->getXright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=nxn-3; i< nxn; i++)
      for (int j=0; j<nyn; j++)
        for (int k=0; k<nzn; k++){

          injFieldsRight->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsRight->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsRight->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsRight->BxITemp[i][j][k]=B0x;
          injFieldsRight->ByITemp[i][j][k]=B0y;
          injFieldsRight->BzITemp[i][j][k]=B0z;

        }

  }

  if (vct->getYleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<3;j++)
        for (int k=0; k<nzn;k++){

          injFieldsBottom->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsBottom->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsBottom->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsBottom->BxITemp[i][j][k]=B0x;
          injFieldsBottom->ByITemp[i][j][k]=B0y;
          injFieldsBottom->BzITemp[i][j][k]=B0z;
        }

  }
  if (vct->getYright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=nyn-3; j<nyn;j++)
        for (int k=0; k<nzn;k++){

          injFieldsTop->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsTop->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsTop->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsTop->BxITemp[i][j][k]=B0x;
          injFieldsTop->ByITemp[i][j][k]=B0y;
          injFieldsTop->BzITemp[i][j][k]=B0z;
        }

  }
  if (vct->getZleft_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<nyn;j++)
        for (int k=0; k<3;k++){

          injFieldsRear->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsRear->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsRear->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsRear->BxITemp[i][j][k]=B0x;
          injFieldsRear->ByITemp[i][j][k]=B0y;
          injFieldsRear->BzITemp[i][j][k]=B0z;
        }

  }

  if (vct->getZright_neighbor() == MPI_PROC_NULL)
  {
    for (int i=0; i< nxn;i++)
      for (int j=0; j<nyn;j++)
        for (int k=nzn-3; k<nzn;k++){

          injFieldsFront->ExITemp[i][j][k]=w_0*B0y-v_0*B0z;
          injFieldsFront->EyITemp[i][j][k]=u_0*B0z-w_0*B0x;
          injFieldsFront->EzITemp[i][j][k]=v_0*B0x-u_0*B0y;

          injFieldsFront->BxITemp[i][j][k]=B0x;
          injFieldsFront->ByITemp[i][j][k]=B0y;
          injFieldsFront->BzITemp[i][j][k]=B0z;
        }
  }

}

void EMfields3D::BoundaryConditionsEImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz, VirtualTopology3D *vct,Grid *grid){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[0][j][k] = vectorX[0][j][k] - injFieldsLeft->ExITemp[0][j][k];
        imageY[0][j][k] = vectorY[0][j][k] - injFieldsLeft->EyITemp[0][j][k];
        imageZ[0][j][k] = vectorZ[0][j][k] - injFieldsLeft->EzITemp[0][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[nx-1][j][k] = vectorX[nx-1][j][k]- injFieldsRight->ExITemp[nx-1][j][k];
        imageY[nx-1][j][k] = vectorY[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];
        imageZ[nx-1][j][k] = vectorZ[nx-1][j][k]- injFieldsRight->EyITemp[nx-1][j][k];

      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][0][k] = vectorX[i][0][k]-injFieldsBottom->ExITemp[i][0][k];
        imageY[i][0][k] = vectorY[i][0][k]-injFieldsBottom->EyITemp[i][0][k];
        imageZ[i][0][k] = vectorZ[i][0][k]-injFieldsBottom->EzITemp[i][0][k];
      }

  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][ny-1][k] = vectorX[i][ny-1][k]-injFieldsTop->ExITemp[i][nx-1][k];
        imageY[i][ny-1][k] = vectorY[i][ny-1][k]-injFieldsTop->EyITemp[i][nx-1][k];
        imageZ[i][ny-1][k] = vectorZ[i][ny-1][k]-injFieldsTop->EzITemp[i][nx-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][0] = vectorX[i][j][0]-injFieldsFront->ExITemp[i][j][0];
        imageY[i][j][0] = vectorY[i][j][0]-injFieldsFront->EyITemp[i][j][0];
        imageZ[i][j][0] = vectorZ[i][j][0]-injFieldsFront->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][nz-1] = vectorX[i][j][nz-1]-injFieldsRear->ExITemp[i][j][nz-1];
        imageY[i][j][nz-1] = vectorY[i][j][nz-1]-injFieldsRear->EyITemp[i][j][nz-1];
        imageZ[i][j][nz-1] = vectorZ[i][j][nz-1]-injFieldsRear->EzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsB(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
        vectorX[0][j][k] = injFieldsLeft->BxITemp[0][j][k];
        vectorY[0][j][k] = injFieldsLeft->ByITemp[0][j][k];
        vectorZ[0][j][k] = injFieldsLeft->BzITemp[0][j][k];

//      vectorX[1][j][k] = injFieldsLeft->BxITemp[1][j][k];
//      vectorY[1][j][k] = injFieldsLeft->ByITemp[1][j][k];
//      vectorZ[1][j][k] = injFieldsLeft->BzITemp[1][j][k];
      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
//      vectorX[nx-2][j][k] = injFieldsRight->BxITemp[nx-2][j][k];
//      vectorY[nx-2][j][k] = injFieldsRight->ByITemp[nx-2][j][k];
//      vectorZ[nx-2][j][k] = injFieldsRight->BzITemp[nx-2][j][k];

        vectorX[nx-1][j][k] = injFieldsRight->BxITemp[nx-1][j][k];
        vectorY[nx-1][j][k] = injFieldsRight->ByITemp[nx-1][j][k];
        vectorZ[nx-1][j][k] = injFieldsRight->BzITemp[nx-1][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][1][k] = injFieldsBottom->BxITemp[i][1][k];
//      vectorY[i][1][k] = injFieldsBottom->ByITemp[i][1][k];
//      vectorZ[i][1][k] = injFieldsBottom->BzITemp[i][1][k];

        vectorX[i][0][k] = injFieldsBottom->BxITemp[i][0][k];
        vectorY[i][0][k] = injFieldsBottom->ByITemp[i][0][k];
        vectorZ[i][0][k] = injFieldsBottom->BzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][ny-2][k] = injFieldsTop->BxITemp[i][ny-2][k];
//      vectorY[i][ny-2][k] = injFieldsTop->ByITemp[i][ny-2][k];
//      vectorZ[i][ny-2][k] = injFieldsTop->BzITemp[i][ny-2][k];

        vectorX[i][ny-1][k] = injFieldsTop->BxITemp[i][ny-1][k];
        vectorY[i][ny-1][k] = injFieldsTop->ByITemp[i][ny-1][k];
        vectorZ[i][ny-1][k] = injFieldsTop->BzITemp[i][ny-1][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][1] = injFieldsRear->BxITemp[i][j][1];
//      vectorY[i][j][1] = injFieldsRear->ByITemp[i][j][1];
//      vectorZ[i][j][1] = injFieldsRear->BzITemp[i][j][1];

        vectorX[i][j][0] = injFieldsRear->BxITemp[i][j][0];
        vectorY[i][j][0] = injFieldsRear->ByITemp[i][j][0];
        vectorZ[i][j][0] = injFieldsRear->BzITemp[i][j][0];
      }
  }


  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][nz-2] = injFieldsFront->BxITemp[i][j][nz-2];
//      vectorY[i][j][nz-2] = injFieldsFront->ByITemp[i][j][nz-2];
//      vectorZ[i][j][nz-2] = injFieldsFront->BzITemp[i][j][nz-2];

        vectorX[i][j][nz-1] = injFieldsFront->BxITemp[i][j][nz-1];
        vectorY[i][j][nz-1] = injFieldsFront->ByITemp[i][j][nz-1];
        vectorZ[i][j][nz-1] = injFieldsFront->BzITemp[i][j][nz-1];
      }
  }

}

void EMfields3D::BoundaryConditionsE(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct){

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
        vectorX[1][j][k] = injFieldsLeft->ExITemp[1][j][k];
        vectorY[1][j][k] = injFieldsLeft->EyITemp[1][j][k];
        vectorZ[1][j][k] = injFieldsLeft->EzITemp[1][j][k];

//      vectorX[0][j][k] = injFieldsLeft->ExITemp[0][j][k];
//      vectorY[0][j][k] = injFieldsLeft->EyITemp[0][j][k];
//      vectorZ[0][j][k] = injFieldsLeft->EzITemp[0][j][k];
      } 
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){

//      vectorX[nx-2][j][k] = injFieldsRight->ExITemp[1][j][k];
//      vectorY[nx-2][j][k] = injFieldsRight->EyITemp[1][j][k];
//      vectorZ[nx-2][j][k] = injFieldsRight->EzITemp[1][j][k];

        vectorX[nx-1][j][k] = injFieldsRight->ExITemp[0][j][k];
        vectorY[nx-1][j][k] = injFieldsRight->EyITemp[0][j][k];
        vectorZ[nx-1][j][k] = injFieldsRight->EzITemp[0][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][1][k] = injFieldsBottom->ExITemp[i][1][k];
//      vectorY[i][1][k] = injFieldsBottom->EyITemp[i][1][k];
//      vectorZ[i][1][k] = injFieldsBottom->EzITemp[i][1][k];

        vectorX[i][0][k] = injFieldsBottom->ExITemp[i][0][k];
        vectorY[i][0][k] = injFieldsBottom->EyITemp[i][0][k];
        vectorZ[i][0][k] = injFieldsBottom->EzITemp[i][0][k];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
//      vectorX[i][ny-2][k] = injFieldsTop->ExITemp[i][1][k];
//      vectorY[i][ny-2][k] = injFieldsTop->EyITemp[i][1][k];
//      vectorZ[i][ny-2][k] = injFieldsTop->EzITemp[i][1][k];

        vectorX[i][ny-1][k] = injFieldsTop->ExITemp[i][0][k];
        vectorY[i][ny-1][k] = injFieldsTop->EyITemp[i][0][k];
        vectorZ[i][ny-1][k] = injFieldsTop->EzITemp[i][0][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][1] = injFieldsRear->ExITemp[i][j][1];
//      vectorY[i][j][1] = injFieldsRear->EyITemp[i][j][1];
//      vectorZ[i][j][1] = injFieldsRear->EzITemp[i][j][1];

        vectorX[i][j][0] = injFieldsRear->ExITemp[i][j][0];
        vectorY[i][j][0] = injFieldsRear->EyITemp[i][j][0];
        vectorZ[i][j][0] = injFieldsRear->EzITemp[i][j][0];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
//      vectorX[i][j][nz-2] = injFieldsFront->ExITemp[i][j][1];
//      vectorY[i][j][nz-2] = injFieldsFront->EyITemp[i][j][1];
//      vectorZ[i][j][nz-2] = injFieldsFront->EzITemp[i][j][1];

        vectorX[i][j][nz-1] = injFieldsFront->ExITemp[i][j][0];
        vectorY[i][j][nz-1] = injFieldsFront->EyITemp[i][j][0];
        vectorZ[i][j][nz-1] = injFieldsFront->EzITemp[i][j][0];
      }
  }
}

/*! get Potential array ** */
double ***EMfields3D::getPHI() {
  return (PHI);
}
/*! get Ex(X,Y,Z) */
double &EMfields3D::getEx(int indexX, int indexY, int indexZ) const {
  return (Ex[indexX][indexY][indexZ]);
}
/*! get Electric field component X array */
double ***EMfields3D::getEx() {
  return (Ex);
}
/*! get Electric Field component X array cell without the ghost cells */
double ***EMfields3D::getExc(Grid3DCU *grid) {
  double ***arr;
  double ***tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr3(double,nxc,nyc,nzc);

  grid->interpN2C(tmp, Ex);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[i][j][k];

  delArr3(tmp,nxc,nyc);
  return arr;
}
/*! get Ey(X,Y,Z) */
double &EMfields3D::getEy(int indexX, int indexY, int indexZ) const {
  return (Ey[indexX][indexY][indexZ]);
}
/*! get Electric field component Y array */
double ***EMfields3D::getEy() {
  return (Ey);
}
/*! get Electric Field component Y array cell without the ghost cells */
double ***EMfields3D::getEyc(Grid3DCU *grid) {
  double ***arr;
  double ***tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr3(double,nxc,nyc,nzc);

  grid->interpN2C(tmp, Ey);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[i][j][k];

  delArr3(tmp,nxc,nyc);
  return arr;
}
/*! get Ez(X,Y,Z) */
double &EMfields3D::getEz(int indexX, int indexY, int indexZ) const {
  return (Ez[indexX][indexY][indexZ]);
}
/*! get Electric field component Z array */
double ***EMfields3D::getEz() {
  return (Ez);
}
/*! get Electric Field component Z array cell without the ghost cells */
double ***EMfields3D::getEzc(Grid3DCU *grid) {
  double ***arr;
  double ***tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr3(double,nxc,nyc,nzc);

  grid->interpN2C(tmp, Ez);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[i][j][k];

  delArr3(tmp,nxc,nyc);
  return arr;
}
/*! get Bx(X,Y,Z) */
double &EMfields3D::getBx(int indexX, int indexY, int indexZ) const {
  return (Bxn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component X array */
double ***EMfields3D::getBx() {
  return (Bxn);
}
/*! get Magnetic Field component X array cell without the ghost cells */
double ***EMfields3D::getBxc() {
  double ***arr;
  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=Bxc[i][j][k];
  return arr;
}
/*! get By(X,Y,Z) */
double &EMfields3D::getBy(int indexX, int indexY, int indexZ) const {
  return (Byn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Y array */
double ***EMfields3D::getBy() {
  return (Byn);
}
/*! get Magnetic Field component Y array cell without the ghost cells */
double ***EMfields3D::getByc() {
  double ***arr;
  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=Byc[i][j][k];
  return arr;
}
/*! get Bz(X,Y,Z) */
double &EMfields3D::getBz(int indexX, int indexY, int indexZ) const {
  return (Bzn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Z array */
double ***EMfields3D::getBz() {
  return (Bzn);
}
/*! get Magnetic Field component Z array cell without the ghost cells */
double ***EMfields3D::getBzc() {
  double ***arr;
  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=Bzc[i][j][k];
  return arr;
}
/*! get rhoc(X,Y,Z) */
double &EMfields3D::getRHOc(int indexX, int indexY, int indexZ) const {
  return (rhoc[indexX][indexY][indexZ]);
} double ***EMfields3D::getRHOc() {
  return (rhoc);
}
/*! get density on node(indexX,indexY,indexZ) */
double &EMfields3D::getRHOn(int indexX, int indexY, int indexZ) const {
  return (rhon[indexX][indexY][indexZ]);
}
/*! get density array defined on nodes */
double ***EMfields3D::getRHOn() {
  return (rhon);
}
/*! get rhos(X,Y,Z) : density for species */
double &EMfields3D::getRHOns(int indexX, int indexY, int indexZ, int is) const {
  return (rhons[is][indexX][indexY][indexZ]);
}
/*! SPECIES: get density array defined on center cells */
double &EMfields3D::getRHOcs(int indexX, int indexY, int indexZ, int is) const {
  return (rhocs[is][indexX][indexY][indexZ]);
}
/*! get density array defined on nodes */
double ****EMfields3D::getRHOns() {
  return (rhons);
}
/*! get species density component X array cell without the ghost cells */
double ***EMfields3D::getRHOcs(Grid3DCU *grid, int is) {
  double ***arr;
  double ****tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr4(double,ns,nxc,nyc,nzc);

  grid->interpN2C(tmp, is, rhons);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[is][i][j][k];

  delArr4(tmp,nxc,nyc,nzc);
  return arr;
}

/*! get Bx_ext(X,Y,Z)  */
double &EMfields3D::getBx_ext(int indexX, int indexY, int indexZ) const{
  return(Bx_ext[indexX][indexY][indexZ]);
}
/*!  get By_ext(X,Y,Z) */
double &EMfields3D::getBy_ext(int indexX, int indexY, int indexZ) const{
  return(By_ext[indexX][indexY][indexZ]);
}
/*!  get Bz_ext(X,Y,Z) */
double &EMfields3D::getBz_ext(int indexX, int indexY, int indexZ) const{
  return(Bz_ext[indexX][indexY][indexZ]);
}

/*! get Bx_ext  */
double ***EMfields3D::getBx_ext() {
  return(Bx_ext);
}
/*!  get By_ext */
double ***EMfields3D::getBy_ext() {
  return(By_ext);
}
/*!  get Bz_ext */
double ***EMfields3D::getBz_ext() {
  return(Bz_ext);
}

/*! SPECIES: get pressure tensor component XX defined on nodes */
double ****EMfields3D::getpXXsn() {
  return (pXXsn);
}
/*! SPECIES: get pressure tensor component XY defined on nodes */
double ****EMfields3D::getpXYsn() {
  return (pXYsn);
}
/*! SPECIES: get pressure tensor component XZ defined on nodes */
double ****EMfields3D::getpXZsn() {
  return (pXZsn);
}
/*! SPECIES: get pressure tensor component YY defined on nodes */
double ****EMfields3D::getpYYsn() {
  return (pYYsn);
}
/*! SPECIES: get pressure tensor component YZ defined on nodes */
double ****EMfields3D::getpYZsn() {
  return (pYZsn);
}
/*! SPECIES: get pressure tensor component ZZ defined on nodes */
double ****EMfields3D::getpZZsn() {
  return (pZZsn);
}
/*! get current -Direction X */
double &EMfields3D::getJx(int indexX, int indexY, int indexZ) const {
  return (Jx[indexX][indexY][indexZ]);
}
/*! get current array X component * */
double ***EMfields3D::getJx() {
  return (Jx);
}
/*! get current -Direction Y */
double &EMfields3D::getJy(int indexX, int indexY, int indexZ) const {
  return (Jy[indexX][indexY][indexZ]);
}
/*! get current array Y component * */
double ***EMfields3D::getJy() {
  return (Jy);
}
/*! get current -Direction Z */
double &EMfields3D::getJz(int indexX, int indexY, int indexZ) const {
  return (Jz[indexX][indexY][indexZ]);
}
/*! get current array Z component * */
double ***EMfields3D::getJz() {
  return (Jz);
}
/*!SPECIES: get current array X component */
double ****EMfields3D::getJxs() {
  return (Jxs);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJxs(int indexX, int indexY, int indexZ, int is) const {
  return (Jxs[is][indexX][indexY][indexZ]);
}
/*! get Magnetic Field component X array species is cell without the ghost cells */
double ***EMfields3D::getJxsc(Grid3DCU *grid, int is) {
  double ***arr;
  double ****tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr4(double,ns,nxc,nyc,nzc);

  grid->interpN2C(tmp, is, Jxs);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[is][i][j][k];

  delArr4(tmp,nxc,nyc,nzc);
  return arr;
}
/*! SPECIES: get current array Y component */
double ****EMfields3D::getJys() {
  return (Jys);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJys(int indexX, int indexY, int indexZ, int is) const {
  return (Jys[is][indexX][indexY][indexZ]);
}
/*! get current component Y array species is cell without the ghost cells */
double ***EMfields3D::getJysc(Grid3DCU *grid, int is) {
  double ***arr;
  double ****tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr4(double,ns,nxc,nyc,nzc);

  grid->interpN2C(tmp, is, Jys);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[is][i][j][k];

  delArr4(tmp,nxc,nyc,nzc);
  return arr;
}
/*!SPECIES: get current array Z component */
double ****EMfields3D::getJzs() {
  return (Jzs);
}
/*! get Jxs(X,Y,Z,is) : density for species */
double &EMfields3D::getJzs(int indexX, int indexY, int indexZ, int is) const {
  return (Jzs[is][indexX][indexY][indexZ]);
}
/*! get current component Z array species is cell without the ghost cells */
double ***EMfields3D::getJzsc(Grid3DCU *grid, int is) {
  double ***arr;
  double ****tmp;

  arr = newArr3(double,nxc-2,nyc-2,nzc-2);
  tmp = newArr4(double,ns,nxc,nyc,nzc);

  grid->interpN2C(tmp, is, Jzs);

  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++)
        arr[i-1][j-1][k-1]=tmp[is][i][j][k];

  delArr4(tmp,nxc,nyc,nzc);
  return arr;
}
/*! get the electric field energy */
double EMfields3D::getEenergy(void) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++)
        localEenergy += .5 * dx * dy * dz * (Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI);

  MPI_Allreduce(&localEenergy, &totalEenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalEenergy);

}
/*! get the magnetic field energy */
double EMfields3D::getBenergy(void) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++){
        Bxt = Bxn[i][j][k]+Bx_ext[i][j][k];
        Byt = Byn[i][j][k]+By_ext[i][j][k];
        Bzt = Bzn[i][j][k]+Bz_ext[i][j][k];
        localBenergy += .5*dx*dy*dz*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalBenergy);
}


/*! Print info about electromagnetic field */
void EMfields3D::print(void) const {
}

/*! destructor: deallocate arrays */
EMfields3D::~EMfields3D() {
  // nodes
  delArr3(Ex, nxn, nyn);
  delArr3(Ey, nxn, nyn);
  delArr3(Ez, nxn, nyn);
  delArr3(Exth, nxn, nyn);
  delArr3(Eyth, nxn, nyn);
  delArr3(Ezth, nxn, nyn);
  delArr3(Bxn, nxn, nyn);
  delArr3(Byn, nxn, nyn);
  delArr3(Bzn, nxn, nyn);
  delArr3(rhon, nxn, nyn);
  delArr3(Jx, nxn, nyn);
  delArr3(Jy, nxn, nyn);
  delArr3(Jz, nxn, nyn);
  delArr3(Jxh, nxn, nyn);
  delArr3(Jyh, nxn, nyn);
  delArr3(Jzh, nxn, nyn);
  // nodes and species
  delArr4(rhons, ns, nxn, nyn);
  delArr4(Jxs, ns, nxn, nyn);
  delArr4(Jys, ns, nxn, nyn);
  delArr4(Jzs, ns, nxn, nyn);
  delArr4(pXXsn, ns, nxn, nyn);
  delArr4(pXYsn, ns, nxn, nyn);
  delArr4(pXZsn, ns, nxn, nyn);
  delArr4(pYYsn, ns, nxn, nyn);
  delArr4(pYZsn, ns, nxn, nyn);
  delArr4(pZZsn, ns, nxn, nyn);
  // central points
  delArr3(PHI, nxc, nyc);
  delArr3(Bxc, nxc, nyc);
  delArr3(Byc, nxc, nyc);
  delArr3(Bzc, nxc, nyc);
  delArr3(rhoc, nxc, nyc);
  delArr3(rhoh, nxc, nyc);
  // various stuff needs to be deallocated too
  delArr3(tempXC, nxc, nyc);
  delArr3(tempYC, nxc, nyc);
  delArr3(tempZC, nxc, nyc);
  delArr3(tempXN, nxn, nyn);
  delArr3(tempYN, nxn, nyn);
  delArr3(tempZN, nxn, nyn);
  delArr3(tempC, nxc, nyc);
  delArr3(tempX, nxn, nyn);
  delArr3(tempY, nxn, nyn);
  delArr3(tempZ, nxn, nyn);
  delArr3(temp2X, nxn, nyn);
  delArr3(temp2Y, nxn, nyn);
  delArr3(temp2Z, nxn, nyn);
  delArr3(imageX, nxn, nyn);
  delArr3(imageY, nxn, nyn);
  delArr3(imageZ, nxn, nyn);
  delArr3(Dx, nxn, nyn);
  delArr3(Dy, nxn, nyn);
  delArr3(Dz, nxn, nyn);
  delArr3(vectX, nxn, nyn);
  delArr3(vectY, nxn, nyn);
  delArr3(vectZ, nxn, nyn);
  delArr3(divC, nxc, nyc);
  for(int i=0;i<sizeMomentsArray;i++)
  {
    delete momentsArray[i];
  }
  delete [] momentsArray;
}
