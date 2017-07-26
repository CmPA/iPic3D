[V,Vx,Vy, Vz,dx,dy,dz]=read_vtk_3d('/data1/gianni/paraview/cazzola_pop/V1_cycle38000.vtk',0);
disp('finished reading')
Vxsm=smooth_ipic2d(Vx);
Vysm=smooth_ipic2d(Vy);
Vzsm=smooth_ipic2d(Vz);
disp('starting writing')
savevtkvector(Vxsm, Vysm,Vzsm, 'V1sm_cycle38000.vtk','V1sm',dx,dy,dz)