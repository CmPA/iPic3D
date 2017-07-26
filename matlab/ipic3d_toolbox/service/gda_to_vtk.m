%
% Must first run preamble to assign variables
%

for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle);

Bx=read_gda(dir,'B_x',cycle,Nx,Ny,Nz);
By=read_gda(dir,'B_y',cycle,Nx,Ny,Nz);
Bz=read_gda(dir,'B_z',cycle,Nx,Ny,Nz);

savevtkvector_bin(Bx,By,Bz,['B_cycle_' ncycle '.vtk'],'B',dx,dy,dz,0.,0.,0.)

Bx=read_gda(dir,'E_x',cycle,Nx,Ny,Nz);
By=read_gda(dir,'E_y',cycle,Nx,Ny,Nz);
Bz=read_gda(dir,'E_z',cycle,Nx,Ny,Nz);

savevtkvector_bin(Bx,By,Bz,['E_cycle_' ncycle '.vtk'],'E',dx,dy,dz,0.,0.,0.)

% Bx=read_gda(dir,'Pe_xx',cycle,Nx,Ny,Nz);
% By=read_gda(dir,'Pe_yy',cycle,Nx,Ny,Nz);
% Bz=read_gda(dir,'Pe_zz',cycle,Nx,Ny,Nz);
% savevtkvector_bin(Bx,By,Bz,['Pdiag_cycle_' ncycle '.vtk'],'Pdiag',dx,dy,dz,0.,0.,0.)
% 
% Bx=read_gda(dir,'Pe_xy',cycle,Nx,Ny,Nz);
% By=read_gda(dir,'Pe_xz',cycle,Nx,Ny,Nz);
% Bz=read_gda(dir,'Pe_yz',cycle,Nx,Ny,Nz);
% savevtkvector_bin(Bx,By,Bz,['Poffdiag_cycle_' ncycle '.vtk'],'Poffdiag',dx,dy,dz,0.,0.,0.)

Bx=read_gda(dir,'Pe_par',cycle,Nx,Ny,Nz);
By=read_gda(dir,'Pe_per1',cycle,Nx,Ny,Nz);
Bz=read_gda(dir,'Pe_per2',cycle,Nx,Ny,Nz);
savevtkvector_bin(Bx,By,Bz,['Pmagcoord_cycle_' ncycle '.vtk'],'Pmc',dx,dy,dz,0.,0.,0.)
end