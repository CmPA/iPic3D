%
% Script to read in all varibles from a file called "dir"
%

%leggo='h5'; 
if(strcmp(leggo,'vtk'))


[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

% 
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,~,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
B2D=sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(B.*B2D);
perp2y=Bz.*By./(B.*B2D);
perp2z=-B2D./B;
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
 [Qbulkex,Qbulkey,Qbulkez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qbulke',cycle);
 [Qenthex,Qenthey,Qenthez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qenthe',cycle);
 [Qbulkix,Qbulkiy,Qbulkiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qbulki',cycle);
 [Qenthix,Qenthiy,Qenthiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qenthi',cycle);
[UdivPe,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPe',cycle);
[UdivPi,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPi',cycle);
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;

elseif(strcmp(leggo,'gda'))

ncycle = num2str(cycle)    
    
file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);

file=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny,Nz);

file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny,Nz);

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);


% file=[dir 'Pi_per1_cycle' ncycle '.gda'];
% fid= fopen(file,'rb');
% Pper1=fread(fid,'real*8');
% fclose(fid);
% Pper1=reshape(Pper1,Nx,Ny,Nz);

file=[dir 'Pe_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pexx=fread(fid,'real*8');
fclose(fid);
Pexx=reshape(Pexx,Nx,Ny,Nz);

file=[dir 'Pe_xy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pexy=fread(fid,'real*8');
fclose(fid);
Pexy=reshape(Pexy,Nx,Ny,Nz);

file=[dir 'Pe_xz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pexz=fread(fid,'real*8');
fclose(fid);
Pexz=reshape(Pexz,Nx,Ny,Nz);

file=[dir 'Pe_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peyy=fread(fid,'real*8');
fclose(fid);
Peyy=reshape(Peyy,Nx,Ny,Nz);

file=[dir 'Pe_yz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peyz=fread(fid,'real*8');
fclose(fid);
Peyz=reshape(Peyz,Nx,Ny,Nz);

file=[dir 'Pe_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pezz=fread(fid,'real*8');
fclose(fid);
Pezz=reshape(Pezz,Nx,Ny,Nz);


B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
B2D=sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(B.*B2D);
perp2y=Bz.*By./(B.*B2D);
perp2z=-B2D./B;
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
 
elseif(strcmp(leggo,'h5'))
    % the next line is specific for HRmaha3D3
     
    namefile = [case_name '-Fields']
    fn=[dir,namefile,'_',ncycle,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3);
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    Ns = hinfo.GroupHierarchy.Groups.Attributes.Value;

    Bx = hdf5read(fn,'/Step#0/Block/Bx/0/');
    By = hdf5read(fn,'/Step#0/Block/By/0/');
    Bz = hdf5read(fn,'/Step#0/Block/Bz/0/');
    Ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    Ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    Ez = hdf5read(fn,'/Step#0/Block/Ez/0/');
    
    
    if(Ns>2)
    rhoe = hdf5read(fn,'/Step#0/Block/rho_0/0/') + hdf5read(fn,'/Step#0/Block/rho_2/0/');
    Jex = hdf5read(fn,'/Step#0/Block/Jx_0/0/') + hdf5read(fn,'/Step#0/Block/Jx_2/0/');
    Jey = hdf5read(fn,'/Step#0/Block/Jy_0/0/') + hdf5read(fn,'/Step#0/Block/Jy_2/0/');
    Jez = hdf5read(fn,'/Step#0/Block/Jz_0/0/') + hdf5read(fn,'/Step#0/Block/Jz_2/0/');
    Pexx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/') + hdf5read(fn,'/Step#0/Block/Pxx_2/0/');
    Peyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/') + hdf5read(fn,'/Step#0/Block/Pyy_2/0/');
    Pezz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/') + hdf5read(fn,'/Step#0/Block/Pzz_2/0/');
    Pexy = hdf5read(fn,'/Step#0/Block/Pxy_0/0/') + hdf5read(fn,'/Step#0/Block/Pxy_2/0/');    
    Pexz = hdf5read(fn,'/Step#0/Block/Pxz_0/0/') + hdf5read(fn,'/Step#0/Block/Pxz_2/0/');
    Peyz = hdf5read(fn,'/Step#0/Block/Pyz_0/0/') + hdf5read(fn,'/Step#0/Block/Pyz_2/0/');
    Qex = hdf5read(fn,'/Step#0/Block/EFx_0/0/') + hdf5read(fn,'/Step#0/Block/EFx_2/0/');
    Qey = hdf5read(fn,'/Step#0/Block/EFy_0/0/') + hdf5read(fn,'/Step#0/Block/EFy_2/0/');
    Qez = hdf5read(fn,'/Step#0/Block/EFz_0/0/') + hdf5read(fn,'/Step#0/Block/EFz_2/0/');
    else
    rhoe = hdf5read(fn,'/Step#0/Block/rho_0/0/') ;
    Jex = hdf5read(fn,'/Step#0/Block/Jx_0/0/') ;
    Jey = hdf5read(fn,'/Step#0/Block/Jy_0/0/') ;
    Jez = hdf5read(fn,'/Step#0/Block/Jz_0/0/') ;
    Pexx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/') ;
    Peyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/') ;
    Pezz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/') ;
    Pexy = hdf5read(fn,'/Step#0/Block/Pxy_0/0/') ;    
    Pexz = hdf5read(fn,'/Step#0/Block/Pxz_0/0/') ;
    Peyz = hdf5read(fn,'/Step#0/Block/Pyz_0/0/') ;
    Qex = hdf5read(fn,'/Step#0/Block/EFx_0/0/') ;
    Qey = hdf5read(fn,'/Step#0/Block/EFy_0/0/') ;
    Qez = hdf5read(fn,'/Step#0/Block/EFz_0/0/') ;
    end
    
    
    if(Ns>3)
    rhoi = hdf5read(fn,'/Step#0/Block/rho_1/0/') + hdf5read(fn,'/Step#0/Block/rho_3/0/');
    Jix = hdf5read(fn,'/Step#0/Block/Jx_1/0/') + hdf5read(fn,'/Step#0/Block/Jx_3/0/');
    Jiy = hdf5read(fn,'/Step#0/Block/Jy_1/0/') + hdf5read(fn,'/Step#0/Block/Jy_3/0/');
    Jiz = hdf5read(fn,'/Step#0/Block/Jz_1/0/') + hdf5read(fn,'/Step#0/Block/Jz_3/0/');
    Pixx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/') + hdf5read(fn,'/Step#0/Block/Pxx_3/0/');
    Piyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/') + hdf5read(fn,'/Step#0/Block/Pyy_3/0/');
    Pizz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/') + hdf5read(fn,'/Step#0/Block/Pzz_3/0/');
    Pixy = hdf5read(fn,'/Step#0/Block/Pxy_1/0/') + hdf5read(fn,'/Step#0/Block/Pxy_3/0/');    
    Pixz = hdf5read(fn,'/Step#0/Block/Pxz_1/0/') + hdf5read(fn,'/Step#0/Block/Pxz_3/0/');
    Piyz = hdf5read(fn,'/Step#0/Block/Pyz_1/0/') + hdf5read(fn,'/Step#0/Block/Pyz_3/0/');
    Qix = hdf5read(fn,'/Step#0/Block/EFx_1/0/') + hdf5read(fn,'/Step#0/Block/EFx_3/0/');
    Qiy = hdf5read(fn,'/Step#0/Block/EFy_1/0/') + hdf5read(fn,'/Step#0/Block/EFy_3/0/');
    Qiz = hdf5read(fn,'/Step#0/Block/EFz_1/0/') + hdf5read(fn,'/Step#0/Block/EFz_3/0/');
    else
    rhoi = hdf5read(fn,'/Step#0/Block/rho_1/0/') ;
    Jix = hdf5read(fn,'/Step#0/Block/Jx_1/0/') ;
    Jiy = hdf5read(fn,'/Step#0/Block/Jy_1/0/') ;
    Jiz = hdf5read(fn,'/Step#0/Block/Jz_1/0/') ;
    Pixx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/') ;
    Piyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/') ;
    Pizz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/') ;
    Pixy = hdf5read(fn,'/Step#0/Block/Pxy_1/0/') ;    
    Pixz = hdf5read(fn,'/Step#0/Block/Pxz_1/0/') ;
    Piyz = hdf5read(fn,'/Step#0/Block/Pyz_1/0/') ;
    Qix = hdf5read(fn,'/Step#0/Block/EFx_1/0/') ;
    Qiy = hdf5read(fn,'/Step#0/Block/EFy_1/0/') ;
    Qiz = hdf5read(fn,'/Step#0/Block/EFz_1/0/') ;
    end
    
    
    B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
    B2D=sqrt(Bx.^2+By.^2);
    perp2x=Bz.*Bx./(B.*B2D);
    perp2y=Bz.*By./(B.*B2D);
    perp2z=-B2D./B;
    Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

    [Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Pepar,Peper1,Peper2]=compute_pressure(Bx,By,Bz,Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez,rhoe, qom);
    [Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Pipar,Piper1,Piper2]=compute_pressure(Bx,By,Bz,Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz,rhoi, 1.0);
   
  
    [Qenthex,Qenthey,Qenthez,Qbulkex,Qbulkey,Qbulkez,Qhfex,Qhfey,Qhfez] = ...
    compute_energy_fluxes(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Qex,Qey,Qez,Jex,Jey,Jez,rhoe, qom);

    [Qenthix,Qenthiy,Qenthiz,Qbulkix,Qbulkiy,Qbulkiz,Qhfix,Qhfiy,Qhfiz] = ...
    compute_energy_fluxes(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Qix,Qiy,Qiz,Jix,Jiy,Jiz,rhoi, 1.0);

end