close all
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
dir='/data1/gianni/7feb09HR/'; %directory where the files are

for cycle=68000:1000:68000

    ncycle = num2str(cycle,'%06d');
leggo=0; 
if(leggo==1)


[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
%[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

% 
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
% B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
% Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
 
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;
end


[JdotE,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'JedotT',cycle,0);
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
%[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);
Lx=dx*Nx;Ly=dy*Ny;Lz=Nz*dz;

[x,y,z]=meshgrid(0:dx:Lx-dx,0:dy:Ly-dy,0:dz:Lz-dz);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

qom_ele = -256;

bufferX=round(Nx/6);
bufferY=round(Ny/6);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
kr=3:Nz-3;
kr=-5:5
kr=kr+round(Nz/2);
Nsm=0

%
% Electrons
%



% Compute J dot E
Nsm=3
tmp=common_image(X(jr,ir),Y(jr,ir),mean(JdotE(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'JeE','JeE',[-1 1]*1e-10, Nsm, 1);

[Sx,Sy,Sz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Poynting',cycle,0);
[divS,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'divS',cycle,0);
%tmp = divergence(x,y,z,permute(Sx,[2 1 3]), permute(Sy, [2 1 3]), permute(Sz, [2,1,3]));
%tmp = permute(tmp, [2 1 3]);
%savevtk_bin(tmp,[dir 'divS' ncycle '.vtk'],'divS',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(divS(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'divS','divS',[-1 1]*1e-9, Nsm, 2);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Sx(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Sx','Sx',[-1 1]*2e-9, Nsm, 3);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Sy(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Sy','Sy',[-1 1]*2e-9, Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Sz(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Sz','Sz',[-1 1]*2e-9, Nsm, 5);

%[curlEx, curlEy, curlEz] = curl(x,y,z,permute(Ex,[2 1 3]), permute(Ey, [2 1 3]), permute(Ez, [2,1,3]));
%curlEx = permute(curlEx, [2 1 3]);
%curlEy = permute(curlEy, [2 1 3]);
%curlEz = permute(curlEz, [2 1 3]);

%dUb=-curlEx.*Bx - curlEy.*By - curlEz.*Bz;
%savevtk_bin(dUb,[dir 'dUb' ncycle '.vtk'],'dUb',dx,dy,dz,0,0,0);
[dUb,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'dUb',cycle,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(dUb(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'dUb','dUb',[-1 1]*1e-9, Nsm, 6);

JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez)
savevtk_bin(JidotE,[dir 'JidotE' ncycle '.vtk'],'JidotE',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(JidotE(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'JidotE','JidotE',[-1 1]*1e-9, Nsm, 7);

end