close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/data2/gianni/nj12/vtk/';

for cycle=0:1000:65000

leggo=1;
if(leggo==1)

[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

[Az,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);

B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

Te=(Pexx+Peyy+Pezz)./(-rhoe);
Ti=(Pixx+Piyy+Pizz)./rhoi;
end

Lx=114/4;Ly=114/2;Lz=1

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

xr=[10 15];
yr=[20 36];
ir=round(xr./dx);ir=min(ir):max(ir);
jr=round(yr./dy);jr=min(jr):max(jr);

color_choice=1
global color_choice symmetric_color labelx labely Ncycle
symmetric_color=1
Nsm=5;
labelx='x/d_i';
labely='y/d_i';
Ncycle =num2str(cycle);

phase=1
if(phase==1)
tmp=common_image(X(jr,ir),Y(jr,ir),Epar(ir,jr), Az(ir,jr),'E_{||}','Epar',[-1 1]*3e-4, Nsm, 1);


% [Jz]=curl(X,Y,permute(Bx,[2 1]),permute(By,[2 1 ])); % Valid only because dx=dy=dz
% [Jy,Jx]=gradient(permute(Bz,[2 1]),dx,dy);
% Jy=-permute(Jy,[2 1])/4/pi;
% Jx=permute(Jx,[2 1])/4/pi;
% Jz=permute(Jz,[2 1 ])/4/pi;

[dBdtz]=curl(X,Y,permute(Ex,[2 1]),permute(Ey,[2 1 ])); % Valid only because dx=dy=dz
[dBdty,dBdtx]=gradient(permute(Ez,[2 1]),dx,dy);
dBdty=-permute(dBdty,[2 1]);
dBdtx=permute(dBdtx,[2 1]);
dBdtz=permute(dBdtz,[2 1 ]);

tmp=common_image(X(jr,ir),Y(jr,ir),dBdtx(ir,jr), Az(ir,jr),'dB/dt_X','dBdtX',[-1 1]*3e-4, Nsm, 2);

color_choice=-1
tmp=common_image(X(jr,ir),Y(jr,ir),Peper1(ir,jr)./Pepar(ir,jr), Az(ir,jr),'T_{\perp}/T_{||}','TperpOTpar',[0.5 1.5], Nsm, 3);

tmp=common_image(X(jr,ir),Y(jr,ir),Peper1(ir,jr)./Peper2(ir,jr), Az(ir,jr),'T_{\perp 1}/T_{\perp 2}','Tperp1OTperp2',[0.5 1.5], Nsm, 4);

%!convert +append -trim *47000.png 47k.png

end
phase=2
if(phase==2)
    symmetric_color=0
    color_choice=-1
tmp=common_image(X(jr,ir),Y(jr,ir),Te(ir,jr), Az(ir,jr),'T_{e}','Te',[0 0], Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),Ti(ir,jr), Az(ir,jr),'T_{i}','Ti',[0 0], Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),Jez(ir,jr), Az(ir,jr),'J_{ez}','Jez',[0 0], Nsm, 4);
   
end    

end