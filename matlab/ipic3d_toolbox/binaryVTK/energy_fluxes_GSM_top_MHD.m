close all
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are



leggo=0; 
if(leggo==1)
HRmaha3D3

    nome='~/feb1508iPIC.034800UT.dat'; %Bow1-5
 
 
% Data structure:
%x(RE) y(RE) z(RE) bx(nT) by(nT) bz(nT) vx(km/s) vy(km/s) vz(km/s)
% density(cm-3) pressure(pPa) jx(nA/m2) jy(nA/m2) jz(nA/m2)
fid=fopen(nome);
 
s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);
xmin=x(1); xmax=x(2); dx=x(3); Nx= x(4);
s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);
ymin=x(1); ymax=x(2); dy=x(3); Ny= x(4);
s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);;
zmin=x(1); zmax=x(2); dz=x(3); Nz= x(4);
 
 
 
a=fscanf(fid,'%f',[14 inf])';
 

x=permute(reshape(a(:,1),Nz,Ny,Nx),[3 2 1]);
y=permute(reshape(a(:,2),Nz,Ny,Nx),[3 2 1]);
z=permute(reshape(a(:,3),Nz,Ny,Nx),[3 2 1]);
 
x2=x(:,:,1);
y2=y(:,:,1);
 
 
Vx=permute(reshape(a(:,7),Nz,Ny,Nx),[3 2 1])*1e3;
Vy=permute(reshape(a(:,8),Nz,Ny,Nx),[3 2 1])*1e3;
Vz=permute(reshape(a(:,9),Nz,Ny,Nx),[3 2 1])*1e3;
 
V=sqrt(Vx.^2+Vy.^2+Vz.^2);
 
Vavgx=-squeeze(mean(Vx,2));
Vavgy=squeeze(mean(Vy,2));
Vavgz=squeeze(mean(Vz,2));
 
Bx=permute(reshape(a(:,4),Nz,Ny,Nx),[3 2 1])*1e-9;
By=permute(reshape(a(:,5),Nz,Ny,Nx),[3 2 1])*1e-9;
Bz=permute(reshape(a(:,6),Nz,Ny,Nx),[3 2 1])*1e-9;
B=sqrt(Bx.^2+By.^2+Bz.^2);
B2D=sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(B.*B2D);
perp2y=Bz.*By./(B.*B2D);
perp2z=-B2D./B;
[Ex  Ey Ez] = cross_prod(-Vx, -Vy, -Vz, Bx, By, Bz);
 
 
n=permute(reshape(a(:,10),Nz,Ny,Nx),[3 2 1])*1e6;

p=permute(reshape(a(:,11),Nz,Ny,Nx),[3 2 1])*1e-12;

Jx=permute(reshape(a(:,12),Nz,Ny,Nx),[3 2 1])*1e-9;
Jy=permute(reshape(a(:,13),Nz,Ny,Nx),[3 2 1])*1e-9;
Jz=permute(reshape(a(:,14),Nz,Ny,Nx),[3 2 1])*1e-9;
end  



bufferX=round(Nx/20);
bufferY=round(Ny/20);
bufferZ=round(Nz/20);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
kr=bufferZ:Nz-bufferZ;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;
[X Z] = meshgrid(0:dx:Lx-dx,0:dz:Lz-dz);

global color_choice symmetric_color labelx labely labelc reversex reversey Ncycle skip
reversex=1;
reversey=1;
symmetric_color=1;
color_choice =3;
labelx ='x/R_E';
labely ='y/R_E';
labelc = 'mW/m^2';

% Compute J dot E
JdotE=dot(Jx,Jy,Jz,Ex,Ey,Ez);



[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
Sx=Sx/mu0;
SZ=Sy/mu0;
Sz=Sz/mu0;

xc=Lx-linspace(0, Lx, Nx);
zc=linspace(0, Lz, Nz);
Wm3 = 1; %4pi is due to the usal division by 4pi of the dencity
nWm3 = 1e9*Wm3;
mWm2= 1e3;


Nsm=5
skip=3
iz=1
poynting=0
if(poynting)

labelc = 'nW/m^3';
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vavgx(ir,kr),Vavgz(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),squeeze(mean(JdotE(ir,jr,kr),2))*nWm3,Vavgx(ir,kr),Vavgz(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1+iz);

labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Sx(ir,jr,kr),2)*1e3,Vavgx(ir,kr),Vavgz(ir,kr) ,['Sx Z=' 'AVG_Z'],'Sx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sy(ir,jr,kr),2)*1e3,Vavgx(ir,kr),Vavgz(ir,kr) ,['Sz Z=' 'AVG_Z'],'Sy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sz(ir,jr,kr),2)*1e3,Vavgx(ir,kr),Vavgz(ir,kr) ,['Sy Z=' 'AVG_Z'],'Sz',[-1 1]*0e-9, Nsm, 4+iz);

%Spar= dot(Sx,Sy,Sz,Bx,By,Bz)./B;
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Spar(ir,jr,kr),2)*1e3,AAz(ir,kr) ,['S_{||} Z=' 'AVG_Z'],'Spar',[-1 1]*0e-9, Nsm, 2+iz);
Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp1(ir,jr,kr),2)*1e3,Vavgx(ir,kr),Vavgz(ir,kr) ,['S \perp_1 Z=' 'AVG_Z'],'Sperp1',[-1 1]*0e-9, Nsm, 2+iz);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp2(ir,jr,kr),2)*1e3,Vavgx(ir,kr),Vavgz(ir,kr) ,['S \perp_2 Z=' 'AVG_Z'],'Sperp2',[-1 1]*0e-9, Nsm, 2+iz);

end

ions=1
if(ions)
labelc = 'mW/m^2'; Nsm=5;

mp= 1.67e-27
Ubulk=.5*n*mp.*(Vx.^2 +Vy.^2 +Vz.^2);
Qbulkx = Ubulk.*Vx;
Qbulky = Ubulk.*Vy;
Qbulkz = Ubulk.*Vz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qbulkx(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulkx Z=' 'AVG_Z'],'Qbulkx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulky(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulkz Z=' 'AVG_Z'],'Qbulky',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkz(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulky Z=' 'AVG_Z'],'Qbulkz',[-1 1]*0e-9, Nsm, 4+iz);

Uth = 3*p/2; %It should be gamma/(gamma-1)*p
Qenthx = Uth.*Vx;
Qenthy = Uth.*Vy;
Qenthz = Uth.*Vz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qenthx(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenthx Z=' 'AVG_Z'],'Qenthx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthy(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenthz Z=' 'AVG_Z'],'Qenthy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthz(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenthy Z=' 'AVG_Z'],'Qenthz',[-1 1]*0e-9, Nsm, 4+iz);


Qenthpar= dot(Qenthx,Qenthy,Qenthz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthpar(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenth || Z=' 'AVG_Z'],'Qenthpar',[-1 1]*0e-9, Nsm, 2+iz);
Qenthperp1=(By.*Qenthx-Bx.*Qenthy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthperp1(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenth \perp_1 Z=' 'AVG_Z'],'Qenthprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qenthperp2=perp2x.*Qenthx+perp2y.*Qenthy+perp2z.*Qenthz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthperp2(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qenth \perp_2 Z=' 'AVG_Z'],'Qenthprp2',[-1 1]*0e-9, Nsm, 2+iz);


Qbulkpar= dot(Qbulkx,Qbulky,Qbulkz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkpar(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulk || Z=' 'AVG_Z'],'Qbulkpar',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkperp1=(By.*Qbulkx-Bx.*Qbulky)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkperp1(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulk \perp_1 Z=' 'AVG_Z'],'Qbulkprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkperp2=perp2x.*Qbulkx+perp2y.*Qbulky+perp2z.*Qbulkz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkperp2(ir,jr,kr),2)*mWm2,Vavgx(ir,kr),Vavgz(ir,kr) ,['Qbulk \perp_2 Z=' 'AVG_Z'],'Qbulkprp2',[-1 1]*0e-9, Nsm, 2+iz);


[fx fy fz]=gradient(p,dx*code_dp,dy*code_dp,dz*code_dp);
UgradP = fx.*Vx+fy.*Vy+fz.*Vz;
Nsm=5
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(UgradP(ir,jr,kr),2)*nWm3,Vavgx(ir,kr),Vavgz(ir,kr),['UgradP Z=' 'AVG_Z'],'UgradP',[-1 1]*0e-9, Nsm, 2+iz);

end



