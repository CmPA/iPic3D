close all
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are

%HRmaha3D3
BOW25

for cycle=8000:1000:8000

    ncycle = num2str(cycle,'%06d');
leggo=0; 
if(leggo==1)


[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

% 
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
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
[UdivPe,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPe',cycle);
[UdivPi,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPi',cycle);
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;
end


Lx=dx*Nx;Ly=dy*Ny;Lz=Nz*dz;

[x,y,z]=meshgrid(0:dx:Lx-dx,0:dy:Ly-dy,0:dz:Lz-dz);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

qom_ele = -256;

bufferX=round(Nx/20);
bufferY=round(Ny/20);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;



%
% Electrons
%


global color_choice symmetric_color labelx labely labelc reversex Ncycle
reversex=1;
symmetric_color=1;
color_choice =3;
labelx ='x/R_E';
labely ='z/R_E';
labelc = 'mW/m^2';

% Compute J dot E
JdotE=dot(Jex,Jey,Jez,Ex,Ey,Ez)+dot(Jix,Jiy,Jiz,Ex,Ey,Ez);
[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
Sx=Sx*code_E*code_B/mu0;
Sy=Sy*code_E*code_B/mu0;
Sz=Sz*code_E*code_B/mu0;

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the dencity
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3


for iz=135
kr=-5:5
kr=kr+round(iz);
Nsm=5



AAz=vecpot(xc,yc,-mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));
labelc = 'nW/m^3';
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JdotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['JE Y=' num2str(gsmz2y(z(1,1,iz)))],'JE',[-1 1]*0e-10, Nsm,1+iz);


labelc = 'mW/m^2';
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Sx(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Sx Y=' num2str(gsmz2y(z(1,1,iz)))],'Sx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sy(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Sz Y=' num2str(gsmz2y(z(1,1,iz)))],'Sy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sz(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Sy Y=' num2str(gsmz2y(z(1,1,iz)))],'Sz',[-1 1]*0e-9, Nsm, 4+iz);

%Spar= dot(Sx,Sy,Sz,Bx,By,Bz)./B;
%tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Spar(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['S_{||} Y=' num2str(gsmz2y(z(1,1,iz)))],'Spar',[-1 1]*0e-9, Nsm, 2+iz);
Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sperp1(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['S \perp_1 Y=' num2str(gsmz2y(z(1,1,iz)))],'Sperp1',[-1 1]*0e-9, Nsm, 2+iz);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sperp2(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['S \perp_2 Y=' num2str(gsmz2y(z(1,1,iz)))],'Sperp2',[-1 1]*0e-9, Nsm, 2+iz);

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qbulkex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkex Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkez Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkey',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkey Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkez',[-1 1]*0e-9, Nsm, 4+iz);

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qenthex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthex Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthez Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthey',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthey Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthez',[-1 1]*0e-9, Nsm, 4+iz);

Qenthepar= dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthepar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthe || Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthepar',[-1 1]*0e-9, Nsm, 2+iz);
Qentheperp1=(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qentheperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenth \perp_1 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qentheprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qentheperp2=perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qentheperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenth \perp_2 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qentheprp2',[-1 1]*0e-9, Nsm, 2+iz);


Qbulkepar= dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkepar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulke || Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkepar',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkeperp1=(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkeperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulk \perp_1 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkeprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkeperp2=perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkeperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulk \perp_2 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkeprp2',[-1 1]*0e-9, Nsm, 2+iz);

kr=-5:5
kr=kr+round(iz);
Nsm=10
labelc = 'nW/m^3';
%tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(UdivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPi Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPi',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(UdivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPe Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPe',[-1 1]*0e-9, Nsm, 2+iz);

end

end