close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are


HRmaha3D3
cycle = 80002  % for h5
%cycle = 80000  % for vtk binary

%BOW25


%TRED77
%cycle=15000
%use for tred60
%for cycle=20010:1000:20010


ncycle = num2str(cycle,'%06d');
poynting=1; ions=1; electrons=1;saveVTK=0;


import_h5_binvtk  



Lx=dx*Nx;Ly=dy*Ny;Lz=Nz*dz;

[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

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
labelx ='x/d_i';
labely ='z/d_i';
labelc = 'mW/m^2';

% Compute J dot E
JedotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);
method='gaussian'
radius=5;
JedotEsm=dot(smooth3(Jex,method,radius),smooth3(Jey,method,radius),smooth3(Jez,method,radius), ...
    smooth3(Ex,method,radius),smooth3(Ey,method,radius),smooth3(Ez,method,radius));


JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JdotE=JedotE+JidotE;


xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the dencity
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3


for iz=round(Nz/2)%135
kr=-5:5
kr=kr+round(iz);
Nsm=5



AAz=vecpot(xc,yc,-mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));

if(poynting)

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);

divS = compute_div(x,y,z,Sx,Sy,Sz)/4/pi;

Sx=Sx*code_E*code_B/mu0;
Sy=Sy*code_E*code_B/mu0;
Sz=Sz*code_E*code_B/mu0;
    
labelc = 'nW/m^3'; range=[-20 20]
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JdotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'JE',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divS(ir,jr,kr),3)*nWm3,AAz(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'divS',range, Nsm,1+iz);


labelc = 'mW/m^2';
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Sx(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sy(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sz',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sz(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sy',[-1 1]*0e-9, Nsm, 4+iz);

%Sperp1=(By.*Sx-Bx.*Sy)./B2D;
%tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sperp1(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sperp1',[-1 1]*0e-9, Nsm, 2+iz);
%Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
%tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sperp2(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sperp2',[-1 1]*0e-9, Nsm, 2+iz);

end

if(electrons)
labelc = 'mW/m^2';
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qbulkex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkex Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkez Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkez',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkey Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkey',[-1 1]*0e-9, Nsm, 4+iz);

divQbulke = compute_div(x,y,z,Qbulkex,Qbulkey,Qbulkez);

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qenthex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthex Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthez Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthez',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthey Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthey',[-1 1]*0e-9, Nsm, 4+iz);

divQenthe = compute_div(x,y,z,Qenthex,Qenthey,Qenthez);

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

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qhfex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfex Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfez Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfez',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfey Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfey',[-1 1]*0e-9, Nsm, 4+iz);

divQhfe = compute_div(x,y,z,Qhfex,Qhfey,Qhfez);

udivPe = compute_udivP(x,y,z,Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez, rhoe);
pe=(Pexx+Peyy+Pezz)/3;
udivPetrace = compute_udivP(x,y,z,pe,pe,pe,0*pe,0*pe,0*pe,Jex,Jey,Jez, rhoe);
udivPedev = compute_udivP(x,y,z,Pexx-pe,Peyy-pe,Pezz-pe,Pexy,Pexz,Peyz,Jex,Jey,Jez, rhoe);

Nsm=10
labelc = 'nW/m^3';labelc = 'nW/m^3'; range=[-10 10]
%tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(UdivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPi Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPi',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPe Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPe',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPetrace(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPetrace Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPetrace',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPedev(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPedev Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPedev',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JedotEsm(ir,jr,kr),3)*nWm3,AAz(ir,jr),['JeE Y=' num2str(gsmz2y(z(1,1,iz)))],'JeE',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQbulke(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQbulke Y=' num2str(gsmz2y(z(1,1,iz)))],'divQbulke',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQenthe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQenthe Y=' num2str(gsmz2y(z(1,1,iz)))],'divQenthe',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQhfe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQhfe Y=' num2str(gsmz2y(z(1,1,iz)))],'divQhfe',range, Nsm,1+iz);
range=[-10 10]
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(-divQhfe(ir,jr,kr)-divQenthe(ir,jr,kr)+udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['dUthe Y=' num2str(gsmz2y(z(1,1,iz)))],'divUthe',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JedotE(ir,jr,kr)-divQbulke(ir,jr,kr)-udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['dUbulke Y=' num2str(gsmz2y(z(1,1,iz)))],'divUbulke',range, Nsm,1+iz);

end

if(ions)
labelc = 'mW/m^2';
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qbulkix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkix Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkix',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkiz Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiz',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkiy Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiy',[-1 1]*0e-9, Nsm, 4+iz);

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qenthix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthix Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthix',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthiz Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiz',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthiy Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiy',[-1 1]*0e-9, Nsm, 4+iz);

divQbulki = compute_div(x,y,z,Qbulkix,Qbulkiy,Qbulkiz);

divQenthi = compute_div(x,y,z,Qenthix,Qenthiy,Qenthiz);

Qenthipar= dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthipar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthi || Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthipar',[-1 1]*0e-9, Nsm, 2+iz);
Qenthiperp1=(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenth \perp_1 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qenthiperp2=perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenth \perp_2 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiprp2',[-1 1]*0e-9, Nsm, 2+iz);


Qbulkipar= dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkipar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulki || Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkipar',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkiperp1=(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulk \perp_1 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkiperp2=perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz;
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulk \perp_2 Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiprp2',[-1 1]*0e-9, Nsm, 2+iz);

tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qhfix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfix Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfix',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfiz Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfiz',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfiy Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfiy',[-1 1]*0e-9, Nsm, 4+iz);

divQhfi = compute_div(x,y,z,Qhfix,Qhfiy,Qhfiz);

udivPi = compute_udivP(x,y,z,Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz, rhoi);
pi=(Pixx+Piyy+Pizz)/3;
udivPitrace = compute_udivP(x,y,z,pi,pi,pi,0*pi,0*pi,0*pi,Jix,Jiy,Jiz, rhoi);
udivPidev = compute_udivP(x,y,z,Pixx-pi,Piyy-pi,Pizz-pi,Pixy,Pixz,Piyz,Jix,Jiy,Jiz, rhoi);

Nsm=10
labelc = 'nW/m^3';range=[-20 20]
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPi Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPi',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPitrace(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPitrace Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPitrace',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivPidev(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPidev Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPidev',range, Nsm, 2+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JidotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['JiE Y=' num2str(gsmz2y(z(1,1,iz)))],'JiE',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQbulki(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQbulki Y=' num2str(gsmz2y(z(1,1,iz)))],'divQbulki',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQenthi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQenthi Y=' num2str(gsmz2y(z(1,1,iz)))],'divQenthi',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQhfi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['divQhfi Y=' num2str(gsmz2y(z(1,1,iz)))],'divQhfi',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(-divQhfi(ir,jr,kr)-divQenthi(ir,jr,kr)+udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['dUthi Y=' num2str(gsmz2y(z(1,1,iz)))],'divUthi',range, Nsm,1+iz);
tmp=common_image(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JidotE(ir,jr,kr)-divQbulki(ir,jr,kr)-udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['dUbulki Y=' num2str(gsmz2y(z(1,1,iz)))],'divUbulki',range, Nsm,1+iz);

end
end

if(saveVTK)
    savevtkvector_bin(Qhfex, Qhfey, Qhfez, [dir 'Qhfe' ncycle '.vtk'],'Qhfe',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qhfix, Qhfiy, Qhfiz, [dir 'Qhfi' ncycle '.vtk'],'Qhfi',dx,dy,dz,0,0,0);
    savevtk_bin(udivPe,[dir 'UdivPe' ncycle '.vtk'],'UdivPe',dx,dy,dz,0,0,0);
    savevtk_bin(udivPi,[dir 'UdivPi' ncycle '.vtk'],'UdivPi',dx,dy,dz,0,0,0);
end

