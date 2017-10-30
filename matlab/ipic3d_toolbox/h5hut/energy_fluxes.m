close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is

%clear all
leggo=2; poynting=1; ions=1; electrons=1;saveVTK=0;

for cycle=40010

if(leggo==2)   
ncycle = num2str(cycle,'%06d');
dir = '~/Dropbox/Science/san_diego/energy-test/'
Lx=25;Ly=50;Lz=1;
qom=[-40,2.5];
    namefile = 'TwoCoils2D-Fields';
    fn=[dir,namefile,'_',ncycle,'.h5'];


    Bx = hdf5read(fn,'/Step#0/Block/Bx/0/')+hdf5read(fn,'/Step#0/Block/Bx_ext/0/');
    [Nx,Ny,Nz]=size(Bx)
    dx=Lx/(Nx-1);
    dy=Ly/(Ny-1);
    dz=Lz/(Nz-1);

    By = hdf5read(fn,'/Step#0/Block/By/0/')+hdf5read(fn,'/Step#0/Block/By_ext/0/');
    Bz = hdf5read(fn,'/Step#0/Block/Bz/0/')+hdf5read(fn,'/Step#0/Block/Bz_ext/0/');
    
    
    Ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    Ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    Ez = hdf5read(fn,'/Step#0/Block/Ez/0/');
    Jex = hdf5read(fn,'/Step#0/Block/Jx_0/0/');%+hdf5read(fn,'/Step#0/Block/Jx_2/0/');
    Jey = hdf5read(fn,'/Step#0/Block/Jy_0/0/');%+hdf5read(fn,'/Step#0/Block/Jy_2/0/');
    Jez = hdf5read(fn,'/Step#0/Block/Jz_0/0/');%+hdf5read(fn,'/Step#0/Block/Jz_2/0/');
    Jix = hdf5read(fn,'/Step#0/Block/Jx_1/0/');%+hdf5read(fn,'/Step#0/Block/Jx_3/0/');
    Jiy = hdf5read(fn,'/Step#0/Block/Jy_1/0/');%+hdf5read(fn,'/Step#0/Block/Jy_3/0/');
    Jiz = hdf5read(fn,'/Step#0/Block/Jz_1/0/');%+hdf5read(fn,'/Step#0/Block/Jz_3/0/');
    
    rhoe = hdf5read(fn,'/Step#0/Block/rho_0/0/');%+hdf5read(fn,'/Step#0/Block/rho_2/0/');
    rhoi = hdf5read(fn,'/Step#0/Block/rho_1/0/');%+hdf5read(fn,'/Step#0/Block/rho_3/0/');

    Pexx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxx_2/0/');
    Peyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/');%+hdf5read(fn,'/Step#0/Block/Pyy_2/0/');
    Pezz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pzz_2/0/');
    Pexy = hdf5read(fn,'/Step#0/Block/Pxy_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxy_2/0/');    
    Pexz = hdf5read(fn,'/Step#0/Block/Pxz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxz_2/0/');
    Peyz = hdf5read(fn,'/Step#0/Block/Pyz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pyz_2/0/');
    
    Pixx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxx_3/0/');
    Piyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/');%+hdf5read(fn,'/Step#0/Block/Pyy_3/0/');
    Pizz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pzz_3/0/');
    Pixy = hdf5read(fn,'/Step#0/Block/Pxy_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxy_3/0/');    
    Pixz = hdf5read(fn,'/Step#0/Block/Pxz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxz_3/0/');
    Piyz = hdf5read(fn,'/Step#0/Block/Pyz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pyz_3/0/');
    B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
    B2D=sqrt(Bx.^2+By.^2);
    perp2x=Bz.*Bx./(B.*B2D);
    perp2y=Bz.*By./(B.*B2D);
    perp2z=-B2D./B;
    Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

    [Pexx,Peyy,Pezz,Pexy,Pexz,Peyz]=compute_pressure(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez,rhoe, qom(1));
    [Pixx,Piyy,Pizz,Pixy,Pixz,Piyz]=compute_pressure(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz,rhoi, qom(2));
    
    Qex = hdf5read(fn,'/Step#0/Block/EFx_0/0/');%+hdf5read(fn,'/Step#0/Block/EFx_2/0/');
    Qey = hdf5read(fn,'/Step#0/Block/EFy_0/0/');%+hdf5read(fn,'/Step#0/Block/EFy_2/0/');
    Qez = hdf5read(fn,'/Step#0/Block/EFz_0/0/');%+hdf5read(fn,'/Step#0/Block/EFz_2/0/');    
    Qix = hdf5read(fn,'/Step#0/Block/EFx_1/0/');%+hdf5read(fn,'/Step#0/Block/EFx_3/0/');
    Qiy = hdf5read(fn,'/Step#0/Block/EFy_1/0/');%+hdf5read(fn,'/Step#0/Block/EFy_3/0/');
    Qiz = hdf5read(fn,'/Step#0/Block/EFz_1/0/');%+hdf5read(fn,'/Step#0/Block/EFz_3/0/'); 
  
    [Qenthex,Qenthey,Qenthez,Qbulkex,Qbulkey,Qbulkez,Qhfex,Qhfey,Qhfez,Ubulke,Uthe] = ...
    compute_energy_fluxes(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Qex,Qey,Qez,Jex,Jey,Jez,rhoe, qom(1));

    [Qenthix,Qenthiy,Qenthiz,Qbulkix,Qbulkiy,Qbulkiz,Qhfix,Qhfiy,Qhfiz,Ubulki,Uthi] = ...
    compute_energy_fluxes(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Qix,Qiy,Qiz,Jix,Jiy,Jiz,rhoi, qom(2));

end

[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

[X Y] = meshgrid(0:dx:Lx,0:dy:Ly);

qom_ele = qom(1);

bufferX=round(Nx/20);
bufferY=round(Ny/20);
ir=1:Nx-bufferX;
jr=bufferY:Ny-bufferY;



%
% Electrons
%


global color_choice symmetric_color labelx labely labelc reversex Ncycle cyl Nsm_div

cyl=1;
reversex=0;signx=1;
symmetric_color=1;
color_choice =3;
labelx ='r/d_i';
labely ='\zeta/d_i';
labelc = 'mW/m^2';
lablec = 'c.u.';
Nsm=10
Nsm_div=Nsm;

divfactor=permute(x,[2 1 3]);
%divfactor=ones(size(divfactor));

% Compute J dot E
%JedotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);

JedotE=dot(smooth3Dnew(Jex,Nsm),smooth3Dnew(Jey,Nsm),smooth3Dnew(Jez,Nsm), ...
    smooth3Dnew(Ex,Nsm),smooth3Dnew(Ey,Nsm),smooth3Dnew(Ez,Nsm));
JedotE = JedotE.*divfactor;

%JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JidotE=dot(smooth3Dnew(Jix,Nsm),smooth3Dnew(Jiy,Nsm),smooth3Dnew(Jiz,Nsm), ...
    smooth3Dnew(Ex,Nsm),smooth3Dnew(Ey,Nsm),smooth3Dnew(Ez,Nsm));
JidotE = JidotE.*divfactor;

JdotE=JedotE+JidotE;


code_E=1;code_J=1;code_dp=1.;code_B=1;mu0=1;eps0=1;

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the dencity
% nWm3 = 1e9*Wm3;
% mWm2= Wm3*code_dp*1e3
nWm3 = Wm3;
mWm2= Wm3*code_dp;


for iz=1%135
kr=0
kr=kr+round(iz);




AAz=vecpot(xc,yc,signx*mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));

if(poynting)

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);

divS = compute_div(x,y,z,Sx,Sy,Sz)/4/pi;

Sx=Sx*code_E*code_B/4/pi;
Sy=Sy*code_E*code_B/4/pi;
Sz=Sz*code_E*code_B/4/pi;
    
labelc = 'c.u.'; range=[-1 1]*.3e-9
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(JdotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r JE Y=' num2str((z(1,1,iz)))],'JE',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divS(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divS Y=' num2str((z(1,1,iz)))],'divS',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),-mean(divS(ir,jr,kr)+JdotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['dEEMF/dt Y=' num2str((z(1,1,iz)))],'dEEMF_dt',range, Nsm,1+iz);


labelc = 'c.u.';range=[-1 1]*0
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Sx(ir,jr,kr),3),AAz(ir,jr) ,['Sx Y=' num2str((z(1,1,iz)))],'Sx',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Sy(ir,jr,kr),3),AAz(ir,jr) ,['Sy Y=' num2str((z(1,1,iz)))],'Sy',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Sz(ir,jr,kr),3),AAz(ir,jr) ,['Sz Y=' num2str((z(1,1,iz)))],'Sz',range, Nsm, 4+iz);

%Spar= dot(Sx,Sy,Sz,Bx,By,Bz)./B;
%tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Spar(ir,jr,kr),3)*1e3,AAz(ir,jr) ,['S_{||} Y=' num2str((z(1,1,iz)))],'Spar',range, Nsm, 2+iz);
Sperp1=divfactor.*(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Sperp1(ir,jr,kr),3),AAz(ir,jr) ,['r S \perp_1 Y=' num2str((z(1,1,iz)))],'Sperp1',range, Nsm, 2+iz);
Sperp2=divfactor.*(perp2x.*Sx+perp2y.*Sy+perp2z.*Sz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Sperp2(ir,jr,kr),3),AAz(ir,jr) ,['r S \perp_2 Y=' num2str((z(1,1,iz)))],'Sperp2',range, Nsm, 2+iz);

end

if(electrons)
labelc = 'c.u.';
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qbulkex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkex Y=' num2str((z(1,1,iz)))],'Qbulkex',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkey Y=' num2str((z(1,1,iz)))],'Qbulkey',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkez Y=' num2str((z(1,1,iz)))],'Qbulkez',range, Nsm, 4+iz);

color_choice =0;symmetric_color=0;range=[0 1]*0
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Ubulke(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Ubulke Y=' num2str((z(1,1,iz)))],'Ubulke',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Uthe(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Uthe Y=' num2str((z(1,1,iz)))],'Uthe',range, Nsm, 2+iz);

small=1e-10;
Tbe=Ubulke./(rhoe*qom(1)+small);
Te=Uthe./(rhoe/qom(1)+small);
Texx=Pexx./(rhoe/qom(1)+small);
Teyy=Peyy./(rhoe/qom(1)+small);
Tezz=Pezz./(rhoe/qom(1)+small);
range=[0 1]*2e-4
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Tbe(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Tbe Y=' num2str((z(1,1,iz)))],'Tbe',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Te(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Te Y=' num2str((z(1,1,iz)))],'Te',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Texx(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Texx Y=' num2str((z(1,1,iz)))],'Texx',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Teyy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Teyy Y=' num2str((z(1,1,iz)))],'Teyy',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Tezz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Tezz Y=' num2str((z(1,1,iz)))],'Tezz',range, Nsm, 2+iz);
color_choice =3;symmetric_color=1;

divQbulke = compute_div(x,y,z,Qbulkex,Qbulkey,Qbulkez);
range=range*0;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qenthex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthex Y=' num2str((z(1,1,iz)))],'Qenthex',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthey Y=' num2str((z(1,1,iz)))],'Qenthey',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthez Y=' num2str((z(1,1,iz)))],'Qenthez',range, Nsm, 4+iz);

%divQenthe = compute_div(x,y,z,smooth3(Qenthex,method,radius),smooth3(Qenthey,method,radius),smooth3(Qenthez,method,radius));
divQenthe = compute_div(x,y,z,Qenthex,Qenthey,Qenthez);

Qenthepar= divfactor.*dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthepar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenthe || Y=' num2str((z(1,1,iz)))],'Qenthepar',range, Nsm, 2+iz);
Qentheperp1=divfactor.*(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qentheperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenth \perp_1 Y=' num2str((z(1,1,iz)))],'Qentheprp1',range, Nsm, 2+iz);
Qentheperp2=divfactor.*(perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qentheperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenth \perp_2 Y=' num2str((z(1,1,iz)))],'Qentheprp2',range, Nsm, 2+iz);


Qbulkepar= divfactor.*dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkepar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulke || Y=' num2str((z(1,1,iz)))],'Qbulkepar',range, Nsm, 2+iz);
Qbulkeperp1=divfactor.*(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkeperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulk \perp_1 Y=' num2str((z(1,1,iz)))],'Qbulkeprp1',range, Nsm, 2+iz);
Qbulkeperp2=divfactor.*(perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkeperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulk \perp_2 Y=' num2str((z(1,1,iz)))],'Qbulkeprp2',range, Nsm, 2+iz);

tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qhfex(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfex Y=' num2str((z(1,1,iz)))],'Qhfex',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfey(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfey Y=' num2str((z(1,1,iz)))],'Qhfey',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfez(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfez Y=' num2str((z(1,1,iz)))],'Qhfez',range, Nsm, 4+iz);

Qhfepar= divfactor.*dot(Qhfex,Qhfey,Qhfez,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfepar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfe || Y=' num2str((z(1,1,iz)))],'Qhfepar',range, Nsm, 2+iz);
Qhfeperp1=divfactor.*(By.*Qhfex-Bx.*Qhfey)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfeperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfe \perp_1 Y=' num2str((z(1,1,iz)))],'Qhfeprp1',range, Nsm, 2+iz);
Qhfeperp2=divfactor.*(perp2x.*Qhfex+perp2y.*Qhfey+perp2z.*Qhfez);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfeperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfe \perp_2 Y=' num2str((z(1,1,iz)))],'Qhfeprp2',range, Nsm, 2+iz);

divQhfe = compute_div(x,y,z,Qhfex,Qhfey,Qhfez);

udivPe = compute_udivP(x,y,z,Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez, rhoe);
pe=(Pexx+Peyy+Pezz)/3;
udivPetrace = compute_udivP(x,y,z,pe,pe,pe,0*pe,0*pe,0*pe,Jex,Jey,Jez, rhoe);
udivPedev = compute_udivP(x,y,z,Pexx-pe,Peyy-pe,Pezz-pe,Pexy,Pexz,Peyz,Jex,Jey,Jez, rhoe);


labelc = 'c.u.';labelc = 'c.u.'; range=[-1 1]*.3e-9
%tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(UdivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['UdivPi Y=' num2str((z(1,1,iz)))],'UdivPi',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPe Y=' num2str((z(1,1,iz)))],'UdivPe',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPetrace(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPetrace Y=' num2str((z(1,1,iz)))],'UdivPetrace',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPedev(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPedev Y=' num2str((z(1,1,iz)))],'UdivPedev',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(JedotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r JeE Y=' num2str((z(1,1,iz)))],'JeE',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQbulke(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQbulke Y=' num2str((z(1,1,iz)))],'divQbulke',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQenthe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQenthe Y=' num2str((z(1,1,iz)))],'divQenthe',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQhfe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQhfe Y=' num2str((z(1,1,iz)))],'divQhfe',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(-divQhfe(ir,jr,kr)-divQenthe(ir,jr,kr)+udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r dUthe Y=' num2str((z(1,1,iz)))],'dUthe',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(JedotE(ir,jr,kr)-divQbulke(ir,jr,kr)-udivPe(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r dUbulke Y=' num2str((z(1,1,iz)))],'dUbulke',range, Nsm,1+iz);

end

if(ions)
labelc = 'c.u.';range=[-1 1]*0;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qbulkix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkix Y=' num2str((z(1,1,iz)))],'Qbulkix',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkiy Y=' num2str((z(1,1,iz)))],'Qbulkiy',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qbulkiz Y=' num2str((z(1,1,iz)))],'Qbulkiz',range, Nsm, 4+iz);

tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Ubulki(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Ubulki Y=' num2str((z(1,1,iz)))],'Ubulki',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Uthi(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Uthi Y=' num2str((z(1,1,iz)))],'Uthi',range, Nsm, 2+iz);

small=1e-10;
Tbi=Ubulke./(rhoi*qom(2)+small);
Ti=Uthe./(rhoi/qom(2)+small);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Tbi(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Tbi Y=' num2str((z(1,1,iz)))],'Tbi',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Ti(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Ti Y=' num2str((z(1,1,iz)))],'Ti',range, Nsm, 2+iz);



tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qenthix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthix Y=' num2str((z(1,1,iz)))],'Qenthix',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthiy Y=' num2str((z(1,1,iz)))],'Qenthiy',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qenthiz Y=' num2str((z(1,1,iz)))],'Qenthiz',range, Nsm, 4+iz);

divQbulki = compute_div(x,y,z,Qbulkix,Qbulkiy,Qbulkiz);

divQenthi = compute_div(x,y,z,Qenthix,Qenthiy,Qenthiz);

Qenthipar= divfactor.*dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthipar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenthi || Y=' num2str((z(1,1,iz)))],'Qenthipar',range, Nsm, 2+iz);
Qenthiperp1=divfactor.*(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthiperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenth \perp_1 Y=' num2str((z(1,1,iz)))],'Qenthiprp1',range, Nsm, 2+iz);
Qenthiperp2=divfactor.*(perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qenthiperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qenth \perp_2 Y=' num2str((z(1,1,iz)))],'Qenthiprp2',range, Nsm, 2+iz);


Qbulkipar= divfactor.*dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkipar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulki || Y=' num2str((z(1,1,iz)))],'Qbulkipar',range, Nsm, 2+iz);
Qbulkiperp1=divfactor.*(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkiperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulk \perp_1 Y=' num2str((z(1,1,iz)))],'Qbulkiprp1',range, Nsm, 2+iz);
Qbulkiperp2=divfactor.*(perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qbulkiperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qbulk \perp_2 Y=' num2str((z(1,1,iz)))],'Qbulkiprp2',range, Nsm, 2+iz);

tmp=common_image((X(jr,ir)),(Y(jr,ir)),signx*mean(Qhfix(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfix Y=' num2str((z(1,1,iz)))],'Qhfix',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfiy(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfiy Y=' num2str((z(1,1,iz)))],'Qhfiy',range, Nsm, 3+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfiz(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['Qhfiz Y=' num2str((z(1,1,iz)))],'Qhfiz',range, Nsm, 4+iz);

Qhfipar= divfactor.*dot(Qhfix,Qhfiy,Qhfiz,Bx,By,Bz)./B;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfipar(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfi || Y=' num2str((z(1,1,iz)))],'Qhfipar',range, Nsm, 2+iz);
Qhfiperp1=divfactor.*(By.*Qhfix-Bx.*Qhfiy)./B2D;
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfiperp1(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfi \perp_1 Y=' num2str((z(1,1,iz)))],'Qhfiprp1',range, Nsm, 2+iz);
Qhfiperp2=divfactor.*(perp2x.*Qhfix+perp2y.*Qhfiy+perp2z.*Qhfiz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(Qhfiperp2(ir,jr,kr),3)*mWm2,AAz(ir,jr) ,['r Qhfi \perp_2 Y=' num2str((z(1,1,iz)))],'Qhfiprp2',range, Nsm, 2+iz);


divQhfi = compute_div(x,y,z,Qhfix,Qhfiy,Qhfiz);

udivPi = compute_udivP(x,y,z,Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz, rhoi);
pres=(Pixx+Piyy+Pizz)/3;
udivPitrace = compute_udivP(x,y,z,pres,pres,pres,0*pres,0*pres,0*pres,Jix,Jiy,Jiz, rhoi);
udivPidev = compute_udivP(x,y,z,Pixx-pres,Piyy-pres,Pizz-pres,Pixy,Pixz,Piyz,Jix,Jiy,Jiz, rhoi);


labelc = 'c.u.';range=[-1 1]*.3e-9
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPi Y=' num2str((z(1,1,iz)))],'UdivPi',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPitrace(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPitrace Y=' num2str((z(1,1,iz)))],'UdivPitrace',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(udivPidev(ir,jr,kr),3)*nWm3,AAz(ir,jr) ,['r UdivPidev Y=' num2str((z(1,1,iz)))],'UdivPidev',range, Nsm, 2+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(JidotE(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r JiE Y=' num2str((z(1,1,iz)))],'JiE',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQbulki(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQbulki Y=' num2str((z(1,1,iz)))],'divQbulki',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQenthi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQenthi Y=' num2str((z(1,1,iz)))],'divQenthi',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(divQhfi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r divQhfi Y=' num2str((z(1,1,iz)))],'divQhfi',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(-divQhfi(ir,jr,kr)-divQenthi(ir,jr,kr)+udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r dUthi Y=' num2str((z(1,1,iz)))],'dUthi',range, Nsm,1+iz);
tmp=common_image((X(jr,ir)),(Y(jr,ir)),mean(JidotE(ir,jr,kr)-divQbulki(ir,jr,kr)-udivPi(ir,jr,kr),3)*nWm3,AAz(ir,jr),['r dUbulki Y=' num2str((z(1,1,iz)))],'dUbulki',range, Nsm,1+iz);

end
end

if(saveVTK)
    savevtkvector_bin(Qhfex, Qhfey, Qhfez, [dir 'Qhfe' ncycle '.vtk'],'Qhfe',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qhfix, Qhfiy, Qhfiz, [dir 'Qhfi' ncycle '.vtk'],'Qhfi',dx,dy,dz,0,0,0);
    savevtk_bin(udivPe,[dir 'UdivPe' ncycle '.vtk'],'UdivPe',dx,dy,dz,0,0,0);
    savevtk_bin(udivPi,[dir 'UdivPi' ncycle '.vtk'],'UdivPi',dx,dy,dz,0,0,0);
end

end