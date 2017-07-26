close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/gianni/Desktop/BOW23/';

global Lx Xgsmrange
global Ly Zgsmrange
global Lz Ygsmrange
Xgsmrange= [7 15];
Zgsmrange= [-2.5 2.5];
Ygsmrange= [-0.5 0.5];


iz=50;
NNcyc=5
for cycle=NNcyc*1000:1000:NNcyc*1000

leggo=0;
if(leggo==1)

[Bx3d,By3d,Bz3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
Bx=squeeze(Bx3d(:,:,iz));
By=squeeze(By3d(:,:,iz));
Bz=squeeze(Bz3d(:,:,iz));
[Ex3d,Ey3d,Ez3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
Ex=squeeze(Ex3d(:,:,iz));
Ey=squeeze(Ey3d(:,:,iz));
Ez=squeeze(Ez3d(:,:,iz));
[Jex3d,Jey3d,Jez3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
Jex=squeeze(Jex3d(:,:,iz));
Jey=squeeze(Jey3d(:,:,iz));
Jez=squeeze(Jez3d(:,:,iz));
[Jix3d,Jiy3d,Jiz3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);
Jix=squeeze(Jix3d(:,:,iz));
Jiy=squeeze(Jiy3d(:,:,iz));
Jiz=squeeze(Jiz3d(:,:,iz));

[Az3d,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
Az=squeeze(Az3d(:,:,iz));

[rhoe3d,rhoi3d,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
rhoe=squeeze(rhoe3d(:,:,iz));
rhoi=squeeze(rhoi3d(:,:,iz));

[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);

B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

%Te=(Pexx+Peyy+Pezz)./(-rhoe);
%Ti=(Pixx+Piyy+Pizz)./rhoi;

% ir=1:Nx;
% jr=1:Ny;
% agyro

Jnorm = 1.2008e-05;
%Enorm = 480.7958*25/1836*1000; %mass ratio correction plus mV/mdir
Enorm = 480.7958*1000;
Bnorm = 1.6038e-06;
end



Lx=32;Ly=20;Lz=4;


dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;


[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);
X=gsmx(X);
Y=gsmy2z(Y);

ir=1:Nx;
jr=20:Ny-20;

%X=Lx-X;

color_choice=1
global color_choice symmetric_color labelx labely Ncycle

Nsm=3;
labelx='x/R_E';
labely='z/R_E';
Ncycle =num2str(cycle);

VexbX = (Ey.*Bz - Ez.*By)./B.^2;
VexbY = (Ez.*Bx - Ex.*Bz)./B.^2;
VexbZ = (Ex.*By - Ey.*Bx)./B.^2;


Va=B./sqrt(4*pi.*rhoi);

Vex= Jex./rhoe ;
Vey= Jey./rhoe ;
Vez= Jez./rhoe ;


Vix= Jix./rhoi ;
Viy= Jiy./rhoi ;
Viz= Jiz./rhoi ;

B2D=sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(B.*B2D);
perp2y=Bz.*By./(B.*B2D);
perp2z=-B2D./B;


Vepar=(Vex.*Bx+Vey.*By+Vez.*Bz)./(B+1e-10);
Veperp1=(By.*Vex-Bx.*Vey)./B2D;

Veperp2=perp2x.*Vex+perp2y.*Vey+perp2z.*Vez;

Eperp1=(By.*Ex-Bx.*Ey)./B2D;

Eperp2=perp2x.*Ex+perp2y.*Ey+perp2z.*Ez;

Vmhdx=(Vix+Vex/25)/(1+1/25);
Vmhdy=(Viy+Vey/25)/(1+1/25);
Vmhdz=(Viz+Vez/25)/(1+1/25);

ohmx = Ex + (Vey.*Bz - Vez.*By);
ohmy = Ey + (Vez.*Bx - Vex.*Bz);
ohmz = Ez + (Vex.*By - Vey.*Bx);

Epx = Ex + (Vmhdy.*Bz - Vmhdz.*By);
Epy = Ey + (Vmhdz.*Bx - Vmhdx.*Bz);
Epz = Ez + (Vmhdx.*By - Vmhdy.*Bx);

JdotEp=(Jex+Jix).*Epx + (Jey+Jiy).*Epy + (Jez+Jiz).*Epz;

Vpar= (Vex.*Bx+Vey.*By+Vez.*Bz)./B;

Veperpx= Vex - Vpar .*Bx./B;
Veperpy= Vey - Vpar .*By./B;
Veperpz= Vez - Vpar .*Bz./B;


Vipar= (Vix.*Bx+Viy.*By+Viz.*Bz)./B;

Viperpx= Vix - Vipar .*Bx./B;
Viperpy= Viy - Vipar .*By./B;
Viperpz= Viz - Vipar .*Bz./B;


xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=vecpot_uniform(xc,yc,-Bx*dy/dx,By);
Ns=3;
Psiz=vecpot_uniform(xc,yc,smooth(Jex./rhoe*dy/dx,Ns),smooth(Jey./rhoe,Ns));


color_choice=-1
symmetric_color=0
%tmp=common_image_due(X(jr,ir),Y(jr,ir),B(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_A','Va',[-1 1]*0e-2, Nsm, 3);

color_choice=-1
symmetric_color=1

Vnorm=5e-2 %nj12
% 
color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Veperpz(ir,jr), AAz(ir,jr),'V_{e\perp Y}/c','VeY',[-1 1]*Vnorm*0, Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),-Viperpz(ir,jr), AAz(ir,jr),'V_{i\perp Y}/c','ViY',[-1 1]*Vnorm/3*0, Nsm, 9);


color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Ex(ir,jr), AAz(ir,jr),'E_{X}','EX',[-1 1]*3e-3, Nsm, 6);
color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Vex(ir,jr), AAz(ir,jr),'V_{eN}/c','VeX',[-1 1]*Vnorm/3, Nsm, 7);
tmp=common_image(X(jr,ir),Y(jr,ir),-Vix(ir,jr), AAz(ir,jr),'V_{iN}/c','ViX',[-1 1]*Vnorm/3, Nsm, 8);
tmp=common_image(X(jr,ir),Y(jr,ir),-Jez(ir,jr)*Jnorm*1e6, AAz(ir,jr),'J_{eY} [\mu A/m^2]','JeY',[-1 1]*0, Nsm, 9);


Nsm=0
for ix=1:Nx
[Y Z] = meshgrid(0:dy:Ly-dy,0:dz:Lz-dz);
Z=gsmz2y(Z);
Y=gsmy2z(Y);
labely='z/R_E';
labelx='y/R_E';
kr=1:Nz;
jr=20:Ny-20;
Vex3d= Jex3d./rhoe3d ; 
Vex2d=squeeze(Vex3d(ix,:,:));
By2d=squeeze(By3d(ix,:,:));
Bz2d=squeeze(Bz3d(ix,:,:));
zc=linspace(0, Lz, Nz);
yc=linspace(0, Ly, Ny);
clear AAz;
AAz=vecpot_uniform(yc,zc,By2d,Bz2d);
tmp=single_image(Z(kr,jr)',Y(kr,jr)',-Vex2d(jr,kr)',['V_{eN}/c  x/R_E=' num2str(gsmx(ix*dx))],['VeXcutYZ' num2str(ix,'%06i')] ,[-1 1]*0, Nsm, 7);
end

end