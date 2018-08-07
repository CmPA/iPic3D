close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/gianni/Desktop/';

cycle=17000;  
ncycle = num2str(cycle,'%06d');

qom=-256;
must_read=false;
if(must_read)
import_h5_binvtk  
end

dx=Lx/(Nx-1);
dy=Ly/(Ny-1);
dz=Lz/(Nz-1);


[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);

qom_ele = -256;


ir=1:round(Nx/4);
jr=round(Ny/2):Ny-1;
%kr=bufferZ:Nz-bufferZ;


VexbX = (Ey.*Bz - Ez.*By)./B.^2;
VexbY = (Ez.*Bx - Ex.*Bz)./B.^2;
VexbZ = (Ex.*By - Ey.*Bx)./B.^2;


Va=B./sqrt(4*pi.*rhoi);

small=1e-10
Vex= Jex./(rhoe-small) ;
Vey= Jey./(rhoe-small) ;
Vez= Jez./(rhoe-small) ;


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

Vmhdx=(Vix+Vex/abs(qom_ele))/(1+1/abs(qom_ele));
Vmhdy=(Viy+Vey/abs(qom_ele))/(1+1/abs(qom_ele));
Vmhdz=(Viz+Vez/abs(qom_ele))/(1+1/abs(qom_ele));

ohmx = Ex + (Vey.*Bz - Vez.*By);
ohmy = Ey + (Vez.*Bx - Vex.*Bz);
ohmz = Ez + (Vex.*By - Vey.*Bx);
ohmpar=(ohmx.*Bx+ohmy.*By+ohmz.*Bz)./B;

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



color_choice=-1
symmetric_color=1

global cyl Nsm_div
Nsm_div=0; cyl=0;

[divPex, divPey, divPez] = compute_divtensor(x,y,z,Pexx,Pexy,Pexz,Peyy,Peyz,Pezz);
[divPix, divPiy, divPiz] = compute_divtensor(x,y,z,Pixx,Pixy,Pixz,Piyy,Piyz,Pizz);



AAz=vecpot(xc,yc,-mean(Bx(:,:,1),3),mean(By(:,:,1),3));
Vnorm=5e-6

color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),ohmx(ir,jr,1), AAz(ir,jr),'ohmx','OHMx',[-1 1]*Vnorm, Nsm, 1);
tmp=common_image(X(jr,ir),Y(jr,ir),ohmy(ir,jr,1), AAz(ir,jr),'ohmy','OHMy',[-1 1]*Vnorm, Nsm, 2);
tmp=common_image(X(jr,ir),Y(jr,ir),ohmz(ir,jr,1), AAz(ir,jr),'ohmz','OHMz',[-1 1]*Vnorm, Nsm, 3);
tmp=common_image(X(jr,ir),Y(jr,ir),ohmpar(ir,jr,1), AAz(ir,jr),'ohmpar','OHMpar',[-1 1]*Vnorm, Nsm, 7);

U=divPey./(rhoe-small);
tmp=common_image(X(jr,ir),Y(jr,ir),U(ir,jr,1), AAz(ir,jr),'divPe/en','OHMdivPe',[-1 1]*Vnorm, Nsm, 4);


tmp=common_image(X(jr,ir),Y(jr,ir),Bz(ir,jr)-mean(Bz(:)), AAz(ir,jr),'\delta Bz','dBz',[-1 1]*Vnorm*0, Nsm, 5);
tmp=common_image(X(jr,ir),Y(jr,ir),rhoe(ir,jr)-mean(rhoe(:)), AAz(ir,jr),'\rho_e','rhoe',[-1 1]*Vnorm*0, Nsm, 6);
