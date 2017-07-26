close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/gianni/Desktop/ddd/nj12/vtk/';

NNcyc=27

for cycle=NNcyc*1000:1000:NNcyc*1000

leggo=0;
if(leggo==1)

[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
%[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

[Az,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);

B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

Te=(Pexx+Peyy+Pezz)./(-rhoe);
%Ti=(Pixx+Piyy+Pizz)./rhoi;

ir=1:Nx;
jr=1:Ny;
agyro
end


dish=6.5;
%dish=1

% nj12
Lx=114/4*dish;Ly=114/2*dish;Lz=1*dish;


Xxpt=80.0
Yxpt=184.4

xr=Xxpt+[-.25 .25]*dish
xr=Xxpt+[-.3 .3]*dish
yr=Yxpt+[-1 1]*dish


dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

X=Lx-X;
xr=Lx-xr;

ir=round(xr./dx);ir=max(1,min(ir)):min(Ny,max(ir));
jr=round(yr./dy);jr=max(1,min(jr)):min(Ny,max(jr));
ir=1:Nx;
jr=1:Ny;
color_choice=1
global color_choice symmetric_color labelx labely Ncycle

Nsm=3;
labelx='N/d_i';
labely='L/d_i';
Ncycle =num2str(cycle);

VexbX = (Ey.*Bz - Ez.*By)./B.^2;
VexbY = (Ez.*Bx - Ex.*Bz)./B.^2;
VexbZ = (Ex.*By - Ey.*Bx)./B.^2;


Va=B./sqrt(4*pi.*rhoi);

Vex= Jex./rhoe ;
Vey= Jey./rhoe ;
Vez= Jez./rhoe ;

ohmx = Ex + (Vey.*Bz - Vez.*By);
ohmy = Ey + (Vez.*Bx - Vex.*Bz);
ohmz = Ez + (Vex.*By - Vey.*Bx);

Vpar= (Vex.*Bx+Vey.*By+Vez.*Bz)./B;

Veperpx= Vex - Vpar .*Bx./B;
Veperpy= Vey - Vpar .*By./B;
Veperpz= Vez - Vpar .*Bz./B;


xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=vecpot_uniform(xc,yc,Bx*dy/dx,By);
Ns=3;
Psiz=vecpot_uniform(xc,yc,smooth(Jex./rhoe*dy/dx,Ns),smooth(Jey./rhoe,Ns));


color_choice=-1
symmetric_color=0
% tmp=common_image(X(jr,ir),Y(jr,ir),B(ir,jr), Az(ir,jr),'V_A','Va',[-1 1]*0e-2, Nsm, 3);

color_choice=-1
symmetric_color=1

Vnorm=5e-2 %nj12
% Vnorm=7e-2 %nj14

% tmp=common_image_due(X(jr,ir),Y(jr,ir),VexbZ(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{ExB}_Z','VexbZ',[-1 1]*Vnorm, Nsm, 1);

[divPx]=divergence(X,Y,permute(Pexx,[2 1]),permute(Pexy,[2 1 ])); % Valid only because dx=dy=dz
[divPy]=divergence(X,Y,permute(Pexy,[2 1]),permute(Peyy,[2 1 ])); % Valid only because dx=dy=dz
[divPz]=divergence(X,Y,permute(Pexz,[2 1]),permute(Peyz,[2 1 ])); % Valid only because dx=dy=dz

divPx=permute(divPx,[2 1])./rhoe;
divPy=permute(divPy,[2 1])./rhoe;
divPz=permute(divPz,[2 1])./rhoe;

VddX = -(divPx.*Bz - divPz.*By)./B.^2;
VddY = -(divPz.*Bx - divPx.*Bz)./B.^2;
VddZ = -(divPx.*By - divPy.*Bx)./B.^2;


% tmp=common_image_due(X(jr,ir),Y(jr,ir),VddZ(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{dd}_Z','VdZ',[-1 1]*Vnorm, Nsm, 2);
% tmp=common_image_due(X(jr,ir),Y(jr,ir),Veperpz(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{\perp Z}','VeZ',[-1 1]*Vnorm, Nsm, 4);



% tmp=common_image_due(X(jr,ir),Y(jr,ir),Ex(ir,jr), Az(ir,jr),Psiz(ir,jr),'E_{X}','Ex',[-1 1]*3e-3, Nsm, 6);
% tmp=common_image_due(X(jr,ir),Y(jr,ir),Vex(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{eX}','Vex',[-1 1]*Vnorm/7, Nsm, 7);
% tmp=common_image_due(X(jr,ir),Y(jr,ir),Jez(ir,jr), Az(ir,jr),Psiz(ir,jr),'J_{ez}','Jez',[-1 1]*0, Nsm, 8);



figure(100)

jxpt=floor(Yxpt/Ly*Ny);
jband=jxpt-10:jxpt+10
%jband=jxpt;

h(1) = subplot(3,1,1); % upper plot
plot(xc(:),mean(ohmx(:,jband),2),xc(:),10*mean(ohmz(:,jband),2))
xlim([min(xc(ir)),max(xc(ir))]);
title(['Y=' num2str(mean(yc(jband)))])
grid on
legend('OHM_x','10* OHM_z','location','NorthWest')
h(2) = subplot(3,1,2); % lower plot
plot(xc(:),mean(Vex(:,jband),2),xc(:),mean(Ex(:,jband),2),xc(:),mean(By(:,jband),2)/10);
xlim([min(xc(ir)),max(xc(ir))]);
legend('V_{ex}',' E_x','By/10','location','NorthWest')
linkaxes(h,'x'); % link the axes in x direction (just for convenience)

grid on
h(3) = subplot(3,1,3)
plot(xc(ir),mean(Agyro_aunai(ir,jband),2),xc(ir),mean(Nongyro_swisdak(ir,jband),2),xc(ir),mean(Agyro(ir,jband),2)/3)
xlim([min(xc(ir)),max(xc(ir))]);
legend('Dng','Q^{1/2}','A\0/3','location','NorthWest')

linkaxes(h,'x'); % link the axes in x direction (just for convenience)

set(h(1),'xticklabel',[]);
set(h(2),'xticklabel',[]);
pos=get(h,'position');
bottom=pos{3}(2);
top=pos{1}(2)+pos{1}(4);
plotspace=top-bottom;
pos{3}(4)=plotspace/3;
pos{2}(4)=plotspace/3;
pos{1}(4)=plotspace/3;
pos{1}(2)=bottom+plotspace/3;
pos{2}(2)=bottom+2*plotspace/3;

set(h(1),'position',pos{1});
set(h(2),'position',pos{2});
set(h(3),'position',pos{3});

print -dpng flythroughs

h=figure(101)
set(h, 'Position', [167 26 515 776])
contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')
hold on
contour(xc(ir),yc(jr),Az(ir,jr)',30,'k')
contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')

print -dpng contours


color_choice=-1
% tmp=common_image_due(X(jr,ir),Y(jr,ir),Agyro(ir,jr), Az(ir,jr),Psiz(ir,jr),'A0_e','Agyro_Scudder',[0 .6], Nsm, 51);
% tmp=common_image_due(X(jr,ir),Y(jr,ir),Agyro_aunai(ir,jr), Az(ir,jr),Psiz(ir,jr),'D_{nge}','Agyro_Aunai',[0 .2], Nsm, 52);

tmp=common_image(X(jr,ir),Y(jr,ir),Nongyro_swisdak(ir,jr), Az(ir,jr),'Q_e','Nongyro',[0 .2], Nsm, 53);

ir=1:Nx;
jr=1:Ny;
tmp=common_image(X(jr,ir),Y(jr,ir),Nongyro_swisdak(ir,jr), Az(ir,jr),'Q_e','Nongyro_without',[0 .15], Nsm, 54);

xr=[60 100];
yr=[140 220];
ir=round(xr./dx);ir=max(1,min(ir)):min(Ny,max(ir));
jr=round(yr./dy);jr=max(1,min(jr)):min(Ny,max(jr));
tmp=common_image(X(jr,ir),Y(jr,ir),Nongyro_swisdak(ir,jr), Az(ir,jr),'Q_e','Nongyro_without2',[0 .15], Nsm, 55);

end