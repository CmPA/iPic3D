close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/giovannilapenta/Desktop/ddd/nj17/vtk/';

NNcyc=27
for cycle=NNcyc*1000:1000:NNcyc*1000

leggo=0;
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
%Ti=(Pixx+Piyy+Pizz)./rhoi;

% ir=1:Nx;
% jr=1:Ny;
% agyro

Jnorm = 1.2008e-05;
Enorm = 480.7958*25/1836*1000; %mass ratio correction plus mV/mdir
Bnorm = 1.6038e-06;
end

dish=6.5;
%dish=1

% nj12
Lx=114/4;Ly=114/2;Lz=1
xr=[10 15];
yr=[20 36];
xr=[11.5 13];
yr=[26.5 30.5];
xr=[12 12.5];
yr=[27.8 29];

% %nj16
Xxpt=2.189*dish
Yxpt=3.989*dish
Lx=4*dish
Ly=8*dish
xr=Xxpt+[-.25 .25]*dish
xr=Xxpt+[-.3 .3]*dish
yr=Yxpt+[-1 1]*dish

%xr=Xxpt+2*[-.3 .3]*dish
%yr=Yxpt+2*[-1 1]*dish

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);


ir=round(xr./dx);ir=max(1,min(ir)):min(Ny,max(ir));
jr=round(yr./dy);jr=max(1,min(jr)):min(Ny,max(jr));

X=Lx-X;

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
% Vnorm=7e-2 %nj14

%tmp=common_image_due(X(jr,ir),Y(jr,ir),VexbZ(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{ExB}_Z','VexbZ',[-1 1]*Vnorm, Nsm, 1);

[divPx]=divergence(X,Y,permute(Pexx,[2 1]),permute(Pexy,[2 1 ])); % Valid only because dx=dy=dz
[divPy]=divergence(X,Y,permute(Pexy,[2 1]),permute(Peyy,[2 1 ])); % Valid only because dx=dy=dz
[divPz]=divergence(X,Y,permute(Pexz,[2 1]),permute(Peyz,[2 1 ])); % Valid only because dx=dy=dz

divPex=permute(divPx,[2 1]);
divPey=permute(divPy,[2 1]);
divPez=permute(divPz,[2 1]);

VddX = -(divPey.*Bz - divPez.*By)./B.^2./rhoe;
VddY = -(divPez.*Bx - divPex.*Bz)./B.^2./rhoe;
VddZ = -(divPex.*By - divPey.*Bx)./B.^2./rhoe;


[divPx]=divergence(X,Y,permute(Pixx,[2 1]),permute(Pixy,[2 1 ])); % Valid only because dx=dy=dz
[divPy]=divergence(X,Y,permute(Pixy,[2 1]),permute(Piyy,[2 1 ])); % Valid only because dx=dy=dz
[divPz]=divergence(X,Y,permute(Pixz,[2 1]),permute(Piyz,[2 1 ])); % Valid only because dx=dy=dz

divPix=permute(divPx,[2 1]);
divPiy=permute(divPy,[2 1]);
divPiz=permute(divPz,[2 1]);

ViddX = -(divPiy.*Bz - divPiz.*By)./B.^2./rhoi;
ViddY = -(divPiz.*Bx - divPix.*Bz)./B.^2./rhoi;
ViddZ = -(divPix.*By - divPiy.*Bx)./B.^2./rhoi;

[gradBx,gradBy]=gradient(permute(B,[2 1]),dx,dy); % Valid only because dx=dy=dz

gradBx=permute(gradBx,[2 1]);
gradBy=permute(gradBy,[2 1]);
VgBX = ( - Bz.*gradBy)./B.^2;
VgBY = (Bz.*gradBx )./B.^2;
VgBZ = (Bx.*gradBy - By.*gradBx)./B.^2;

VgBZ = -Te./B.*(Bx.*gradBy - By.*gradBx)./B.^2;

% tmp=common_image(X(jr,ir),Y(jr,ir),VddZ(ir,jr), Az(ir,jr),'V_{dde}_Z/c','VedZ',[-1 1]*Vnorm/2.5, Nsm, 2);
% tmp=common_image(X(jr,ir),Y(jr,ir),ViddZ(ir,jr), Az(ir,jr),'V_{ddi}_Z/c','VidZ',[-1 1]*Vnorm/2.5, Nsm*2, 3);
% 
color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Veperpz(ir,jr), AAz(ir,jr),'V_{e\perp M}/c','VeM',[-1 1]*Vnorm*0, Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),-Viperpz(ir,jr), AAz(ir,jr),'V_{i\perp M}/c','ViM',[-1 1]*Vnorm/3*0, Nsm, 9);


color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Ex(ir,jr), AAz(ir,jr),'E_{N}','EN',[-1 1]*3e-3, Nsm, 6);
color_choice=3
tmp=common_image(X(jr,ir),Y(jr,ir),-Vex(ir,jr), AAz(ir,jr),'V_{eN}/c','VeN',[-1 1]*Vnorm/7, Nsm, 7);
tmp=common_image(X(jr,ir),Y(jr,ir),-Vix(ir,jr), AAz(ir,jr),'V_{iN}/c','ViN',[-1 1]*Vnorm/7/2, Nsm, 7);
tmp=common_image(X(jr,ir),Y(jr,ir),-Jez(ir,jr)*Jnorm*1e6, AAz(ir,jr),'J_{eM} [\mu A/m^2]','JeM',[-1 1]*0, Nsm, 8);



figure(100)

jxpt=floor(Yxpt/Ly*Ny);
jband=jxpt-5:jxpt+5
%jband=jxpt;

h(1) = subplot(3,1,1); % upper plot
plot(xc(:),-mean(ohmx(:,jband),2)*Enorm,xc(:),-mean(ohmz(:,jband),2)*Enorm*10);xlim([min(xc(ir)),max(xc(ir))]);
title(['Y=' num2str(mean(yc(jband)))])
grid on
legend('OHM_N [mV/m]','10 OHM_M [mV/m]','location','NorthEast')
h(2) = subplot(3,1,2); % lower plot
plot(xc(:),-mean(Vex(:,jband),2)*3e8/1e3,xc(:),-mean(Vix(:,jband),2)*3e8/1e3,xc(:),-mean(VddX(:,jband),2)*3e8/1e3,xc(:),-mean(VexbX(:,jband),2)*3e8/1e3,xc(:),-mean(Ex(:,jband),2)*Enorm*100,xc(:),mean(By(:,jband),2)*Bnorm*1e9*30,'k');
xlim([min(xc(ir)),max(xc(ir))]);
ylim([-.5 .5]*1e4)
legend('V_{eN} [Km/s]','V_{iN} [Km/s]','V_{\nabla pe,N} [Km/s]','V_{E\times B,N} [Km/s]',' E_N*100 [mV/m]','B_L*30 [nT]','location','NorthEast')
linkaxes(h,'x'); % link the axes in x direction (just for convenience)

grid on
h(3) = subplot(3,1,3)
plot(xc(ir),mean(Agyro_aunai(ir,jband),2),xc(ir),mean(Nongyro_swisdak(ir,jband),2),xc(ir),mean(Agyro(ir,jband),2)/3,'k')
xlim([min(xc(ir)),max(xc(ir))]);
legend('Dng','Q^{1/2}','A\0/3','location','NorthEast')
xlabel('N/d_i','fontsize',[12])
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

figure(102)
plot(xc(:),-mean(Vex(:,jband),2)*3e8/1e3,xc(:),-mean(VexbX(:,jband),2)*3e8/1e3,xc(:),-mean(VddX(:,jband),2)*3e8/1e3,xc(:),-mean(VgBX(:,jband),2)*3e8/1e3)
;xlim([min(xc(ir)),max(xc(ir))]);
xlim([min(xc(ir)),max(xc(ir))]);
ylim([-.5 .5]*1e4)
% h=figure(101)
% set(h, 'Position', [167 26 515 776])
% contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')
% hold on
% contour(xc(ir),yc(jr),Az(ir,jr)',30,'k')
% contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')
% 
% print -dpng contours


% figure(102)
% 
% plot(xc(:),mean(VddZ(:,jband),2),xc(:),mean(Vez(:,jband),2),xc(:),mean(VexbZ(:,jband),2));xlim([min(xc(ir)),max(xc(ir))]);
% 
% vtot=VddZ+VgBZ+VexbZ;
% plot(xc(:),mean(VddZ(:,jband),2),xc(:),mean(VgBZ(:,jband),2),'g',xc(:),mean(VexbZ(:,jband),2),xc(:),mean(Veperpz(:,jband),2)); xlim([min(xc(ir)),max(xc(ir))]);
% 
% figure(103)
% plot(xc(:),mean(ViddZ(:,jband),2),xc(:),mean(Viz(:,jband),2),xc(:),mean(VexbZ(:,jband),2));xlim([min(xc(ir)),max(xc(ir))]);
% 
% vtot=ViddZ-VgBZ+VexbZ;
% plot(xc(:),mean(ViddZ(:,jband),2),xc(:),mean(-VgBZ(:,jband),2),'g',xc(:),mean(VexbZ(:,jband),2),xc(:),mean(Viperpz(:,jband),2)); xlim([min(xc(ir)),max(xc(ir))]);


color_choice=-1
%tmp=common_image_due(X(jr,ir),Y(jr,ir),Agyro(ir,jr), Az(ir,jr),Psiz(ir,jr),'A0_e','Agyro_Scudder',[0 .6], Nsm, 51);
%tmp=common_image_due(X(jr,ir),Y(jr,ir),Agyro_aunai(ir,jr), Az(ir,jr),Psiz(ir,jr),'D_{nge}','Agyro_Aunai',[0 .2], Nsm, 52);

tmp=common_image(X(jr,ir),Y(jr,ir),Nongyro_swisdak(ir,jr), AAz(ir,jr),'Q_e','Nongyro',[0 .2], Nsm, 53);

save('Agyro.mat','Agyro','Agyro_aunai','Nongyro_swisdak')
end