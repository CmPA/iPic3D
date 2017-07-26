close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/gianni/Desktop/ddd/nj12/vtk/';

NNcyc=63
for cycle=NNcyc*1000:1000:NNcyc*1000

leggo=1;
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
end

% nj12
Lx=114/4;Ly=114/2;Lz=1
xr=[10 15];
yr=[20 36];
xr=[11.5 13];
yr=[26.5 30.5];
xr=[12 12.5];
yr=[27.8 29];

% %nj14
% Lx=4
% Ly=8
% xr=[1.25 1.65]
% yr=[3.7 4.3]

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);


ir=round(xr./dx);ir=min(ir):max(ir);
jr=round(yr./dy);jr=min(jr):max(jr);

color_choice=1
global color_choice symmetric_color labelx labely Ncycle

Nsm=0;
labelx='x/d_i';
labely='y/d_i';
Ncycle =num2str(cycle);

VexbX = (Ey.*Bz - Ez.*By)./B.^2;
VexbY = (Ez.*Bx - Ex.*Bz)./B.^2;
VexbZ = (Ex.*By - Ey.*Bx)./B.^2;

Va=B./sqrt(4*pi.*rhoi);

Vex= Jex./rhoe ;
Vey= Jey./rhoe ;
Vez= Jez./rhoe ;

Vpar= (Vex.*Bx+Vey.*By+Vez.*Bz)./B;

Veperpx= Vex - Vpar .*Bx./B;
Veperpy= Vey - Vpar .*By./B;
Veperpz= Vez - Vpar .*Bz./B;


xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=vecpot_uniform(xc,yc,Bx*dy/dx,By);
Ns=3;
Psiz=vecpot_uniform(xc,yc,smooth(Jex./rhoe*dy/dx,Ns),smooth(Jey./rhoe,Ns));


color_choice=-1
symmetric_color=0
tmp=common_image_due(X(jr,ir),Y(jr,ir),B(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_A','Va',[-1 1]*0e-2, Nsm, 3);

symmetric_color=1

Vnorm=5e-2 %nj12
% Vnorm=7e-2 %nj14

tmp=common_image_due(X(jr,ir),Y(jr,ir),VexbZ(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{ExB}_Z','VexbZ',[-1 1]*Vnorm, Nsm, 1);

[divPx]=divergence(X,Y,permute(Pexx,[2 1]),permute(Pexy,[2 1 ])); % Valid only because dx=dy=dz
[divPy]=divergence(X,Y,permute(Pexy,[2 1]),permute(Peyy,[2 1 ])); % Valid only because dx=dy=dz
[divPz]=divergence(X,Y,permute(Pexz,[2 1]),permute(Peyz,[2 1 ])); % Valid only because dx=dy=dz

divPx=permute(divPx,[2 1])./rhoe;
divPy=permute(divPy,[2 1])./rhoe;
divPz=permute(divPz,[2 1])./rhoe;

VddX = -(divPx.*Bz - divPz.*By)./B.^2;
VddY = -(divPz.*Bx - divPx.*Bz)./B.^2;
VddZ = -(divPx.*By - divPy.*Bx)./B.^2;

color_choice=-1
tmp=common_image_due(X(jr,ir),Y(jr,ir),VddZ(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{dd}_Z','VdZ',[-1 1]*Vnorm, Nsm, 2);
tmp=common_image_due(X(jr,ir),Y(jr,ir),Veperpz(ir,jr), Az(ir,jr),Psiz(ir,jr),'V_{\perp Z}','VeZ',[-1 1]*Vnorm, Nsm, 4);



tmp=common_image_due(X(jr,ir),Y(jr,ir),Ex(ir,jr), Az(ir,jr),Psiz(ir,jr),'E_{X}','Ex',[-1 1]*0, Nsm, 6);


figure(100)
jxpt=1594;
jband=jxpt-5:jxpt+5
%jband=jxpt;

h(1) = subplot(2,1,1); % upper plot
plot(xc(:),mean(Veperpz(:,jband),2),xc(:),mean(VexbZ(:,jband),2),xc(:),mean(VddZ(:,jband),2))
xlim([12 12.5]);
title(['Y=' num2str(mean(yc(jband)))])
legend('V_{e\perp z}','V_{ExBz}','V_{DDz}')
h(2) = subplot(2,1,2); % lower plot
plot(xc(:),mean(Vex(:,jband),2),xc(:),mean(VddX(:,jband)+VexbX(:,jband),2),xc(:),mean(Ex(:,jband),2),xc(:),mean(By(:,jband),2)/10);
xlim([12 12.5]);
legend('V_{ex}','V_{ExBx}+V_{DDx}',' E_x','By/10','location','SouthWest')
grid on

linkaxes(h,'x'); % link the axes in x direction (just for convenience)
xlim([12 12.5]);
set(h(1),'xticklabel',[]);
pos=get(h,'position');
bottom=pos{2}(2);
top=pos{1}(2)+pos{1}(4);
plotspace=top-bottom;
pos{2}(4)=plotspace/2;
pos{1}(4)=plotspace/2;
pos{1}(2)=bottom+plotspace/2;

set(h(1),'position',pos{1});
set(h(2),'position',pos{2});

print -dpng flythroughs

h=figure(101)
set(h, 'Position', [167 26 515 776])
contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')
hold on
contour(xc(ir),yc(jr),Az(ir,jr)',30,'k')
contour(xc(ir),yc(jr),Psiz(ir,jr)',30,'m')

print -dpng contours
end