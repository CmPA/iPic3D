close all

!rm *.png
addpath(genpath('../../ipic3d_toolbox'));

must_read=true;
leggo='h5';
if(must_read)

sim_name='tred81'
switch sim_name
case 'tred77'
TRED77;
case_name='GEM';
cycle = 15000;
zcode = Lz/2;
case 'tred81'
tred81;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;
    case 'tred82'
tred82;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;
case 'AH'
generic;
case_name='AH';
cycle =4000;
zcode = Lz/2;
case 'HRmaha3D3'
HRmaha3D3;
leggo='gda';
    case_name='GEM';
dir='/data1/gianni/HRmaha3D3/h5/'; cycle= 80002; ncycle = num2str(cycle,'%06d');
cycle = 80002;  % for h5
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
case '7feb09'
FEB09;
cycle=18000
case_name='MHDUCLA'
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
otherwise
print('no recognised case selected')
end

% Prepare string
ntime = num2str(cycle,'%06d');
ncycle = num2str(cycle,'%06d');


import_h5_binvtk
end

global color_choice symmetric_color labelx labely labelc reversex Ncycle
color_choice = 5;
symmetric_color = 0;
labelx = 'x/d_i'
labely = 'y/d_i'
labelc= ''
reversex = 0;
Ncycle = num2str(cycle);


[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx/4/pi, By/4/pi, Bz/4/pi);
S=sqrt(Sx.^2+Sy.^2+Sz.^2);
[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);
radius=0.01;
divS = compute_div(x,y,z,Sx,Sy,Sz,radius, 0.);
divS=divS/B0^3;

ir=5:Nx-5;
jr=5:Ny-5;
kr=5:Nz-5;

xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=vecpot(xc,yc,mean(Bx,3),mean(By,3));
[X Y] = meshgrid(0:dx:Lx,0:dy:Ly);
%imagesc([0 Lx], [0 Ly],squeeze(log10(moment(divS,2,3)')));colorbar
value=log10(moment(divS,2,3));
tmp= common_image(X(jr,ir),Y(jr,ir),value(ir,jr),AAz(ir,jr) ,'Zavg','ZdivS',[-1 1]*0e-9, radius, 4);
AAz=vecpot(xc,yc,mean(Bx,3),mean(By,3));
value=log10(sqrt(mean(Qex,3).^2+mean(Qey,3).^2+mean(Qez,3).^2));
tmp= common_image(X(jr,ir),Y(jr,ir),value(ir,jr),AAz(ir,jr) ,'Zavg','ZQe',[-1 1]*0e-9, radius, 4);
value = mean(Jex./rhoe,3);
tmp= common_image(X(jr,ir),Y(jr,ir),value(ir,jr),AAz(ir,jr) ,'Zavg','ZVex',[-1 1]*0e-9, radius, 4);
JedotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);
%value = JedotE(:,:,round(Nz/2));
value = mean(JedotE,3);
%value = sign(value).*log10(abs(value)/1e-18+1);
color_choice = 3;
tmp= common_image(X(jr,ir),Y(jr,ir),value(ir,jr),AAz(ir,jr) ,'Zavg','ZJedotE',[-.5 .5]*1e-8, radius, 4);
value=mean((Pexx+Pixx+Peyy+Piyy+Pezz+Pizz)./(B.^2/8/pi),3);
value=real(log10(value));
color_choice = 5;
tmp= common_image(X(jr,ir),Y(jr,ir),value(ir,jr),AAz(ir,jr) ,'Zavg','Zlog_beta',[-1 1]*0e-8, radius, 4);


labely = 'z/d_i'
xc=linspace(0, Lx, Nx);
zc=linspace(0, Lz, Nz);
iy=round(Ny/2);
AAz=vecpot(xc,zc,squeeze(Bx(:,iy,:)),squeeze(Bz(:,iy,:)));
[X Z] = meshgrid(0:dx:Lx,0:dz:Lz);
value=squeeze(log10(moment(divS(:,iy-5:iy+5,:),2,2)));
tmp= common_image(X(kr,ir),Z(kr,ir),value(ir,kr),AAz(ir,kr)*0 ,'Yavg','YdivS',[-1 1]*0e-9, radius, 4);
