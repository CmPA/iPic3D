%
% Energy plots in the XZ plane averaging on the whole range in y
%

close all
addpath(genpath('../')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are

addpath(genpath('../../ipic3d_toolbox'));

must_read=true;
leggo='h5';


sim_name='tred82'
switch sim_name    
    case 'tred74'
tred74;
case_name='GEM';
cycle = 18000;
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
case 'tred77'
    TRED77;
    case_name='GEM';
    cycle = 15000;
    zcode = Lz/2;
case 'AH'
    generic;
    case_name='AH';
    cycle =5000;
    zcode = Lz/2;
case 'HRmaha3D3'
    HRmaha3D3;
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
    ygsm=7.05;%1.8;
    zcode = (ygsm - Ygsmrange(1)) * Lz/(Ygsmrange(2)-Ygsmrange(1));
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
    ygsm=7.05;%1.8;
    zcode = (ygsm - Ygsmrange(1)) * Lz/(Ygsmrange(2)-Ygsmrange(1));
otherwise
    print('no recognised case selected')
end



cyl = 0; % cartesian geometry

% Prepare string
ntime = num2str(cycle,'%06d');
ncycle = num2str(cycle,'%06d');


if(must_read)
import_h5_binvtk   
end

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

bufferX=round(Nx/20);
bufferY=round(Ny/20);

%bufferX=round(Nx/2.5);
%bufferY=round(Ny/2.3);

ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;


radius=1; %radius=4


global color_choice symmetric_color labelx labely labelc reversex reversey Ncycle skip

reversey=1;
symmetric_color=1;
color_choice =3;
switch sim_name
case {'tred81','tred82','tred77','AH'}
    labelx ='x/d_i';
    labely ='y/d_i';
    labelc_flux = 'c.u.';
    labelc_power = 'c.u.';
    signx = 1;
    Wm3 = 1; %4pi is due to the usal division by 4pi of the density
    nWm3 = 1;
    mWm2= 1;
    reversex=0;
otherwise
    labelx ='x/R_E';
    labely ='y/R_E';
    labelc_flux = 'mW/m^2';
    labelc_power = 'nW/m^3';
    signx = -1;
    Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the density
    nWm3 = 1e9*Wm3;
    mWm2= Wm3*code_dp*1e3
    reversex=1;
end

skip=10;

xcoord = gsmx(X(jr,ir));
ycoord = gsmy2z(Y(jr,ir));
    xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=vecpot(xc,yc,mean(Bx,3),mean(By,3));
[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);
radius=0.01;

[cp1, cp2, cp3] = cross_prod(mean(Jex,3), mean(Jey,3), mean(Jez,3), mean(Bx,3), mean(By,3), mean(Bz,3));
labelc = labelc_power;
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMJxB_X',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp2(ir,jr),AAz(ir,jr),['Yavg'],'OHMJxB_Y',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp3(ir,jr),AAz(ir,jr),['Yavg'],'OHMJxB_Z',[-1 1]*0e-10, radius,1);

cp1 = mean(rhoe,3).*mean(Ex,3)
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMJnE_X',[-1 1]*0e-10, radius,1);
cp1 = mean(rhoe,3).*mean(Ey,3)
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMJnE_Y',[-1 1]*0e-10, radius,1);
cp1 = mean(rhoe,3).*mean(Ez,3)
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMJnE_Z',[-1 1]*0e-10, radius,1);

[tx, ty] = gradient(imgaussfilt3(permute(mean(Jex./rhoe,3),[2 1]),radius), dx, dy);
tx=permute(tx,[2 1]);ty=permute(ty,[2 1]);
cp1 = tx.*mean(Jex,3) + ty.*mean(Jey,3);

[tx, ty] = gradient(imgaussfilt3(permute(mean(Jey./rhoe,3),[2 1]),radius), dx, dy);
tx=permute(tx,[2 1]);ty=permute(ty,[2 1]);
cp1 = tx.*mean(Jex,3) + ty.*mean(Jey,3);

[tx, ty] = gradient(imgaussfilt3(permute(mean(Jez./rhoe,3),[2 1]),radius), dx, dy);
tx=permute(tx,[2 1]);ty=permute(ty,[2 1]);
cp1 = tx.*mean(Jex,3) + ty.*mean(Jey,3);

tmp=common_image(xcoord,ycoord,cp1(ir,jr)/qom,AAz(ir,jr),['Yavg'],'OHMinert_X',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp2(ir,jr)/qom,AAz(ir,jr),['Yavg'],'OHMinert_Y',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp3(ir,jr)/qom,AAz(ir,jr),['Yavg'],'OHMinert_Z',[-1 1]*0e-10, radius,1);

[cp1, cp2, cp3] = compute_divtensor(x,y,z,mean(Pexx,3),mean(Pexy,3),mean(Pexz,3),mean(Peyy,3), ...
    mean(Peyz,3),mean(Pezz,3),radius,0);
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMdivP_X',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp2(ir,jr),AAz(ir,jr),['Yavg'],'OHMdivP_Y',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,cp3(ir,jr),AAz(ir,jr),['Yavg'],'OHMdivP_Z',[-1 1]*0e-10, radius,1);


[cp1, cp2, cp3] = cross_prod(fluct(Jex), fluct(Jey), fluct(Jez), fluct(Bx), fluct(By), fluct(Bz));
tmp=common_image(xcoord,ycoord,mean(cp1(ir,jr,:),3),AAz(ir,jr),['Yavg'],'OHMdJxdB_X',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,mean(cp2(ir,jr,:),3),AAz(ir,jr),['Yavg'],'OHMdJxdB_Y',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,mean(cp3(ir,jr,:),3),AAz(ir,jr),['Yavg'],'OHMdJxdB_Z',[-1 1]*0e-10, radius,1);

cp1 = mean(fluct(rhoe).*fluct(Ex),3);
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMdndE_X',[-1 1]*0e-10, radius,1);
cp1 = mean(fluct(rhoe).*fluct(Ey),3);
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMdndE_Y',[-1 1]*0e-10, radius,1);
cp1 = mean(fluct(rhoe).*fluct(Ez),3);
tmp=common_image(xcoord,ycoord,cp1(ir,jr),AAz(ir,jr),['Yavg'],'OHMdndE_Z',[-1 1]*0e-10, radius,1);


[tx, ty, tz] = gradient(imgaussfilt3(permute(fluct(Jex./rhoe),[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
cp1 = tx.*fluct(Jex) + ty.*fluct(Jey) + tz.*fluct(Jez);

[tx, ty, tz] = gradient(imgaussfilt3(permute(fluct(Jey./rhoe),[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
cp2 = tx.*fluct(Jex) + ty.*fluct(Jey) + tz.*fluct(Jez);

[tx, ty, tz] = gradient(imgaussfilt3(permute(fluct(Jez./rhoe),[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
cp3 = tx.*fluct(Jex) + ty.*fluct(Jey) + tz.*fluct(Jez);

tmp=common_image(xcoord,ycoord,mean(cp1(ir,jr,:),3)/qom,AAz(ir,jr),['Yavg'],'OHMdinert_X',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,mean(cp2(ir,jr,:),3)/qom,AAz(ir,jr),['Yavg'],'OHMdinert_Y',[-1 1]*0e-10, radius,1);
tmp=common_image(xcoord,ycoord,mean(cp3(ir,jr,:),3)/qom,AAz(ir,jr),['Yavg'],'OHMdinert_Z',[-1 1]*0e-10, radius,1);


function [dp] = fluct(p)
p_avg=mean(p,3);
[Nx Ny Nz]=size(p);
dp=p;
for k=1:Nz
    dp(:,:,k)=p(:,:,k)-p_avg;
end
end