%
% Energy plots in the XZ plane averaging on a partial range in Y (code coordinates)
%

close all
unix('rm *.png')
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are

sim_name='7feb09'
switch sim_name
case 'tred77'
TRED77;
case_name='GEM';
cycle = 15000;
zcode = Lz/2;
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
case 'AH'
generic;
case_name='AH';
cycle =4000;
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
case '7feb09'
FEB09;
cycle=18000
%cycle=1
case_name='MHDUCLA'
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
leggo='h5'
otherwise
print('no recognised case selected')
end

poynting=true;
electrons=true;
ions=true;
saveVTK=false;
agyro=false;

cyl = 0; % cartesian geometry

% Prepare string
ntime = num2str(cycle,'%06d');
ncycle = num2str(cycle,'%06d');

must_read=true;
if(must_read)
import_h5_binvtk
end

[X Z] = meshgrid(0:dx:Lx-dx,0:dz:Lz-dz);

bufferX=round(Nx/20);
bufferZ=round(Nz/20);
ir=bufferX:Nx-bufferX;
kr=bufferZ:Nz-bufferZ;

jr=-3:3; % span around the maximum plane for integration

radius=1; %radius=4

e=1.6e-19;

global color_choice symmetric_color labelx labely labelc reversex reversey Ncycle skip

reversey=1;
symmetric_color=1;
color_choice =3;
switch sim_name
case {'tred77','AH'}
labelx ='x/d_i';
labely ='z/d_i';
labelc_flux = 'c.u.';
labelc_power = 'c.u.';
signx = 1;
Wm3 = 1; %4pi is due to the usal division by 4pi of the density
nWm3 = 1;
mWm2= 1;
reversex=0;
otherwise
labelx ='x/R_E';
labely ='z/R_E';
labelc_flux = 'mW/m^2';
labelc_power = 'nW/m^3';
signx = -1;
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the density
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3
reversex=1;
end

skip=10;


% Call the heavy lifting
energy_workhorse

%
% Define plane iy(ix,iz) in cell indexes
%

iy_axis=1:Ny;
for ix=1:Nx
for iz=1:Nz
w=Pipar(ix,:,iz);

%w=1.0./sqrt(Bx(ix,:,iz).^2+1e-10);

iy_plane(ix,iz) = round(sum(w.*iy_axis)./sum(w));

Vex_plane(ix,iz) = signx*Vex(ix,iy_plane(ix,iz),iz);
Vez_plane(ix,iz) = Vez(ix,iy_plane(ix,iz),iz);
Vix_plane(ix,iz) = signx*Vix(ix,iy_plane(ix,iz),iz);
Viz_plane(ix,iz) = Viz(ix,iy_plane(ix,iz),iz);
end
end

%
%Smooth velocities
%

%
%Elm energy
%
if(poynting)
Sx=Sx*code_E*code_B/mu0;
Sy=Sy*code_E*code_B/mu0;
Sz=Sz*code_E*code_B/mu0;
divS = divS*nWm3/4/pi;
% the poynting flux is not in W/m^2 that is in SI unit.
%
% I verified that if instead one computes from code units, then divides it
% by 4pi and rescale it with mWm2 like all other fluxes the result is
% identical.
%

%
% n, J and p from the code need to be multiplied by 4pi and then
% renormalized because of the 4pi division in the code from the MHD density
% By the same token, all particle energy fluxes need the 4 pi
% multiplication, but not the Poynting flux that is based on the fields.
%

[dBx,dBy,dBz]=compute_curl(x,y,z,Ex,Ey,Ez);
dBdt=-(Bx.*dBx+By.*dBy+dBz.*Bz)/4/pi*nWm3;

labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JdotE(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','JE',[-1 1]*1e-1, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JdotEp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','JEp',[-1 1]*1e-1, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divS(ir,:,kr),Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','divS',[-1 1]*1e-1, radius, 4);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),dBdt(ir,:,kr),Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','dW_Bdt',[-1 1]*1e-1,radius, 4);

% The poynting flux is in SI units, W/m^3 so we need multiplication by 1e3
% to have it in mW/m^2
labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Sx(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Sx',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sy(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Sy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sz(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Sz',[-1 1]*0e-9, radius, 4);


Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sperp1(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Sperp1',[-1 1]*0e-9, radius, 2);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sperp2(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Sperp2',[-1 1]*0e-9, radius, 2);

end

if(electrons)

labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JedotE(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','JeE',[-1 1]*6e-2, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQbulk(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','divQbulke',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQenth(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','divQenthe',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQhf(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','divQhfe',[-1 1]*0e-10, radius,1);
divQe = compute_div(x,y,z,Qex,Qey,Qez, radius, cyl);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQe(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','divQe',[-1 1]*6e-2, radius,1);

labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Qex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Qez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qbulkex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Qbulkex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qbulkey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'MAXP','Qbulkez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qenthex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qenthex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qenthey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qenthez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qhfex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qhfex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qhfey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qhfez',[-1 1]*0e-9, radius, 4);


Qenthepar= dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthepar(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qenthepar',[-1 1]*0e-9, radius, 2);
Qentheperp1=(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qentheperp1(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qentheprp1',[-1 1]*0e-9, radius, 2);
Qentheperp2=perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qentheperp2(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'MAXP','Qentheprp2',[-1 1]*0e-9, radius, 2);


Qbulkepar= dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkepar(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','Qbulkepar',[-1 1]*0e-9, radius, 2);
Qbulkeperp1=(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkeperp1(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr),'MAXP','Qbulkeprp1',[-1 1]*0e-9, radius, 2);
Qbulkeperp2=perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkeperp2(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','Qbulkeprp2',[-1 1]*0e-9, radius, 2);


labelc = labelc_power;

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divUP(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','divUPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','UdivPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),ugradp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','Ugradpe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)-ugradp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','offUdivPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','PgradVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),pdivv(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','pdivVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)-pdivv(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','offPgradVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','UbulkdivVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','UthdivVe',[-1 1]*0e-9, radius, 2);


tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUbulkDt(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','DUbulkeDt',[-1 1]*6e-2, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUthDt(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','DUtheDt',[-1 1]*6e-2, radius, 2);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','Ubulke',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','Uthe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),2/3*Uth(ir,:,kr)./abs(rhoe(ir,:,kr))*code_T/e,Vex_plane(ir,kr),Vez_plane(ir,kr), 'MAXP','Te',[-1 1]*0e-9, radius, 2);

end

if(ions)
[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV, ugradp, pdivv, divUP] = compute_energy_balance( ...
    rhoi, Jix, Jiy, Jiz,... 
    Qbulkix, Qbulkiy, Qbulkiz, Qenthix, Qenthiy, Qenthiz, Qhfix, Qhfiy, Qhfiz, ...
    Pixx, Piyy, Pizz, Pixy, Pixz, Piyz, x, y, z, dx, dy, dz, 1.0, radius,cyl);

DUbulkDt = JidotE - Ubulk.*divVi - udivP;
DUthDt = -Uth.*divVi -PgradV;

end

if(ions)

labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JidotE(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr),'MAXP','JiE',[-1 1]*1e-1, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQbulk(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','divQbulki',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQenth(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr),'MAXP','divQenthi',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQhf(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','divQhfi',[-1 1]*0e-10, radius,1);
divQi = compute_div(x,y,z,Qix,Qiy,Qiz, radius, cyl);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQi(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr),'MAXP','divQi',[-1 1]*1e-1, radius,1);

labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qbulkix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qbulkix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qbulkiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qbulkiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qenthix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qenthix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qenthiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qenthiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qhfix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qhfix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qhfiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qhfiz',[-1 1]*0e-9, radius, 4);


Qenthipar= dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthipar(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qenthipar',[-1 1]*0e-9, radius, 2);
Qenthiperp1=(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiperp1(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qenthiprp1',[-1 1]*0e-9, radius, 2);
Qenthiperp2=perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiperp2(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qenthiprp2',[-1 1]*0e-9, radius, 2);


Qbulkipar= dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkipar(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'MAXP','Qbulkipar',[-1 1]*0e-9, radius, 2);
Qbulkiperp1=(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiperp1(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qbulkiprp1',[-1 1]*0e-9, radius, 2);
Qbulkiperp2=perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiperp2(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'MAXP','Qbulkiprp2',[-1 1]*0e-9, radius, 2);


labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divUP(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','divUPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','UdivPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),ugradp(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Ugradpi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)-ugradp(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','offUdivPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','PgradVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),pdivv(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','pdivVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)-pdivv(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','offPgradVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','UbulkdivVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','UthdivVi',[-1 1]*0e-9, radius, 2);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUbulkDt(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','DUbulkiDt',[-1 1]*1e-1, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUthDt(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','DUthiDt',[-1 1]*1e-1, radius, 2);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Ubulki',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Uthi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),2/3.*Uth(ir,:,kr)./(rhoi(ir,:,kr))*code_T/e,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Ti',[-1 1]*0e-9, radius, 2);

end

if(saveVTK)
    savevtkvector_bin(Bx*code_B, By*code_B, Bz*code_B, [dir 'B' ncycle '.vtk'],'B',dx,dy,dz,0,0,0);
    radius=1;
    savevtkvector_bin(imgaussfilt3(Qix,radius)*mWm2, imgaussfilt3(Qiy,radius)*mWm2, imgaussfilt3(Qiz,radius)*mWm2, [dir 'Qi' ncycle '.vtk'],'Qi',dx,dy,dz,0,0,0);
    savevtkvector_bin(imgaussfilt3(Qex,radius)*mWm2, imgaussfilt3(Qey,radius)*mWm2, imgaussfilt3(Qez,radius)*mWm2, [dir 'Qe' ncycle '.vtk'],'Qe',dx,dy,dz,0,0,0);
    savevtk_bin(JedotE*nWm3,[dir 'JedotE' ncycle '.vtk'],'JedotE',dx,dy,dz,0,0,0);
    savevtk_bin(JidotE*nWm3,[dir 'JidotE' ncycle '.vtk'],'JidotE',dx,dy,dz,0,0,0);
    savevtk_bin(divS,[dir 'divS' ncycle '.vtk'],'divS',dx,dy,dz,0,0,0);
    savevtk_bin(divS,[dir 'divQe' ncycle '.vtk'],'divQe',dx,dy,dz,0,0,0);
    savevtk_bin(divS,[dir 'divQi' ncycle '.vtk'],'divQi',dx,dy,dz,0,0,0);
    savevtk_bin(rhoe*code_n,[dir 'rhoe' ncycle '.vtk'],'rhoe',dx,dy,dz,0,0,0);

end

if(agyro)
    
    symmetric_color=0;
    color_choice =-1;
    Agyro_aunai=hdf5read(fn,'/Step#0/Block/Agyro_aunai/0/');
    Agyro=hdf5read(fn,'/Step#0/Block/Agyro/0/');
    Nongyro_swisdak=hdf5read(fn,'/Step#0/Block/Nongyro_swisdak/0/');
    Nongyro_aunai=hdf5read(fn,'/Step#0/Block/Nongyro_aunai/0/');
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Agyro(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Agyro',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Agyro_aunai(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Agyro-aunai',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Nongyro_swisdak(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Nongyro-swisdak',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Nongyro_aunai(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'MAXP','Nongyro-aunai',[-1 1]*0e-9, radius, 2);
end



%!/usr/local/bin/convert \( PgradVe.png -trim pdivVe.png -trim offPgradVe.png -trim -append \)  \( UdivPe.png -trim Ugradpe.png -trim offUdivPe.png -trim -append \) \( divUPe.png -trim JeE.png -trim JEp.png -trim -append \) \( Agyro.png -trim Agyro-aunai.png -trim Nongyro-swisdak.png -trim -append \) +append comboe.png

%!/usr/local/bin/convert \( PgradVi.png -trim pdivVi.png -trim offPgradVi.png -trim -append \)  \( UdivPi.png -trim Ugradpi.png -trim offUdivPi.png -trim -append \) divUPi.png -trim +append comboi.png

unix('convert \( PgradVe.png -trim pdivVe.png -trim offPgradVe.png -trim -append \)  \( UdivPe.png -trim Ugradpe.png -trim offUdivPe.png -trim -append \) \( divUPe.png -trim JeE.png -trim JEp.png -trim -append \) \( Agyro.png -trim Agyro-aunai.png -trim Nongyro-swisdak.png -trim -append \) +append comboe.png')

unix('convert \( PgradVi.png -trim pdivVi.png -trim offPgradVi.png -trim -append \)  \( UdivPi.png -trim Ugradpi.png -trim offUdivPi.png -trim -append \) divUPi.png -trim +append comboi.png')


