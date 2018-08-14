%
% Energy plots in the XZ plane averaging on a partial range in Y (code coordinates)
%

close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are

sim_name='tred77'
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

poynting=true;
electrons=true;
ions=true;
saveVTK=false;


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
divS = divS*nWm3;
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



labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JdotE(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','JE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JdotEp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','JEp',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divS(ir,:,kr),Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','divS',[-1 1]*0e-9, radius, 4);


% The poynting flux is in SI units, W/m^3 so we need multiplication by 1e3
% to have it in mW/m^2
labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Sx(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','Sx',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sy(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Sy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sz(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','Sz',[-1 1]*0e-9, radius, 4);


Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sperp1(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Sperp1',[-1 1]*0e-9, radius, 2);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Sperp2(ir,:,kr)*1e3,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','Sperp2',[-1 1]*0e-9, radius, 2);

end

if(electrons)

labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JedotE(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','JeE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQbulk(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','divQbulke',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQenth(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','divQenthe',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQhf(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','divQhfe',[-1 1]*0e-10, radius,1);

labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qbulkex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','Qbulkex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qbulkey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) ,'AVG_Z','Qbulkez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qenthex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qenthex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qenthey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qenthez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qhfex(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qhfex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfey(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qhfey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfez(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qhfez',[-1 1]*0e-9, radius, 4);


Qenthepar= dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthepar(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qenthepar',[-1 1]*0e-9, radius, 2);
Qentheperp1=(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qentheperp1(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qentheprp1',[-1 1]*0e-9, radius, 2);
Qentheperp2=perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qentheperp2(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr) , 'AVG_Z','Qentheprp2',[-1 1]*0e-9, radius, 2);


Qbulkepar= dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkepar(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','Qbulkepar',[-1 1]*0e-9, radius, 2);
Qbulkeperp1=(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkeperp1(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr),'AVG_Z','Qbulkeprp1',[-1 1]*0e-9, radius, 2);
Qbulkeperp2=perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkeperp2(ir,:,kr)*mWm2,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','Qbulkeprp2',[-1 1]*0e-9, radius, 2);


labelc = labelc_power;

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divUP(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','divUPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','UdivPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),ugradp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','Ugradpe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)-ugradp(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','offUdivPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','PgradVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),pdivv(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','pdivVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)-pdivv(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','offPgradVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','UbulkdivVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','UthdivVe',[-1 1]*0e-9, radius, 2);


tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUbulkDt(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','DUbulkeDt',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUthDt(ir,:,kr)*nWm3,Vex_plane(ir,kr),Vez_plane(ir,kr), 'AVG_Z','DUtheDt',[-1 1]*0e-9, radius, 2);

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
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),JidotE(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr),'AVG_Z','JiE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQbulk(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','divQbulki',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQenth(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr),'AVG_Z','divQenthi',[-1 1]*0e-10, radius,1);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divQhf(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','divQhfi',[-1 1]*0e-10, radius,1);

labelc = labelc_flux;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qbulkix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qbulkix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qbulkiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qbulkiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qenthix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qenthix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qenthiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qenthiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),signx*Qhfix(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qhfix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfiy(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qhfiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qhfiz(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qhfiz',[-1 1]*0e-9, radius, 4);


Qenthipar= dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthipar(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qenthipar',[-1 1]*0e-9, radius, 2);
Qenthiperp1=(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiperp1(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qenthiprp1',[-1 1]*0e-9, radius, 2);
Qenthiperp2=perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qenthiperp2(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qenthiprp2',[-1 1]*0e-9, radius, 2);


Qbulkipar= dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkipar(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) , 'AVG_Z','Qbulkipar',[-1 1]*0e-9, radius, 2);
Qbulkiperp1=(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiperp1(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qbulkiprp1',[-1 1]*0e-9, radius, 2);
Qbulkiperp2=perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Qbulkiperp2(ir,:,kr)*mWm2,Vix_plane(ir,kr),Viz_plane(ir,kr) ,'AVG_Z','Qbulkiprp2',[-1 1]*0e-9, radius, 2);


labelc = labelc_power;
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),divUP(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','divUPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','UdivPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),ugradp(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','Ugradpi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),udivP(ir,:,kr)-ugradp(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','offUdivPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','PgradVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),pdivv(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','pdivVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),PgradV(ir,:,kr)-pdivv(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','offPgradVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Ubulk(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','UbulkdivVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Uth(ir,:,kr).*divVe(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','UthdivVi',[-1 1]*0e-9, radius, 2);

tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUbulkDt(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','DUbulkiDt',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),DUthDt(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','DUthiDt',[-1 1]*0e-9, radius, 2);

end

agyro=true
if(agyro)
    symmetric_color=0;
    color_choice =-1;
    Agyro_aunai=hdf5read(fn,'/Step#0/Block/Agyro_aunai/0/');
    Agyro=hdf5read(fn,'/Step#0/Block/Agyro/0/');
    Nongyro_swisdak=hdf5read(fn,'/Step#0/Block/Nongyro_swisdak/0/');
    Nongyro_aunai=hdf5read(fn,'/Step#0/Block/Nongyro_aunai/0/');
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Agyro(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','Agyro',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Agyro_aunai(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','Agyro-aunai',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Nongyro_swisdak(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','Nongyro-swisdak',[-1 1]*0e-9, radius, 2);
    tmp=common_image_vel_maxp(iy_plane(ir,kr),jr,gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),Nongyro_aunai(ir,:,kr)*nWm3,Vix_plane(ir,kr),Viz_plane(ir,kr), 'AVG_Z','Nongyro-aunai',[-1 1]*0e-9, radius, 2);
end



!/usr/local/bin/convert \( PgradVe.png -trim pdivVe.png -trim offPgradVe.png -trim -append \)  \( UdivPe.png -trim Ugradpe.png -trim offUdivPe.png -trim -append \) \( divUPe.png -trim JeE.png -trim JEp.png -trim -append \) \( Agyro.png -trim Agyro-aunai.png -trim Nongyro-swisdak.png -trim -append \) +append comboe.png

!/usr/local/bin/convert \( PgradVi.png -trim pdivVi.png -trim offPgradVi.png -trim -append \)  \( UdivPi.png -trim Ugradpi.png -trim offUdivPi.png -trim -append \) divUPi.png -trim +append comboi.png

unix('convert \( PgradVe.png -trim pdivVe.png -trim offPgradVe.png -trim -append \)  \( UdivPe.png -trim Ugradpe.png -trim offUdivPe.png -trim -append \) \( divUPe.png -trim JeE.png -trim JEp.png -trim -append \) \( Agyro.png -trim Agyro-aunai.png -trim Nongyro-swisdak.png -trim -append \) +append comboe.png')

unix('convert \( PgradVi.png -trim pdivVi.png -trim offPgradVi.png -trim -append \)  \( UdivPi.png -trim Ugradpi.png -trim offUdivPi.png -trim -append \) divUPi.png -trim +append comboi.png')


