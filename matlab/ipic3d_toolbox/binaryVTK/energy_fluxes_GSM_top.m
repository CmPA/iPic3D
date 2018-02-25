%
% Energy plots in the XZ plane averaging on the whole range in y
%

close all
addpath(genpath('~/iPic3D-github/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are


HRmaha3D3

dir='/data1/gianni/HRmaha3D3/h5/'; cycle= 80002; ncycle = num2str(cycle,'%06d');

cycle = 80002  % for h5
%cycle = 80000  % for vtk binary

% for HRmaha3D1:
 time=60*(cycle/75000.0*Dt/.125) %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D

%ADD initial time of the RUN
time=time+initial_time%(03*60+48)*60
% Prepare string
ntime = datestr(time/86400,'HH:MM:SS UT')


ncycle = num2str(cycle,'%06d');

import_h5_binvtk   


%Lx=dx*Nx;LZ=dy*Ny;Lz=Nz*dz;

%[x,y,z]=meshgrid(0:dx:Lx-dx,0:dy:Ly-dy,0:dz:Lz-dz);

[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

[X Z] = meshgrid(0:dx:Lx-dx,0:dz:Lz-dz);

qom_ele = -256;

bufferX=round(Nx/20);
bufferY=round(Ny/20);
bufferZ=round(Nz/20);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
kr=bufferZ:Nz-bufferZ;

Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the density
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3


cyl=1 % this means it is cartesian
Nsm_div=5;
global cyl Nsm_div


global color_choice symmetric_color labelx labely labelc reversex reversey Ncycle skip
reversex=1;
reversey=1;
symmetric_color=1;
color_choice =3;
labelx ='x/R_E';
labely ='y/R_E';
labelc = 'mW/m^2';
skip=10


% Compute J dot E
JedotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);
method='gaussian'
radius=5;
JedotEsm=dot(smooth3(Jex,method,radius),smooth3(Jey,method,radius),smooth3(Jez,method,radius), ...
    smooth3(Ex,method,radius),smooth3(Ey,method,radius),smooth3(Ez,method,radius));


JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JdotE=JedotE+JidotE;

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
divS = compute_div(x,y,z,smooth3(Sx,'box',5),smooth3(Sy,'box',5),smooth3(Sz,'box',5));

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
xc=Lx-linspace(0, Lx, Nx);
zc=linspace(0, Lz, Nz);

%
% n, J and p from the code need to be multiplied by 4pi and then
% renormalized because of the 4pi division in the code from the MHD density
% By the same token, all particle energy fluxes need the 4 pi
% multiplication, but not the Poynting flux that is based on the fields.
%
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the dencity
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3

%
% Electrons
%


%for iz=135
%kr=-5:5
%kr=kr+round(iz);
Nsm=5


% Vix=Jix./rhoi;Viz=Jiz./rhoi;
% AAzi=vecpot(xc,zc,-squeeze(mean(Vix(:,jr,:),2)),squeeze(mean(Viz(:,jr,:),2)));

Vex=Jex./rhoe;Vez=Jez./rhoe;
Vex=-smoothbc(squeeze(mean(Vex(:,jr,:),2)),Nsm);
Vez=smoothbc(squeeze(mean(Vez(:,jr,:),2)),Nsm);
AAze=vecpot(xc,zc,Vex,Vez);
divVe = compute_div(x,y,z,Vex,Vey,Vez);

Vix=Jix./rhoi;Viz=Jiz./rhoi;
Vix=-smoothbc(squeeze(mean(Vix(:,jr,:),2)),Nsm);
Viz=smoothbc(squeeze(mean(Viz(:,jr,:),2)),Nsm);
AAze=vecpot(xc,zc,Vix,Viz);
divVi = compute_div(x,y,z,Vix,Viy,Viz);

poynting=1
if(poynting)

labelc = 'nW/m^3';
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divS(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr) ,['divS Z=' 'AVG_Z'],'divS',[-1 1]*0e-9, Nsm, 4);


% The poynting flux is in SI units, W/m^3 so we need multiplication by 1e3
% to have it in mW/m^2
labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Sx(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sx Z=' 'AVG_Z'],'Sx',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sy(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sz Z=' 'AVG_Z'],'Sy',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sz(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sy Z=' 'AVG_Z'],'Sz',[-1 1]*0e-9, Nsm, 4);


%Spar= dot(Sx,Sy,Sz,Bx,By,Bz)./B;
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Spar(ir,jr,kr),2)*1e3,AAz(ir,kr) ,['S_{||} Z=' 'AVG_Z'],'Spar',[-1 1]*0e-9, Nsm, 2);
Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp1(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['S \perp_1 Z=' 'AVG_Z'],'Sperp1',[-1 1]*0e-9, Nsm, 2);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp2(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['S \perp_2 Z=' 'AVG_Z'],'Sperp2',[-1 1]*0e-9, Nsm, 2);

end

electrons=1
if(electrons)


[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV] = compute_energy_balance( ...
    rhoe, Jex, Jey, Jez,... 
    Qbulkex, Qbulkey, Qbulkez, Qenthex, Qenthey, Qenthez, Qhfex, Qhfey, Qhfez, ...
    Pexx, Peyy, Pezz, Pexy, Pexz, Peyz, dx, dy, dz, qom_ele)

labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JedotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JeE Z=' 'AVG_Z'],'JeE',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQbulk(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['divQbulke Z=' 'AVG_Z'],'divQbulke',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQenth(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['divQenthe Z=' 'AVG_Z'],'divQenthe',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQhf(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['divQhfe Z=' 'AVG_Z'],'divQhfe',[-1 1]*0e-10, Nsm,1);

labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qbulkex(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkex Z=' 'AVG_Z'],'Qbulkex',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkey(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkez Z=' 'AVG_Z'],'Qbulkey',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkez(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkey Z=' 'AVG_Z'],'Qbulkez',[-1 1]*0e-9, Nsm, 4);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qenthex(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthex Z=' 'AVG_Z'],'Qenthex',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthey(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthez Z=' 'AVG_Z'],'Qenthey',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthez(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthey Z=' 'AVG_Z'],'Qenthez',[-1 1]*0e-9, Nsm, 4);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qhfex(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qhfex Z=' 'AVG_Z'],'Qhfex',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qhfey(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qhfez Z=' 'AVG_Z'],'Qhfey',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qhfez(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qhfey Z=' 'AVG_Z'],'Qhfez',[-1 1]*0e-9, Nsm, 4);


Qenthepar= dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthepar(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthe || Z=' 'AVG_Z'],'Qenthepar',[-1 1]*0e-9, Nsm, 2);
Qentheperp1=(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qentheperp1(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenth \perp_1 Z=' 'AVG_Z'],'Qentheprp1',[-1 1]*0e-9, Nsm, 2);
Qentheperp2=perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qentheperp2(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenth \perp_2 Z=' 'AVG_Z'],'Qentheprp2',[-1 1]*0e-9, Nsm, 2);


Qbulkepar= dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkepar(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulke || Z=' 'AVG_Z'],'Qbulkepar',[-1 1]*0e-9, Nsm, 2);
Qbulkeperp1=(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkeperp1(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulk \perp_1 Z=' 'AVG_Z'],'Qbulkeprp1',[-1 1]*0e-9, Nsm, 2);
Qbulkeperp2=perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkeperp2(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulk \perp_2 Z=' 'AVG_Z'],'Qbulkeprp2',[-1 1]*0e-9, Nsm, 2);

Nsm=10
labelc = 'nW/m^3';

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(udivP(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UdivPe Z=' 'AVG_Z'],'UdivPe',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(PgradV(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['PgradVe Z=' 'AVG_Z'],'PgradVe',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Ubulk(ir,jr,kr).*divVe(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UbulkdivVe Z=' 'AVG_Z'],'UbulkdivVe',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Uth(ir,jr,kr).*divVe(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UthdivVe Z=' 'AVG_Z'],'UthdivVe',[-1 1]*0e-9, Nsm, 2);


end

ions=1
if(ions)
[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV] = compute_energy_balance( ...
    rhoi, Jix, Jiy, Jiz,... 
    Qbulkix, Qbulkiy, Qbulkiz, Qenthix, Qenthiy, Qenthiz, Qhfix, Qhfiy, Qhfiz, ...
    Pixx, Piyy, Pizz, Pixy, Pixz, Piyz, dx, dy, dz, 1.0)

Nsm=5;labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JidotE(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['JiE Z=' 'AVG_Z'],'JiE',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQbulk(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['divQbulki Z=' 'AVG_Z'],'divQbulki',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQenth(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['divQenthi Z=' 'AVG_Z'],'divQenthi',[-1 1]*0e-10, Nsm,1);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQhf(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['divQhfi Z=' 'AVG_Z'],'divhfi',[-1 1]*0e-10, Nsm,1);

labelc = 'mW/m^2'; 
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qbulkix(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkix Z=' 'AVG_Z'],'Qbulkix',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiy(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkiz Z=' 'AVG_Z'],'Qbulkiy',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiz(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkiy Z=' 'AVG_Z'],'Qbulkiz',[-1 1]*0e-9, Nsm, 4);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qenthix(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthix Z=' 'AVG_Z'],'Qenthix',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiy(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthiz Z=' 'AVG_Z'],'Qenthiy',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiz(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthiy Z=' 'AVG_Z'],'Qenthiz',[-1 1]*0e-9, Nsm, 4);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qhfix(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qhfix Z=' 'AVG_Z'],'Qhfix',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(hfiy(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qhfiz Z=' 'AVG_Z'],'Qhfiy',[-1 1]*0e-9, Nsm, 3);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(hfiz(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qhfiy Z=' 'AVG_Z'],'Qhfiz',[-1 1]*0e-9, Nsm, 4);


Qenthipar= dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthipar(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi || Z=' 'AVG_Z'],'Qenthipar',[-1 1]*0e-9, Nsm, 2);
Qenthiperp1=(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiperp1(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi \perp_1 Z=' 'AVG_Z'],'Qenthiprp1',[-1 1]*0e-9, Nsm, 2);
Qenthiperp2=perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiperp2(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi \perp_2 Z=' 'AVG_Z'],'Qenthiprp2',[-1 1]*0e-9, Nsm, 2);


Qbulkipar= dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkipar(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki || Z=' 'AVG_Z'],'Qbulkipar',[-1 1]*0e-9, Nsm, 2);
Qbulkiperp1=(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiperp1(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki \perp_1 Z=' 'AVG_Z'],'Qbulkiprp1',[-1 1]*0e-9, Nsm, 2);
Qbulkiperp2=perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiperp2(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki \perp_2 Z=' 'AVG_Z'],'Qbulkiprp2',[-1 1]*0e-9, Nsm, 2);

Nsm=5
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(udivP(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr) ,['UdivPi Z=' 'AVG_Z'],'UdivPi',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(PgradV(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr) ,['PgradVi Z=' 'AVG_Z'],'PgradVi',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Ubulk(ir,jr,kr).*divVi(ir,jr,kr),2)*nWm3,Vex(ir,kr),Viz(ir,kr) ,['UbulkdivVi Z=' 'AVG_Z'],'UbulkdivVi',[-1 1]*0e-9, Nsm, 2);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Uth(ir,jr,kr).*divVi(ir,jr,kr),2)*nWm3,Vex(ir,kr),Viz(ir,kr) ,['UthdivVi Z=' 'AVG_Z'],'UthdivVi',[-1 1]*0e-9, Nsm, 2);

end


%end


