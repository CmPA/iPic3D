%
% Energy plots in the XZ plane averaging on the whole range in y
%

close all
addpath(genpath('~/iPic3D-github/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are


UHRmaha3D3

%dir='/data1/gianni/HRmaha3D3/h5/'; cycle= 80002; ncycle = num2str(cycle,'%06d');

cycle = 80002  % for h5 in HRmaha3D3
cycle = 50004  % for h5 in UHRmaha3D3
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


cyl = 0 % Cartesian Geometry

poynting=true
electrons=true
ions=true
saveVTK=true


[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

bufferX=round(Nx/20);
bufferY=round(Ny/20);
%bufferZ=round(Nz/20);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
%kr=bufferZ:Nz-bufferZ;

ygsm=1.8
zcode = (ygsm - Ygsmrange(1)) * Lz/(Ygsmrange(2)-Ygsmrange(1))
iz= round(zcode/dz)%round(Nz/2)%135
if (Nz<11)
    kr=0;
else
    kr=-5:5
end    
kr=kr+round(iz);

Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the density
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3


radius=4;



global color_choice symmetric_color labelx labely labelc reversex reversey Ncycle skip
reversex=1;
reversey=1;
symmetric_color=1;
color_choice =3;
labelx ='x/R_E';
labely ='y/R_E';
labelc = 'mW/m^2';
skip=10


% Call the heavy lifting
energy_workhorse

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


%
% Elm energy
%



if(poynting)

labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JdotE(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'JE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divS(ir,jr,kr),3),Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'divS',[-1 1]*0e-9, radius, 4);


% The poynting flux is in SI units, W/m^3 so we need multiplication by 1e3
% to have it in mW/m^2
labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Sx(ir,jr,kr),3)*1e3,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sx',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sy(ir,jr,kr),3)*1e3,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Sz(ir,jr,kr),3)*1e3,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Sz',[-1 1]*0e-9, radius, 4);


if(saveVTK)
    savevtkvector_bin(Sx, Sy, Sz, [dir 'S' ncycle '.vtk'],'S',dx,dy,dz,0,0,0);
    savevtkvector_bin(Bx*code_B, By*code_B, Bz*code_B, [dir 'B' ncycle '.vtk'],'B',dx,dy,dz,0,0,0);
    savevtkvector_bin(Ex*code_E, Ey*code_E, Ez*code_E, [dir 'E' ncycle '.vtk'],'E',dx,dy,dz,0,0,0);
end

end

%
% Electrons
%

if(electrons)


labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JedotE(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'JeE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQbulk(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQbulke',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQenth(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQenthe',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQhf(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQhfe',[-1 1]*0e-10, radius,1);

labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qbulkex(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkey(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkez(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qenthex(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthey(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthez(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthez',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qhfex(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfex',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfey(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfey',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfez(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfez',[-1 1]*0e-9, radius, 4);


labelc = 'nW/m^3';

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivP(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(PgradV(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'PgradVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Ubulk(ir,jr,kr).*divVe(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'UbulkdivVe',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Uth(ir,jr,kr).*divVe(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'UthdivVe',[-1 1]*0e-9, radius, 2);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(DUbulkDt(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'DUbulkeDt',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(DUthDt(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'DUtheDt',[-1 1]*0e-9, radius, 2);

if(saveVTK)
    savevtkvector_bin(Qhfex*mWm2, Qhfey*mWm2, Qhfez*mWm2, [dir 'Qhfe' ncycle '.vtk'],'Qhfe',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qenthex*mWm2, Qenthey*mWm2, Qenthez*mWm2, [dir 'Qenthe' ncycle '.vtk'],'Qenthe',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qbulkex*mWm2, Qbulkey*mWm2, Qbulkez*mWm2, [dir 'Qbulke' ncycle '.vtk'],'Qbulke',dx,dy,dz,0,0,0);
    savevtk_bin(DUbulkDt*nWm3,[dir 'DUbulkedt' ncycle '.vtk'],'DUbulkedt',dx,dy,dz,0,0,0);
    savevtk_bin(DUthDt*nWm3,[dir 'DUthedt' ncycle '.vtk'],'DUthedt',dx,dy,dz,0,0,0);
end

end


if(ions)
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(JidotE(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'JiE',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQbulk(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQbulki',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQenth(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr),['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQenthi',[-1 1]*0e-10, radius,1);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(divQhf(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'divQhfi',[-1 1]*0e-10, radius,1);

labelc = 'mW/m^2'; 
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qbulkix(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiy(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qbulkiz(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qbulkiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qenthix(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiy(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qenthiz(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qenthiz',[-1 1]*0e-9, radius, 4);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),-mean(Qhfix(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfix',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfiy(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfiy',[-1 1]*0e-9, radius, 3);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Qhfiz(ir,jr,kr),3)*mWm2,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'Qhfiz',[-1 1]*0e-9, radius, 4);


labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(udivP(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr) ,['Y=' num2str(gsmz2y(z(1,1,iz)))],'UdivPi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(PgradV(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'PgradVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Ubulk(ir,jr,kr).*divVi(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'UbulkdivVi',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(Uth(ir,jr,kr).*divVi(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr) , ['Y=' num2str(gsmz2y(z(1,1,iz)))],'UthdivVi',[-1 1]*0e-9, radius, 2);

tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(DUbulkDt(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'DUbulkiDt',[-1 1]*0e-9, radius, 2);
tmp=common_image_vel(gsmx(X(jr,ir)),gsmy2z(Y(jr,ir)),mean(DUthDt(ir,jr,kr),3)*nWm3,Vx(ir,jr),Vy(ir,jr), ['Y=' num2str(gsmz2y(z(1,1,iz)))],'DUthiDt',[-1 1]*0e-9, radius, 2);

if(saveVTK)
    savevtkvector_bin(Qhfix*mWm2, Qhfiy*mWm2, Qhfiz*mWm2, [dir 'Qhfi' ncycle '.vtk'],'Qhfi',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qenthix*mWm2, Qenthiy*mWm2, Qenthiz*mWm2, [dir 'Qenthi' ncycle '.vtk'],'Qenthi',dx,dy,dz,0,0,0);
    savevtkvector_bin(Qbulkix*mWm2, Qbulkiy*mWm2, Qbulkiz*mWm2, [dir 'Qbulki' ncycle '.vtk'],'Qbulki',dx,dy,dz,0,0,0);
    savevtk_bin(DUbulkDt*nWm3,[dir 'DUbulkidt' ncycle '.vtk'],'DUbulkidt',dx,dy,dz,0,0,0);
    savevtk_bin(DUthDt*nWm3,[dir 'DUthidt' ncycle '.vtk'],'DUthidt',dx,dy,dz,0,0,0);
end

end




