close all
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are

HRmaha3D3

for cycle=80000:1000:80000
%for cycle=118000:1000:118000

% for HRmaha3D1:
 time=60*(cycle/75000.0*Dt/.125) %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D

%ADD initial time of the RUN
time=time+initial_time%(03*60+48)*60
% Prepare string
ntime = datestr(time/86400,'HH:MM:SS UT')



    ncycle = num2str(cycle,'%06d');
leggo=0; 
if(leggo==1)


[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

% 
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
B2D=sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(B.*B2D);
perp2y=Bz.*By./(B.*B2D);
perp2z=-B2D./B;
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
 [Qbulkex,Qbulkey,Qbulkez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qbulke',cycle);
 [Qenthex,Qenthey,Qenthez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qenthe',cycle);
 [Qbulkix,Qbulkiy,Qbulkiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qbulki',cycle);
 [Qenthix,Qenthiy,Qenthiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Qenthi',cycle);
[UdivPe,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPe',cycle);
[UdivPi,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'UdivPi',cycle);
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;
end


Lx=dx*Nx;LZ=dy*Ny;Lz=Nz*dz;

[x,y,z]=meshgrid(0:dx:Lx-dx,0:dy:Ly-dy,0:dz:Lz-dz);

[X Z] = meshgrid(0:dx:Lx-dx,0:dz:Lz-dz);

qom_ele = -256;

bufferX=round(Nx/20);
bufferY=round(Ny/20);
bufferZ=round(Nz/20);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
kr=bufferZ:Nz-bufferZ;


jr=round(Ny/2)+[-5:5];

Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the density
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3
%
% Electrons
%


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
radius=10;
JedotEsm=dot(smooth3Dlike(Jex,method,radius),smooth3Dlike(Jey,method,radius),smooth3Dlike(Jez,method,radius), ...
    smooth3Dlike(Ex,method,radius),smooth3Dlike(Ey,method,radius),smooth3Dlike(Ez,method,radius));
JidotEsm=dot(smooth3Dlike(Jix,method,radius),smooth3Dlike(Jiy,method,radius),smooth3Dlike(Jiz,method,radius), ...
    smooth3Dlike(Ex,method,radius),smooth3Dlike(Ey,method,radius),smooth3Dlike(Ez,method,radius));


JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JdotE=JedotEsm+JidotEsm;

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);

divS = compute_div(x,y,z,smooth3Dlike(Sx,method,radius),smooth3Dlike(Sy,method,radius),smooth3Dlike(Sz,method,radius));

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


for iz=135
%kr=-5:5
%kr=kr+round(iz);
Nsm=10


% Vix=Jix./rhoi;Viz=Jiz./rhoi;
% AAzi=vecpot(xc,zc,-squeeze(mean(Vix(:,jr,:),2)),squeeze(mean(Viz(:,jr,:),2)));

Vex=Jex./rhoe;Vez=Jez./rhoe;
Vex=-smoothbc(squeeze(mean(Vex(:,jr,:),2)),Nsm);
Vez=smoothbc(squeeze(mean(Vez(:,jr,:),2)),Nsm);
AAze=vecpot(xc,zc,Vex,Vez);

Vix=Jix./rhoi;Viz=Jiz./rhoi;
Vix=-smoothbc(squeeze(mean(Vix(:,jr,:),2)),Nsm);
Viz=smoothbc(squeeze(mean(Viz(:,jr,:),2)),Nsm);
AAze=vecpot(xc,zc,Vix,Viz);

poynting=1
if(poynting)

labelc = 'nW/m^3';
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*1.5e-2, Nsm,1+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divS(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr) ,['divS Z=' 'AVG_Z'],'divS',[-1 1]*0, Nsm, 2+iz);

end

electrons=1
if(electrons)
    
Tepar = 3 * code_T * Pepar./(-rhoe) /e/1e3;   
Teperp = 3 * code_T * .5*(Peper1+Peper2)./(-rhoe)/e/1e3; 
Te = code_T *(Pexx+Peyy+Pezz)./(-rhoe) /e/1e3;
Ubulke=.5*(Jex.*Jex+Jey.*Jey+Jez.*Jez)./(-rhoe).^2 *code_T/e/1e3/(-qom_ele);
 
divQbulke = compute_div(x,y,z,smooth3Dlike(Qbulkex,method,radius),smooth3Dlike(Qbulkey,method,radius),smooth3Dlike(Qbulkez,method,radius));    
divQenthe = compute_div(x,y,z,smooth3Dlike(Qenthex,method,radius),smooth3Dlike(Qenthey,method,radius),smooth3Dlike(Qenthez,method,radius));


labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JedotEsm(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JeE Z=' 'AVG_Z'],'JeE',[-1 1]*3e-2, Nsm,3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQbulke(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['divQbulke Z=' 'AVG_Z'],'divQbulke',[-1 1]*3e-2, Nsm,4+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQenthe(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['divQenthe Z=' 'AVG_Z'],'divQenthe',[-1 1]*3e-2, Nsm,5+iz);

symmetric_color=0;
color_choice =0;
labelc = 'keV';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Tepar(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Tepar Z=' 'AVG_Z'],'Tepar',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Teperp(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Teperp Z=' 'AVG_Z'],'Teperp',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Te(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Te Z=' 'AVG_Z'],'Te',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Ubulke(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Ubulke Z=' 'AVG_Z'],'Ubulke',[-1 1]*0e-10, Nsm,5+iz);


symmetric_color=1;
color_choice =3;

labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(UdivPe(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UdivPe Z=' 'AVG_Z'],'UdivPe',[-1 1]*3e-2, Nsm, 6+iz);

newsmooth=1
if (newsmooth)
%radius=5
%method='gaussian';
Vx=smooth3Dlike(Jex./rhoe,method,radius);
Vy=smooth3Dlike(Jey./rhoe,method,radius);
Vz=smooth3Dlike(Jez./rhoe,method,radius);
tmp = divergence(x,y,z,smooth3Dlike(permute(Pexx,[2 1 3]),method,radius), smooth3Dlike(permute(Pexy, [2 1 3]),method,radius), smooth3Dlike(permute(Pexz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = tmp.* Vx;
tmp = divergence(x,y,z,smooth3Dlike(permute(Pexy,[2 1 3]),method,radius), smooth3Dlike(permute(Peyy, [2 1 3]),method,radius), smooth3Dlike(permute(Peyz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vy;
tmp = divergence(x,y,z,smooth3Dlike(permute(Pexz,[2 1 3]),method,radius), smooth3Dlike(permute(Peyz, [2 1 3]),method,radius), smooth3Dlike(permute(Pezz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vz;

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(udivP(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UdivPe Z=' 'AVG_Z'],'UdivPe2',[-1 1]*3e-2, Nsm, 7+iz);
end
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),(-mean(divQenthe(ir,jr,kr),2)+mean(UdivPe(ir,jr,kr),2))*nWm3,Vex(ir,kr),Vez(ir,kr) ,['Source Z=' 'AVG_Z'],'Source',[-1 1]*3e-2, Nsm, 7+iz);
end

end

ions=1
if(ions)
divQbulki = compute_div(x,y,z,smooth3Dlike(Qbulkix,method,radius),smooth3Dlike(Qbulkiy,method,radius),smooth3Dlike(Qbulkiz,method,radius));

divQenthi = compute_div(x,y,z,smooth3Dlike(Qenthix,method,radius),smooth3Dlike(Qenthiy,method,radius),smooth3Dlike(Qenthiz,method,radius));


Tipar = 3 * code_T * Pipar./(rhoi) /e/1e3;   
Tiperp = 3 *code_T * .5*(Piper1+Piper2)./(rhoi)/e/1e3; 
Ti = code_T *(Pixx+Piyy+Pizz)./(rhoi) /e/1e3;
Ubulki=.5*(Jix.*Jix+Jiy.*Jiy+Jiz.*Jiz)./(rhoi).^2 *code_T/e/1e3;
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JidotEsm(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['JiE Z=' 'AVG_Z'],'JiE',[-1 1]*0e-10, Nsm,8+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQbulki(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['divQbulki Z=' 'AVG_Z'],'divQbulki',[-1 1]*0e-10, Nsm,9+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(divQenthi(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['divQenthi Z=' 'AVG_Z'],'divQenthi',[-1 1]*0e-10, Nsm,10+iz);

symmetric_color=0;
color_choice =0;
labelc = 'keV';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Tipar(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Tipar Z=' 'AVG_Z'],'Tipar',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Tiperp(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Tiperp Z=' 'AVG_Z'],'Tiperp',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Ti(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Ti Z=' 'AVG_Z'],'Ti',[-1 1]*0e-10, Nsm,5+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Ubulki(ir,jr,kr),2),Vex(ir,kr),Vez(ir,kr),['Ubulki Z=' 'AVG_Z'],'Ubulki',[-1 1]*0e-10, Nsm,5+iz);

symmetric_color=1;
color_choice =3;

labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(UdivPi(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['UdivPi Z=' 'AVG_Z'],'UdivPi',[-1 1]*0e-9, Nsm, 11+iz);

end


end

end
