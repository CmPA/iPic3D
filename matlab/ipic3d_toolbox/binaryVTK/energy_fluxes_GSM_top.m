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
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
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
radius=5;
JedotEsm=dot(smooth3(Jex,method,radius),smooth3(Jey,method,radius),smooth3(Jez,method,radius), ...
    smooth3(Ex,method,radius),smooth3(Ey,method,radius),smooth3(Ez,method,radius));


JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JdotE=JedotE+JidotE;

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
Sx=Sx*code_E*code_B/mu0;
SZ=Sy*code_E*code_B/mu0;
Sz=Sz*code_E*code_B/mu0;

xc=Lx-linspace(0, Lx, Nx);
zc=linspace(0, Lz, Nz);
Wm3 = code_E*code_J*4*pi; %4pi is due to the usal division by 4pi of the dencity
nWm3 = 1e9*Wm3;
mWm2= Wm3*code_dp*1e3


for iz=135
%kr=-5:5
%kr=kr+round(iz);
Nsm=5


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

poynting=0
if(poynting)

labelc = 'nW/m^3';
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JdotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JE Z=' 'AVG_Z'],'JE',[-1 1]*0e-10, Nsm,1+iz);

labelc = 'mW/m^2';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Sx(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sx Z=' 'AVG_Z'],'Sx',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sy(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sz Z=' 'AVG_Z'],'Sy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sz(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['Sy Z=' 'AVG_Z'],'Sz',[-1 1]*0e-9, Nsm, 4+iz);

%Spar= dot(Sx,Sy,Sz,Bx,By,Bz)./B;
%tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Spar(ir,jr,kr),2)*1e3,AAz(ir,kr) ,['S_{||} Z=' 'AVG_Z'],'Spar',[-1 1]*0e-9, Nsm, 2+iz);
Sperp1=(By.*Sx-Bx.*Sy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp1(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['S \perp_1 Z=' 'AVG_Z'],'Sperp1',[-1 1]*0e-9, Nsm, 2+iz);
Sperp2=perp2x.*Sx+perp2y.*Sy+perp2z.*Sz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Sperp2(ir,jr,kr),2)*1e3,Vex(ir,kr),Vez(ir,kr) ,['S \perp_2 Z=' 'AVG_Z'],'Sperp2',[-1 1]*0e-9, Nsm, 2+iz);

end

electrons=1
if(electrons)

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JedotE(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr),['JeE Z=' 'AVG_Z'],'JeE',[-1 1]*0e-10, Nsm,1+iz);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qbulkex(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkex Z=' 'AVG_Z'],'Qbulkex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkey(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkez Z=' 'AVG_Z'],'Qbulkey',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkez(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulkey Z=' 'AVG_Z'],'Qbulkez',[-1 1]*0e-9, Nsm, 4+iz);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qenthex(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthex Z=' 'AVG_Z'],'Qenthex',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthey(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthez Z=' 'AVG_Z'],'Qenthey',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthez(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthey Z=' 'AVG_Z'],'Qenthez',[-1 1]*0e-9, Nsm, 4+iz);


Qenthepar= dot(Qenthex,Qenthey,Qenthez,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthepar(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenthe || Z=' 'AVG_Z'],'Qenthepar',[-1 1]*0e-9, Nsm, 2+iz);
Qentheperp1=(By.*Qenthex-Bx.*Qenthey)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qentheperp1(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenth \perp_1 Z=' 'AVG_Z'],'Qentheprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qentheperp2=perp2x.*Qenthex+perp2y.*Qenthey+perp2z.*Qenthez;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qentheperp2(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qenth \perp_2 Z=' 'AVG_Z'],'Qentheprp2',[-1 1]*0e-9, Nsm, 2+iz);


Qbulkepar= dot(Qbulkex,Qbulkey,Qbulkez,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkepar(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulke || Z=' 'AVG_Z'],'Qbulkepar',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkeperp1=(By.*Qbulkex-Bx.*Qbulkey)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkeperp1(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulk \perp_1 Z=' 'AVG_Z'],'Qbulkeprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkeperp2=perp2x.*Qbulkex+perp2y.*Qbulkey+perp2z.*Qbulkez;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkeperp2(ir,jr,kr),2)*mWm2,Vex(ir,kr),Vez(ir,kr) ,['Qbulk \perp_2 Z=' 'AVG_Z'],'Qbulkeprp2',[-1 1]*0e-9, Nsm, 2+iz);

Nsm=10
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(UdivPe(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UdivPe Z=' 'AVG_Z'],'UdivPe',[-1 1]*0e-9, Nsm, 2+iz);

newsmooth=0
if (newsmooth)
radius=5
method='gaussian';
Vx=smooth3(Jex./rhoe,method,radius);
Vy=smooth3(Jey./rhoe,method,radius);
Vz=smooth3(Jez./rhoe,method,radius);
tmp = divergence(x,y,z,smooth3(permute(Pexx,[2 1 3]),method,radius), smooth3(permute(Pexy, [2 1 3]),method,radius), smooth3(permute(Pexz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = tmp.* Vx;
tmp = divergence(x,y,z,smooth3(permute(Pexy,[2 1 3]),method,radius), smooth3(permute(Peyy, [2 1 3]),method,radius), smooth3(permute(Peyz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vy;
tmp = divergence(x,y,z,smooth3(permute(Pexz,[2 1 3]),method,radius), smooth3(permute(Peyz, [2 1 3]),method,radius), smooth3(permute(Pezz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vz;

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(udivP(ir,jr,kr),2)*nWm3,Vex(ir,kr),Vez(ir,kr) ,['UdivPe Z=' 'AVG_Z'],'UdivPe2',[-1 1]*0e-9, Nsm, 2+iz);
end

end

ions=0
if(ions)
labelc = 'mW/m^2'; Nsm=5;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(JidotE(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['JiE Z=' 'AVG_Z'],'JiE',[-1 1]*0e-10, Nsm,1+iz);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qbulkix(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkix Z=' 'AVG_Z'],'Qbulkix',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiy(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkiz Z=' 'AVG_Z'],'Qbulkiy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiz(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulkiy Z=' 'AVG_Z'],'Qbulkiz',[-1 1]*0e-9, Nsm, 4+iz);

tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),-mean(Qenthix(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthix Z=' 'AVG_Z'],'Qenthix',[-1 1]*0e-9, Nsm, 2+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiy(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthiz Z=' 'AVG_Z'],'Qenthiy',[-1 1]*0e-9, Nsm, 3+iz);
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiz(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthiy Z=' 'AVG_Z'],'Qenthiz',[-1 1]*0e-9, Nsm, 4+iz);


Qenthipar= dot(Qenthix,Qenthiy,Qenthiz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthipar(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi || Z=' 'AVG_Z'],'Qenthipar',[-1 1]*0e-9, Nsm, 2+iz);
Qenthiperp1=(By.*Qenthix-Bx.*Qenthiy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiperp1(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi \perp_1 Z=' 'AVG_Z'],'Qenthiprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qenthiperp2=perp2x.*Qenthix+perp2y.*Qenthiy+perp2z.*Qenthiz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qenthiperp2(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qenthi \perp_2 Z=' 'AVG_Z'],'Qenthiprp2',[-1 1]*0e-9, Nsm, 2+iz);


Qbulkipar= dot(Qbulkix,Qbulkiy,Qbulkiz,Bx,By,Bz)./B;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkipar(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki || Z=' 'AVG_Z'],'Qbulkipar',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkiperp1=(By.*Qbulkix-Bx.*Qbulkiy)./B2D;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiperp1(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki \perp_1 Z=' 'AVG_Z'],'Qbulkiprp1',[-1 1]*0e-9, Nsm, 2+iz);
Qbulkiperp2=perp2x.*Qbulkix+perp2y.*Qbulkiy+perp2z.*Qbulkiz;
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(Qbulkiperp2(ir,jr,kr),2)*mWm2,Vix(ir,kr),Viz(ir,kr) ,['Qbulki \perp_2 Z=' 'AVG_Z'],'Qbulkiprp2',[-1 1]*0e-9, Nsm, 2+iz);

Nsm=5
labelc = 'nW/m^3';
tmp=common_image_vel(gsmx(X(kr,ir)),gsmz2y(Z(kr,ir)),mean(UdivPi(ir,jr,kr),2)*nWm3,Vix(ir,kr),Viz(ir,kr),['UdivPi Z=' 'AVG_Z'],'UdivPi',[-1 1]*0e-9, Nsm, 2+iz);

end


end

end
