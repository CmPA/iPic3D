close all
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
dir='/shared/gianni/tred70/'; %directory where the files are

for cycle=21000:1000:21000

    ncycle = num2str(cycle,'%06d');
leggo=0; part1 =1; part2 =1;
if(leggo==1)


[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
%[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

% 
[Az,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Az',cycle);
[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
% B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
% Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
 
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;
end

Lx=dx*Nx;Ly=dy*Ny;Lz=Nz*dz;

[x,y,z]=meshgrid(0:dx:Lx-dx,0:dy:Ly-dy,0:dz:Lz-dz);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

qom_ele = -256;

bufferX=round(Nx/6);
bufferY=round(Ny/6);
ir=bufferX:Nx-bufferX;
jr=bufferY:Ny-bufferY;
kr=3:Nz-3;
Nsm=0

%
% Electrons
%


if(part1)
% Compute J dot E

JdotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);
savevtk_bin(JdotE,[dir 'JedotT' ncycle '.vtk'],'JedotE',dx,dy,dz,0,0,0);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(JdotE(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'JeE','JeE',[0 0], Nsm, 22);

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
savevtkvector_bin(Sx, Sy, Sz, [dir 'Poynting' ncycle '.vtk'],'Poynting',dx,dy,dz,0,0,0);

clear Sx Sy Sz
% compute bulk energy

Ubulk = 0.5 .* (Jex.^2 + Jey.^2 + Jez.^2) ./ rhoe /qom_ele;


savevtk_bin(Ubulk,[dir 'Ubulke' ncycle '.vtk'],'Ubulke',dx,dy,dz,0,0,0);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(Ubulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Ubulke','Ubulke',[0 0], Nsm, 1);

% Compute thermal energy

Uth = 0.5 .* (Pexx + Peyy + Pezz);
savevtk_bin(Uth,[dir 'Uthe' ncycle '.vtk'],'Uthe',dx,dy,dz,0,0,0);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(Uth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Uthe','Uthe',[0 0], Nsm, 2);

% Compute bulk energy fluxes

Qxbulk = Ubulk.*Jex ./rhoe;
Qybulk = Ubulk.*Jey ./rhoe;
Qzbulk = Ubulk.*Jez ./rhoe;

savevtkvector_bin(Qxbulk, Qybulk, Qzbulk, [dir 'Qbulke' ncycle '.vtk'],'Qbulke',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qxbulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qxbulke','Qxbulke',[0 0], Nsm, 3);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qybulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qybulke','Qybulke',[0 0], Nsm, 4);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qzbulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qzbulke','Qzbulke',[0 0], Nsm, 5);

% Compute thermal enrgy fluxes and enthalpy fluxes

Qxth = dot( Jex, Jey, Jez, Pexx, Pexy, Pexz, rhoe);
Qxenth = Qxth + Uth.*Jex ./rhoe;
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qxth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qxthe','Qxteh',[0 0], Nsm, 6);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qxenth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qxenthe','Qxenthe',[0 0], Nsm, 7);

Qyth = dot( Jex, Jey, Jez, Pexy, Peyy, Peyz, rhoe);
Qyenth = Qyth + Uth.*Jey ./rhoe;
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qyth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qythe','Qythe',[0 0], Nsm, 8);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qyenth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qyenthe','Qyenthe',[0 0], Nsm, 9);

Qzth = dot( Jex, Jey, Jez, Pexz, Peyz, Pezz, rhoe);
Qzenth = Qzth + Uth.*Jez ./rhoe;

savevtkvector_bin(Qxth, Qyth, Qzth, [dir 'Qthe' ncycle '.vtk'],'Qthe',dx,dy,dz,0,0,0);
savevtkvector_bin(Qxenth, Qyenth, Qzenth, [dir 'Qenthe' ncycle '.vtk'],'Qenthe',dx,dy,dz,0,0,0);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qzth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qzthe','Qzthe',[0 0], Nsm, 10);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Qzenth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Qzenthe','Qzenthe',[0 0], Nsm, 11);


tmp = divergence(x,y,z,permute(Qxth,[2 1 3]), permute(Qyth, [2 1 3]), permute(Qzth, [2,1,3]));
tmp = permute(tmp, [2 1 3]);
EULth = -tmp;

savevtk_bin(tmp,[dir 'divQthe' ncycle '.vtk'],'divQthe',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(tmp(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'divQthe','divQthe',[0 0], Nsm, 12);



clear Qxth Qyth Qzth Qxenth Qyenth Qzenth 


tmp = divergence(x,y,z,permute(Qxbulk,[2 1 3]), permute(Qybulk, [2 1 3]), permute(Qzbulk, [2,1,3]));
tmp = permute(tmp, [2 1 3]);
EULbulk = JdotE - tmp; 

savevtk_bin(tmp,[dir 'divQbulke' ncycle '.vtk'],'divQbulke',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(tmp(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'divQbulke','divQbulke',[0 0], Nsm, 13);



clear Qxbulk Qybulk Qzbulk

Vx=Jex./rhoe;
Vy=Jey./rhoe;
Vz=Jez./rhoe;
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Vx(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Vex','Vex',[0 0], Nsm, 14);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Vy(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Vey','Vey',[0 0], Nsm, 15);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Vz(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Vez','Vez',[0 0], Nsm, 16);

% Computing div u
divu = divergence(x,y,z,permute(Vx,[2 1 3]), permute(Vy, [2 1 3]), permute(Vz, [2,1,3]));
divu=permute(divu,[2 1 3]);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(divu(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'divUe','divUe',[0 0], Nsm, 17);
Udivu=divu.*Ubulk;

LAGbulk = - Udivu + JdotE ;


tmp=common_image(X(jr,ir),Y(jr,ir),mean(Udivu(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Ubulkdivue','Ubulkdivue',[0 0], Nsm, 18);
Udivu=divu.*Uth;
tmp=common_image(X(jr,ir),Y(jr,ir),mean(Udivu(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'Uthdivue','Uthdivue',[0 0], Nsm, 19);

LAGth = - Udivu;

clear Uth Ubulk Udivu divu JdotE
end


if(part2)
% Computing u . div(P)

tmp = divergence(x,y,z,permute(Pexx,[2 1 3]), permute(Pexy, [2 1 3]), permute(Pexz, [2,1,3]));
tmp=permute(tmp,[2 1 3]);
udivP = tmp.* Vx;
tmp = divergence(x,y,z,permute(Pexy,[2 1 3]), permute(Peyy, [2 1 3]), permute(Peyz, [2,1,3]));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vy;
tmp = divergence(x,y,z,permute(Pexz,[2 1 3]), permute(Peyz, [2 1 3]), permute(Pezz, [2,1,3]));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vz;

savevtk_bin(udivP,[dir 'udivPe' ncycle '.vtk'],'udivPe',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(udivP(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'UdivPe','UdivPe',[0 0], Nsm, 20);

LAGbulk = LAGbulk - udivP ;
EULbulk = EULbulk - udivP ;
EULth = EULth + udivP;

clear udivP

% Computing P : grad(U)

[tmpx tmpy tmpz] = gradient(permute(Vx,[2 1 3]),dx,dy,dz);
tmp=dot(permute(tmpx,[2 1 3]),permute(tmpy,[2 1 3]),permute(tmpz,[2 1 3]),Pexx, Pexy, Pexz);
[tmpx tmpy tmpz] = gradient(permute(Vy, [2 1 3]),dx,dy,dz);
tmp=tmp+dot(permute(tmpx,[2 1 3]),permute(tmpy,[2 1 3]),permute(tmpz,[2 1 3]),Pexy, Peyy, Peyz);
[tmpx tmpy tmpz] = gradient(permute(Vz, [2 1 3]),dx,dy,dz);
tmp=tmp+dot(permute(tmpx,[2 1 3]),permute(tmpy,[2 1 3]),permute(tmpz,[2 1 3]),Pexz, Peyz, Pezz);
LAGth = LAGth - tmp;

savevtk_bin(tmp,[dir 'Pgradue' ncycle '.vtk'],'Pgradue',dx,dy,dz,0,0,0);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(tmp(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'PgradUe','PgradUe',[0 0], Nsm, 21);


end
savevtk_bin(EULth,[dir 'EULthe' ncycle '.vtk'],'EULthe',dx,dy,dz,0,0,0);
savevtk_bin(EULbulk,[dir 'EULbulke' ncycle '.vtk'],'EULbulke',dx,dy,dz,0,0,0);
savevtk_bin(LAGth,[dir 'LAGthe' ncycle '.vtk'],'LAGthe',dx,dy,dz,0,0,0);
savevtk_bin(LAGbulk,[dir 'LAGbulke' ncycle '.vtk'],'LAGbulke',dx,dy,dz,0,0,0);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(EULth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'EULthe','EULthe',[-4 4]*1e-8, Nsm, 23);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(EULbulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'EULbulke','EULbulke',[-4 4]*1e-8, Nsm, 24);

tmp=common_image(X(jr,ir),Y(jr,ir),mean(LAGth(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'LAGthe','LAGthe',[-4 4]*1e-8, Nsm, 25);
tmp=common_image(X(jr,ir),Y(jr,ir),mean(LAGbulk(ir,jr,kr),3), mean(Az(ir,jr,kr),3),'LAGbulke','LAGbulke',[-4 4]*1e-8, Nsm, 26);

clear EULth EULbulk LAGth LAGbulk 



end