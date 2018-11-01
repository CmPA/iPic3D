close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are



sim_name='tred81'
switch sim_name
case 'tred81'
tred81;
case_name='GEM';
cycle = 22000;
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




    ncycle = num2str(cycle,'%06d');
leggo=true
if(leggo)
    
    
    namefile = [case_name '-Fields'];
    fn=[dir,namefile,'_',ncycle,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3);
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    

    Bx = hdf5read(fn,'/Step#0/Block/Bx/0/');
    By = hdf5read(fn,'/Step#0/Block/By/0/');
    Bz = hdf5read(fn,'/Step#0/Block/Bz/0/');
    Ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    Ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    Ez = hdf5read(fn,'/Step#0/Block/Ez/0/');
    Jex = hdf5read(fn,'/Step#0/Block/Jx_0/0/')+hdf5read(fn,'/Step#0/Block/Jx_2/0/');
    Jey = hdf5read(fn,'/Step#0/Block/Jy_0/0/')+hdf5read(fn,'/Step#0/Block/Jy_2/0/');
    Jez = hdf5read(fn,'/Step#0/Block/Jz_0/0/')+hdf5read(fn,'/Step#0/Block/Jz_2/0/');
    Jix = hdf5read(fn,'/Step#0/Block/Jx_1/0/')+hdf5read(fn,'/Step#0/Block/Jx_3/0/');
    Jiy = hdf5read(fn,'/Step#0/Block/Jy_1/0/')+hdf5read(fn,'/Step#0/Block/Jy_3/0/');
    Jiz = hdf5read(fn,'/Step#0/Block/Jz_1/0/')+hdf5read(fn,'/Step#0/Block/Jz_3/0/');
    
    rhoe = hdf5read(fn,'/Step#0/Block/rho_0/0/')+hdf5read(fn,'/Step#0/Block/rho_2/0/');
    rhoi = hdf5read(fn,'/Step#0/Block/rho_1/0/')+hdf5read(fn,'/Step#0/Block/rho_3/0/');

    Pexx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/')+hdf5read(fn,'/Step#0/Block/Pxx_2/0/');
    Peyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/')+hdf5read(fn,'/Step#0/Block/Pyy_2/0/');
    Pezz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/')+hdf5read(fn,'/Step#0/Block/Pzz_2/0/');
    Pexy = hdf5read(fn,'/Step#0/Block/Pxy_0/0/')+hdf5read(fn,'/Step#0/Block/Pxy_2/0/');    
    Pexz = hdf5read(fn,'/Step#0/Block/Pxz_0/0/')+hdf5read(fn,'/Step#0/Block/Pxz_2/0/');
    Peyz = hdf5read(fn,'/Step#0/Block/Pyz_0/0/')+hdf5read(fn,'/Step#0/Block/Pyz_2/0/');
    
    Pixx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/')+hdf5read(fn,'/Step#0/Block/Pxx_3/0/');
    Piyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/')+hdf5read(fn,'/Step#0/Block/Pyy_3/0/');
    Pizz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/')+hdf5read(fn,'/Step#0/Block/Pzz_3/0/');
    Pixy = hdf5read(fn,'/Step#0/Block/Pxy_1/0/')+hdf5read(fn,'/Step#0/Block/Pxy_3/0/');    
    Pixz = hdf5read(fn,'/Step#0/Block/Pxz_1/0/')+hdf5read(fn,'/Step#0/Block/Pxz_3/0/');
    Piyz = hdf5read(fn,'/Step#0/Block/Pyz_1/0/')+hdf5read(fn,'/Step#0/Block/Pyz_3/0/');
    B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
    B2D=sqrt(Bx.^2+By.^2);
    perp2x=Bz.*Bx./(B.*B2D);
    perp2y=Bz.*By./(B.*B2D);
    perp2z=-B2D./B;
    Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

    [Pexx,Peyy,Pezz,Pexy,Pexz,Peyz]=compute_pressure(Bx, By, Bz, Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez,rhoe, qom);
    [Pixx,Piyy,Pizz,Pixy,Pixz,Piyz]=compute_pressure(Bx, By, Bz, Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz,rhoi, 1.0);
    
    Qex = hdf5read(fn,'/Step#0/Block/EFx_0/0/')+hdf5read(fn,'/Step#0/Block/EFx_2/0/');
    Qey = hdf5read(fn,'/Step#0/Block/EFy_0/0/')+hdf5read(fn,'/Step#0/Block/EFy_2/0/');
    Qez = hdf5read(fn,'/Step#0/Block/EFz_0/0/')+hdf5read(fn,'/Step#0/Block/EFz_2/0/');    
    Qix = hdf5read(fn,'/Step#0/Block/EFx_1/0/')+hdf5read(fn,'/Step#0/Block/EFx_3/0/');
    Qiy = hdf5read(fn,'/Step#0/Block/EFy_1/0/')+hdf5read(fn,'/Step#0/Block/EFy_3/0/');
    Qiz = hdf5read(fn,'/Step#0/Block/EFz_1/0/')+hdf5read(fn,'/Step#0/Block/EFz_3/0/'); 
  
    [Qenthex,Qenthey,Qenthez,Qbulkex,Qbulkey,Qbulkez,Qhfex,Qhfey,Qhfez] = ...
    compute_energy_fluxes(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Qex,Qey,Qez,Jex,Jey,Jez,rhoe, qom);

    [Qenthix,Qenthiy,Qenthiz,Qbulkix,Qbulkiy,Qbulkiz,Qhfix,Qhfiy,Qhfiz] = ...
    compute_energy_fluxes(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Qix,Qiy,Qiz,Jix,Jiy,Jiz,rhoi, 1.0);

end

[nx ny nz]= size(Pexx)

for i=1:nx
for iy=1:ny
for k=1:nz

p(1,1)=(Pexx(i,iy,k));
p(1,2)=(Pexy(i,iy,k));
p(1,3)=(Pexz(i,iy,k));
p(2,2)=(Peyy(i,iy,k));
p(2,3)=(Peyz(i,iy,k));
p(3,3)=(Pezz(i,iy,k));
p(2,1)=p(1,2);
p(3,1)=p(1,3);
p(3,2)=p(2,3);

b(1)=(Bx(i,iy,k));
b(2)=(By(i,iy,k));
b(3)=(Bz(i,iy,k));

b=reshape(b,3,1);
b=b./sqrt(sum(b.^2));

%%%%%%%%%%
% Scudder
%%%%%%%%%%

for l=1:3
N1(l,:)=cross(b,p(l,:));
end

for l=1:3
N(:,l)=cross(b,N1(:,l));
end

lamb=sort(eig(N));
lambda1(i,iy,k)=lamb(1);
lambda2(i,iy,k)=lamb(2);
lambda3(i,iy,k)=lamb(3);

Agyro(i,iy,k)= 2*(lamb(3)-lamb(2))/(lamb(3)+lamb(2));


%%%%%%%%%%%
% Aunai
%%%%%%%%%%%

Tr=(p(1,1)+p(2,2)+p(3,3));
Ppar=b'*p*b;
Pper=(Tr-Ppar)/2;
G=eye(3,3)*Pper+(Ppar-Pper)*kron(b,b');
N=p-G;
Agyro_aunai(i,iy,k)=sqrt(sum(N(:).^2))./Tr;

%%%%%%%%%%%
% Swisdak
%%%%%%%%%%%

I2=p(1,1)*p(2,2)+p(1,1)*p(3,3)+p(2,2)*p(3,3);
I2=I2-(p(1,2).^2+p(1,3).^2+p(2,3).^2);
Q=1-4*I2./((Tr-Ppar).*(Tr+3*Ppar));
Nongyro_swisdak(i,iy,k)=sqrt(Q);
% The following formula form Swisdak paper is actually wrong
%Nongyro_aunai(i,k)=sqrt(8*(p(1,2).^2+p(1,3).^2+p(2,3).^2))./(Ppar+2*Pper);



end
end
end
Agyro_sm=smooth3D(Agyro,6);
Agyro_aunai_sm=smooth3D(Agyro_aunai,6);
Nongyro_swisdak_sm=smooth3D(Nongyro_swisdak,6);

savevtk_bin(Agyro_sm,['Agyro_xyz_cycle' ncycle '.vtk'],'Agyro',dx,dy,dz)
savevtk_bin(Agyro_aunai_sm,['Agyro_aunai_xyz_cycle' ncycle '.vtk'],'Agyro',dx,dy,dz)
savevtk_bin(Nongyro_swisdak_sm,['Nongyro_swisdak_xyz_cycle' ncycle '.vtk'],'Agyro',dx,dy,dz)

opath=fn
h5create(opath,'/Step#0/Block/Agyro/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Agyro/0',Agyro)
h5create(opath,'/Step#0/Block/Agyro_aunai/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Agyro_aunai/0',Agyro_aunai)
h5create(opath,'/Step#0/Block/Nongyro_swisdak/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Nongyro_swisdak/0',Nongyro_swisdak)
h5create(opath,'/Step#0/Block/Nongyro_aunai/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Nongyro_aunai/0',Nongyro_swisdak)

