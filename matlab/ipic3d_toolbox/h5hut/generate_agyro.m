close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is
%dir='/data1/gianni/HRmaha3D3/vtk/'; %directory where the files are


must_read=true;
leggo='h5'
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

