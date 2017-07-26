
for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);

B2=sqrt(Bx.^2+By.^2+Bz.^2);
Bx=Bx./B2;
By=By./B2;
Bz=Bz./B2;

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);

file=[dir 'Pe_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeXX=fread(fid,'real*8');
fclose(fid);
PeXX=reshape(PeXX,Nx,Ny,Nz);

file=[dir 'Pe_xy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeXY=fread(fid,'real*8');
fclose(fid);
PeXY=reshape(PeXY,Nx,Ny,Nz);

file=[dir 'Pe_xz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeXZ=fread(fid,'real*8');
fclose(fid);
PeXZ=reshape(PeXZ,Nx,Ny,Nz);

file=[dir 'Pe_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeYY=fread(fid,'real*8');
fclose(fid);
PeYY=reshape(PeYY,Nx,Ny,Nz);

file=[dir 'Pe_yz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeYZ=fread(fid,'real*8');
fclose(fid);
PeYZ=reshape(PeYZ,Nx,Ny,Nz);

file=[dir 'Pe_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PeZZ=fread(fid,'real*8');
fclose(fid);
PeZZ=reshape(PeZZ,Nx,Ny,Nz);

end

[nx ny nz]= size(Pper1)

for i=1:nx
for iy=1:ny
for k=1:nz

p(1,1)=(PeXX(i,iy,k));
p(1,2)=(PeXY(i,iy,k));
p(1,3)=(PeXZ(i,iy,k));
p(2,2)=(PeYY(i,iy,k));
p(2,3)=(PeYZ(i,iy,k));
p(3,3)=(PeZZ(i,iy,k));
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
savevtk_bin(Agyro_sm,['Agyro' ncycle1 '.vtk'],'Agyro',dx,dy,dz)
savevtk_bin(Agyro_aunai_sm,['Agyro_aunai' ncycle1 '.vtk'],'Agyro',dx,dy,dz)
savevtk_bin(Nongyro_swisdak_sm,['Nongyro_swisdak' ncycle1 '.vtk'],'Agyro',dx,dy,dz)


%savevtk(lambda1,['lambda1' ncycle1 '.vtk'],'lambda1',dx,dy,dz)
%savevtk(lambda2,['lambda2' ncycle1 '.vtk'],'lambda2',dx,dy,dz)
%savevtk(lambda3,['lambda3' ncycle1 '.vtk'],'lambda3',dx,dy,dz)


%for iz=1:Nz
%immagine(x,y,Agyro_sm(:,:,iz),['Agyro_' num2str(iz) '_' ncycle1],[0 .4],3,ncycle1, Ygsm)
%end


end
        

