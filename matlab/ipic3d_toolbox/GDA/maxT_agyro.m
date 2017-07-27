
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
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
%w=1.0./sqrt(B2(i,:,k).^2+1e-10);

ymax(i,k)=round(sum(w.*y)./sum(w));

iy=round(ymax(i,k))-5:round(ymax(i,k))+5;
p(1,1)=mean(PeXX(i,iy,k));
p(1,2)=mean(PeXY(i,iy,k));
p(1,3)=mean(PeXZ(i,iy,k));
p(2,2)=mean(PeYY(i,iy,k));
p(2,3)=mean(PeYZ(i,iy,k));
p(3,3)=mean(PeZZ(i,iy,k));
p(2,1)=p(1,2);
p(3,1)=p(1,3);
p(3,2)=p(2,3);

b(1)=mean(Bx(i,iy,k));
b(2)=mean(By(i,iy,k));
b(3)=mean(Bz(i,iy,k));
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
lambda1(i,k)=lamb(1);
lambda2(i,k)=lamb(2);
lambda3(i,k)=lamb(3);

Agyro(i,k)= 2*(lamb(3)-lamb(2))/(lamb(3)+lamb(2));


%%%%%%%%%%%
% Aunai
%%%%%%%%%%%

Tr=(p(1,1)+p(2,2)+p(3,3));
Ppar=b'*p*b;
Pper=(Tr-Ppar)/2;
G=eye(3,3)*Pper+(Ppar-Pper)*kron(b,b');
N=p-G;
Agyro_aunai(i,k)=sqrt(sum(N(:).^2))./Tr;

%%%%%%%%%%%
% Swisdak
%%%%%%%%%%%

I2=p(1,1)*p(2,2)+p(1,1)*p(3,3)+p(2,2)*p(3,3);
I2=I2-(p(1,2).^2+p(1,3).^2+p(2,3).^2);
Q=1-4*I2./((Tr-Ppar).*(Tr+3*Ppar));
Nongyro_swisdak(i,k)=sqrt(Q);
% The following formula form Swisdak paper is actually wrong
%Nongyro_aunai(i,k)=sqrt(8*(p(1,2).^2+p(1,3).^2+p(2,3).^2))./(Ppar+2*Pper);
                                 
end
end

                      
global color_choice symmetric_color initial_time Dt
symmetric_color=0;
color_choice=0;
                                 
l1sm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),lambda1,['Lambda_1' ncycle1],[0 0],6,ncycle1,[],3,'x/R_E','y/R_E','Lambda_1');
l2sm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),lambda2,['Lambda_2' ncycle1],[0 3]*1e-7,6,ncycle1,[],3,'x/R_E','y/R_E','Lambda_2');
l3sm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),lambda3,['Lambda_3' ncycle1],[0 3]*1e-7,6,ncycle1,[],3,'x/R_E','y/R_E','Lambda_3');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),l3sm-l2sm,['LambdaDiff' ncycle1],[0 1]*1e-7,0,ncycle1,[],3,'x/R_E','y/R_E','LambdaDiff');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),2*(l3sm-l2sm)./(l3sm+l2sm),['LambdaRatio' ncycle1],[.1 .3],0,ncycle1,[],3,'x/R_E','y/R_E','LambdaRatio');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Agyro,['Agyro' ncycle1],[0 .4],3,ncycle1,[],3,'x/R_E','y/R_E','Agyro');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Agyro_aunai,['Agyro_Aunai' ncycle1],[0 .2],3,ncycle1,[],3,'x/R_E','y/R_E','Agyro');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Nongyro_swisdak,['Nongyro_Swisdak' ncycle1],[0 .2],3,ncycle1,[],3,'x/R_E','y/R_E','NonGyro');
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Nongyro_aunai,['Nongyro_Aunai' ncycle1],[0 .2],3,ncycle1,[],3,'x/R_E','y/R_E','NonGyro');

end
        

