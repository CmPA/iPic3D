
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

[nx ny nz]= size(PeXX)

for Ygsm=Ygsmrange(1)+1:.3:Ygsmrange(2)-1
kc = round(Nz * (Ygsm-Ygsmrange(1)) / (Ygsmrange(2)-Ygsmrange(1)) )
navg=10;
k=kc-navg:kc+navg;

for i=1:nx
for iy=1:ny


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
lambda1(i,iy)=lamb(1);
lambda2(i,iy)=lamb(2);
lambda3(i,iy)=lamb(3);

Agyro(i,iy)= 2*(lamb(3)-lamb(2))/(lamb(3)+lamb(2));


%%%%%%%%%%%
% Aunai
%%%%%%%%%%%

Tr=(p(1,1)+p(2,2)+p(3,3));
Ppar=b'*p*b;
Pper=(Tr-Ppar)/2;
G=eye(3,3)*Pper+(Ppar-Pper)*kron(b,b');
N=p-G;
Agyro_aunai(i,iy)=sqrt(sum(N(:).^2))./Tr;


%%%%%%%%%%%
% Swisdak
%%%%%%%%%%%

I2=p(1,1)*p(2,2)+p(1,1)*p(3,3)+p(2,2)*p(3,3);
I2=I2-(p(1,2).^2+p(1,3).^2+p(2,3).^2);
Q=1-4*I2./((Tr-Ppar).*(Tr+3*Ppar));
Nongyro_swisdak(i,iy)=sqrt(Q);

end
end

%Agyro_sm=smooth(Agyro,3);

x=fliplr(Xgsmrange);
y=Zgsmrange;

global blowup contours
blowup=0;
contours=1;


immagine(x,y,Agyro,['Agyro_' ncycle1 '_Ygsm_' num2str(Ygsm)],[0 .2],3,ncycle1, Ygsm)

if(contours)
Bxm=mean(Bx(:,:,kc),3);
Bym=mean(By(:,:,kc),3);

Ay=zeros(size(Bx));

Ay=vecpot_uniform(xc,yc,Bxm,Bym);
hold on
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
end
  
print('-dpng','-r300',['Agyro_combo_' ncycle1 '_Ygsm_' num2str(Ygsm) '.png'])
 
 
immagine(x,y,Agyro_aunai,['Agyro_Aunai_' ncycle1 '_Ygsm_' num2str(Ygsm)],[0 .1],3,ncycle1, Ygsm)
hold on
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
print('-dpng','-r300',['Agyro_Aunai_combo_' ncycle1 '_Ygsm_' num2str(Ygsm) '.png'])
  
immagine(x,y,Nongyro_swisdak,['Nongyro_Swisdak_' ncycle1 '_Ygsm_' num2str(Ygsm)],[0 .1],3,ncycle1, Ygsm)
hold on
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
print('-dpng','-r300',['Nongyro_Swisdak_combo_' ncycle1 '_Ygsm_' num2str(Ygsm) '.png'])
  
end 
   
end
        

