%maxp_common

for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

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

B2=Bx.^2+By.^2+Bz.^2;

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny,Nz);

file=[dir 'Pi_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pxx=reshape(V,Nx,Ny,Nz);
Txx=Pxx./(rhoe);

file=[dir 'Pi_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pyy=reshape(V,Nx,Ny,Nz);
Tyy=Pyy./(rhoe);

file=[dir 'Pi_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pzz=reshape(V,Nx,Ny,Nz);
Tzz=Pzz./(rhoe);


file=[dir 'Pi_xy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pxy=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_xz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pxz=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_yz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pyz=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Ppar=reshape(V,Nx,Ny,Nz);
Tpar=Ppar./rhoe;

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(V,Nx,Ny,Nz);
Tper1=Pper1./rhoe;

file=[dir 'Pi_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Pper2=reshape(V,Nx,Ny,Nz);
Tper2=Pper2./rhoe;

[nx ny nz]= size(Bx)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pxx(i,:,k)+Pyy(i,:,k)+Pzz(i,:,k);

ymax(i,k)=round(sum(w.*y)./sum(w));

navg=10;
iy=round(ymax(i,k))-navg:round(ymax(i,k))+navg;

Txx_max(i,k)=mean(Txx(i,iy,k));
Tyy_max(i,k)=mean(Tyy(i,iy,k));
Tzz_max(i,k)=mean(Tzz(i,iy,k));

Bx_max(i,k)=mean(Bx(i,iy,k));
By_max(i,k)=mean(By(i,iy,k));
Bz_max(i,k)=mean(Bz(i,iy,k));
B_max(i,k)=sqrt(mean(Bx(i,iy,k).^2+By(i,iy,k).^2+Bz(i,iy,k).^2));

TperM=.5*mean(Tper1(i,iy,k)+Tper2(i,iy,k));
TM=(Txx_max(i,k)+Tyy_max(i,k)+Tzz_max(i,k))/3;


Ani_max(i,k)=mean(Tpar(i,iy,k)-TperM)./TM;
Agy_max(i,k)=mean(Tper1(i,iy,k)-Tper2(i,iy,k))/TperM;

p(1,1)=mean(Pxx(i,iy,k));
p(1,2)=mean(Pxy(i,iy,k));
p(1,3)=mean(Pxz(i,iy,k));
p(2,2)=mean(Pyy(i,iy,k));
p(2,3)=mean(Pyz(i,iy,k));
p(3,3)=mean(Pzz(i,iy,k));
p(2,1)=p(1,2);
p(3,1)=p(1,3);
p(3,2)=p(2,3);

[v,e]=eig(p);

lamb=diag(e);
lambda1(i,k)=lamb(1);
lambda2(i,k)=lamb(2);
lambda3(i,k)=lamb(3);

An_max(i,k)=log10(lamb(3)./lamb(1));
A2n_max(i,k)=log10(lamb(2)./lamb(1));
Pr_max(i,k)=log10(lamb(3).*lamb(1)./lamb(2).^2);
Bi=Pr_max(i,k)./An_max(i,k);

w1=v(:,1);
w2=v(:,2);
w3=v(:,3);

b=[Bx_max(i,k); By_max(i,k); Bz_max(i,k)];

Mis_max(i,k)=min([norm(cross(w1,b))./norm(b)./norm(w1) norm(cross(w2,b))./norm(b)./norm(w2) norm(cross(w3,b))./norm(b)./norm(w3) ]);
end
end

global color_choice

color_choice=0

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Txx_max*code_T,['Ti_x' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E','T_{iX}[keV]');

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tyy_max*code_T,['Ti_y' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E','T_{iY}[keV]');

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tzz_max*code_T,['Ti_z' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E', 'T_{iZ}[keV]');

color_choice=1

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Ani_max,['Ani' ncycle1],[-1 1]*0.3,3,ncycle1,[],3,'x/R_E','y/R_E', 'Anisotropy');

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Agy_max,['Agy' ncycle1],[-1 1]*0.3,3,ncycle1,[],3,'x/R_E','y/R_E', 'Agyrotropy');

color_choice=0

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),An_max,['Lamda3_spread' ncycle1],[-1 1]*0,3,ncycle1,[],3,'x/R_E','y/R_E', '\lambda_3 / \lambda_1');

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),A2n_max,['Lamda2_spread' ncycle1],[-1 1]*0,3,ncycle1,[],3,'x/R_E','y/R_E', '\lambda_2 / \lambda_1');

color_choice=1

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Pr_max,['PrObla' ncycle1],[-1 1]*0.1,3,ncycle1,[],3,'x/R_E','y/R_E', 'Pro/Obla');


color_choice=0
d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Mis_max,['Misalignment' ncycle1],[0 1],3,ncycle1,[],3,'x/R_E','y/R_E', 'Misalignment');

end
        

