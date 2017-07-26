maxp_common

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

B2=Bx.^2+By.^2+Bz.^2;

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);


file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny,Nz);

file=[dir 'Pi_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper1=fread(fid,'real*8');
fclose(fid);
Peper1=reshape(Peper1,Nx,Ny,Nz);
Teper1=Peper1./(rhoe);

file=[dir 'Pi_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper2=fread(fid,'real*8');
fclose(fid);
Peper2=reshape(Peper2,Nx,Ny,Nz);
Teper2=Peper2./(rhoe);

file=[dir 'Pi_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pepar=fread(fid,'real*8');
fclose(fid);
Pepar=reshape(Pepar,Nx,Ny,Nz);
Tepar=Pepar./(rhoe);

file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny,Nz);

file=[dir 'Ji_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vy=fread(fid,'real*8');
fclose(fid);
Vy=reshape(Vy,Nx,Ny,Nz);

file=[dir 'Ji_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vz=fread(fid,'real*8');
fclose(fid);
Vz=reshape(Vz,Nx,Ny,Nz);

Vx=Vx./rhoe;
Vy=Vy./rhoe;
Vz=Vz./rhoe;

adiabatic=Teper1./B2;

E= 4*pi*rhoe.*.5.*(Vx.^2+Vy.^2+Vz.^2)./qom;
E= E./(4*pi*rhoe);
maxE= max(E(:))*code_T

work=0
if(work)
file=[dir 'PdVi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_x=fread(fid,'real*8');
fclose(fid);
PdVi_x=reshape(PdVi_x,Nx,Ny,Nz);

file=[dir 'PdVi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_y=fread(fid,'real*8');
fclose(fid);
PdVi_y=reshape(PdVi_y,Nx,Ny,Nz);

file=[dir 'PdVi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_z=fread(fid,'real*8');
fclose(fid);
PdVi_z=reshape(PdVi_z,Nx,Ny,Nz);


file=[dir 'PxBi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_x=fread(fid,'real*8');
fclose(fid);
PxBi_x=reshape(PxBi_x,Nx,Ny,Nz);

<<<<<<< HEAD
navg=1;
iy=round(ymax(i,k))-navg:round(ymax(i,k))+navg;
=======
file=[dir 'PxBi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_y=fread(fid,'real*8');
fclose(fid);
PxBi_y=reshape(PxBi_y,Nx,Ny,Nz);

file=[dir 'PxBi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_z=fread(fid,'real*8');
fclose(fid);
PxBi_z=reshape(PxBi_z,Nx,Ny,Nz);
>>>>>>> e59a8656edc19a4f737953264983d475c21a3237


file=[dir 'divPui_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
divPui_x=fread(fid,'real*8');
fclose(fid);
divPui_x=reshape(divPui_x,Nx,Ny,Nz);

file=[dir 'divPui_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
divPui_y=fread(fid,'real*8');
fclose(fid);
divPui_y=reshape(divPui_y,Nx,Ny,Nz);

file=[dir 'divPui_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
divPui_z=fread(fid,'real*8');
fclose(fid);
divPui_z=reshape(divPui_z,Nx,Ny,Nz);


file=[dir 'WorkPi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_x=fread(fid,'real*8');
fclose(fid);
WorkPi_x=reshape(WorkPi_x,Nx,Ny,Nz);

file=[dir 'WorkPi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_y=fread(fid,'real*8');
fclose(fid);
WorkPi_y=reshape(WorkPi_y,Nx,Ny,Nz);

file=[dir 'WorkPi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_z=fread(fid,'real*8');
fclose(fid);
WorkPi_z=reshape(WorkPi_z,Nx,Ny,Nz);

<<<<<<< HEAD
An_max(i,k)=log10(lamb(3)./lamb(1));
A2n_max(i,k)=log10(lamb(2)./lamb(1));
Pr=log10(lamb(3).*lamb(1)./lamb(2).^2);
Bi_max(i,k)=Pr./An_max(i,k);
=======
>>>>>>> e59a8656edc19a4f737953264983d475c21a3237

file=[dir 'Enthi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_x=fread(fid,'real*8');
fclose(fid);
Enthi_x=reshape(Enthi_x,Nx,Ny,Nz);

file=[dir 'Enthi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_y=fread(fid,'real*8');
fclose(fid);
Enthi_y=reshape(Enthi_y,Nx,Ny,Nz);

file=[dir 'Enthi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_z=fread(fid,'real*8');
fclose(fid);
Enthi_z=reshape(Enthi_z,Nx,Ny,Nz);
end

end

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Teper1(i,:,k)+Tepar(i,:,k)+Teper2(i,:,k);
%w=Pper1(i,:,k);
ymax(i,k)=round(sum(w.*y)./sum(w));
Tepar_max(i,k)=Tepar(i,round(ymax(i,k)),k);
Teper1_max(i,k)=Teper1(i,round(ymax(i,k)),k);
Teper2_max(i,k)=Teper2(i,round(ymax(i,k)),k);
E_max(i,k)=E(i,round(ymax(i,k)),k);
B2_max(i,k)=B2(i,round(ymax(i,k)),k);
By_max(i,k)=By(i,round(ymax(i,k)),k);

if(work)
PdVi_x_max(i,k) = PdVi_x(i,round(ymax(i,k)),k);
PdVi_y_max(i,k) = PdVi_y(i,round(ymax(i,k)),k);
PdVi_z_max(i,k) = PdVi_z(i,round(ymax(i,k)),k);

PxBi_x_max(i,k) = PxBi_x(i,round(ymax(i,k)),k);
PxBi_y_max(i,k) = PxBi_y(i,round(ymax(i,k)),k);
PxBi_z_max(i,k) = PxBi_z(i,round(ymax(i,k)),k);

divPui_x_max(i,k) = divPui_x(i,round(ymax(i,k)),k);
divPui_y_max(i,k) = divPui_y(i,round(ymax(i,k)),k);
divPui_z_max(i,k) = divPui_z(i,round(ymax(i,k)),k);

WorkPi_x_max(i,k) = WorkPi_x(i,round(ymax(i,k)),k);
WorkPi_y_max(i,k) = WorkPi_y(i,round(ymax(i,k)),k);
WorkPi_z_max(i,k) = WorkPi_z(i,round(ymax(i,k)),k);

Enthi_x_max(i,k) = Enthi_x(i,round(ymax(i,k)),k);
Enthi_y_max(i,k) = Enthi_y(i,round(ymax(i,k)),k);
Enthi_z_max(i,k) = Enthi_z(i,round(ymax(i,k)),k);
end

end
end

adiabatic=(Teper1_max+Teper2_max)./sqrt(B2_max);
badiabatic=By_max./Teper1_max;


maxE= max(E_max(:))*code_T

<<<<<<< HEAD
d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Ani_max,['Ani' ncycle1],[-1 1]*0.4,3,ncycle1,[],3,'x/R_E','y/R_E', 'Anisotropy');
=======
nsmooth=8;
val_range=[-.5 .5]*1e-9;
>>>>>>> e59a8656edc19a4f737953264983d475c21a3237

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tepar_max*code_T,['Ti_x' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E','T_{iX}[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper1_max*code_T,['Ti_y' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E','T_{iY}[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper2_max*code_T,['Ti_z' ncycle1],[20 80],3,ncycle1,[],3,'x/R_E','y/R_E', 'T_{iZ}[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),E_max*code_T,['Ei' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','Tibulk[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(adiabatic),['Adiabatici' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','Adiabatic')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),badiabatic,['AdiBz' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','Adiabatic2')

if(work)
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-PdVi_x_max,['PdVi_x' ncycle1],val_range,nsmooth,ncycle1,[],5,'x/R_E','y/R_E','PdVi_{x} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-PdVi_y_max,['PdVi_y' ncycle1],val_range,nsmooth,ncycle1,[],5,'x/R_E','y/R_E','PdVi_{y} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-PdVi_z_max,['PdVi_z' ncycle1],val_range,nsmooth,ncycle1,[],5,'x/R_E','y/R_E','PdVi_{z} [code]')

<<<<<<< HEAD
d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),An_max,['Lamda3_spread' ncycle1],[-1 1]*0,3,ncycle1,[],3,'x/R_E','y/R_E', 'log_10(\lambda_3 / \lambda_1)');

d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),A2n_max,['Lamda2_spread' ncycle1],[-1 1]*0,3,ncycle1,[],3,'x/R_E','y/R_E', 'log_10(\lambda_2 / \lambda_1)');
=======
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PxBi_x_max,['PxBi_x' ncycle1],val_range,6,ncycle1,[],5,'x/R_E','y/R_E','PxBi_{x} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PxBi_y_max,['PxBi_y' ncycle1],val_range,6,ncycle1,[],5,'x/R_E','y/R_E','PxBi_{y} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PxBi_z_max,['PxBi_z' ncycle1],val_range,6,ncycle1,[],5,'x/R_E','y/R_E','PxBi_{z} [code]')

>>>>>>> e59a8656edc19a4f737953264983d475c21a3237


<<<<<<< HEAD
d=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bi_max,['Biaxiality' ncycle1],[-1 1],3,ncycle1,[],3,'x/R_E','y/R_E', 'Biaxiality');
=======
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-divPui_x_max,['divPui_x' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','divPui_{x} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-divPui_y_max,['divPui_y' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','divPui_{y} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-divPui_z_max,['divPui_z' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','divPui_{z} [code]')
>>>>>>> e59a8656edc19a4f737953264983d475c21a3237

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Enthi_x_max,['Enthi_x' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','Enthi_{x} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Enthi_y_max,['Enthi_y' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','Enthi_{y} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Enthi_z_max,['Enthi_z' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','Enthi_{z} [code]')

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),WorkPi_x_max,['WorkPi_x' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','WorkPi_{x} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),WorkPi_y_max,['WorkPi_y' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','WorkPi_{y} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),WorkPi_z_max,['WorkPi_z' ncycle1],val_range,nsmooth,ncycle1,[],15,'x/R_E','y/R_E','WorkPi_{z} [code]')
end

end
        

