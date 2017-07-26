%maxp_common

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


file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny,Nz);

file=[dir 'Pe_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper1=fread(fid,'real*8');
fclose(fid);
Peper1=reshape(Peper1,Nx,Ny,Nz);
Teper1=Peper1./(-rhoe);

file=[dir 'Pe_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper2=fread(fid,'real*8');
fclose(fid);
Peper2=reshape(Peper2,Nx,Ny,Nz);
Teper2=Peper2./(-rhoe);

file=[dir 'Pe_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pepar=fread(fid,'real*8');
fclose(fid);
Pepar=reshape(Pepar,Nx,Ny,Nz);
Tepar=Pepar./(-rhoe);

file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny,Nz);

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vy=fread(fid,'real*8');
fclose(fid);
Vy=reshape(Vy,Nx,Ny,Nz);

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vz=fread(fid,'real*8');
fclose(fid);
Vz=reshape(Vz,Nx,Ny,Nz);

Vx=Vx./rhoe;
Vy=Vy./rhoe;
Vz=Vz./rhoe;

adiabatic=Teper1./B2;

E= 4*pi*rhoe.*.5.*(Vx.^2+Vy.^2+Vz.^2)./qom;
E= E./(-4*pi*rhoe);
maxE= max(E(:))*code_T
end

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Teper1(i,:,k);
ymax(i,k)=round(sum(w.*y)./sum(w));
Tepar_max(i,k)=Tepar(i,round(ymax(i,k)),k);
Teper1_max(i,k)=Teper1(i,round(ymax(i,k)),k);
Teper2_max(i,k)=Teper2(i,round(ymax(i,k)),k);
E_max(i,k)=E(i,round(ymax(i,k)),k);
B2_max(i,k)=B2(i,round(ymax(i,k)),k);
By_max(i,k)=By(i,round(ymax(i,k)),k);
end
end

adiabatic=Teper1_max./sqrt(B2_max);
badiabatic=By_max./Teper1_max;


maxE= max(E_max(:))*code_T


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tepar_max*code_T,['Tepar' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E','T_{e||}[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper1_max*code_T,['Teper1' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E','T_{e\perp 1}[keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper2_max*code_T,['Teper2' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E','T_{e\perp 2}[keV]')
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(E_max*code_T),['E' ncycle1],[-1 0],3,ncycle1,[],3,'x/R_E','y/R_E')
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(adiabatic),['Adiabatic' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E')
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),badiabatic,['AdBz' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E')

end
        

