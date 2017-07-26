%maxp_common

for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

qom=-256;

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);


file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny,Nz);

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoi=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(rhoi,Nx,Ny,Nz);

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

E= 4*pi*rhoe.*.5.*(Vx.^2+Vy.^2+Vz.^2)./qom;
E= E./(-4*pi*rhoe);
maxE= max(E(:))*code_T


file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vix=fread(fid,'real*8');
fclose(fid);
Vix=reshape(Vix,Nx,Ny,Nz);

file=[dir 'Ji_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Viy=fread(fid,'real*8');
fclose(fid);
Viy=reshape(Viy,Nx,Ny,Nz);

file=[dir 'Ji_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Viz=fread(fid,'real*8');
fclose(fid);
Viz=reshape(Viz,Nx,Ny,Nz);

Vix=Vix./rhoi;
Viy=Viy./rhoi;
Viz=Viz./rhoi;

Ei= 4*pi*rhoi.*.5.*(Vix.^2+Viy.^2+Viz.^2);
Ei= Ei./(4*pi*rhoi);
maxEi= max(Ei(:))*code_T

file=[dir 'PdVi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_par=fread(fid,'real*8');
fclose(fid);
PdVi_par=reshape(PdVi_par,Nx,Ny,Nz);

file=[dir 'PdVi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_per1=fread(fid,'real*8');
fclose(fid);
PdVi_per1=reshape(PdVi_per1,Nx,Ny,Nz);

file=[dir 'PdVi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PdVi_per2=fread(fid,'real*8');
fclose(fid);
PdVi_per2=reshape(PdVi_per2,Nx,Ny,Nz);



file=[dir 'Pe_xy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXY=reshape(V,Nx,Ny,Nz);

file=[dir 'Pe_xz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXZ=reshape(V,Nx,Ny,Nz);

file=[dir 'Pe_yz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PYZ=reshape(V,Nx,Ny,Nz);


file=[dir 'Pe_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXX=reshape(V,Nx,Ny,Nz);

file=[dir 'Pe_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PYY=reshape(V,Nx,Ny,Nz);

file=[dir 'Pe_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PZZ=reshape(V,Nx,Ny,Nz);

P=(PXX+PYY+PZZ)/3;

[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

[FX]=divergence(X,Y,Z,permute(PXX,[2 1 3]),permute(PXY,[2 1 3]),permute(PXZ,[2 1 3])); % Valid only because dx=dy=dz
FX=permute(FX,[2 1 3])./rhoe;

[FY]=divergence(X,Y,Z,permute(PXY,[2 1 3]),permute(PYY,[2 1 3]),permute(PYZ,[2 1 3])); % Valid only because dx=dy=dz
FY=permute(FY,[2 1 3])./rhoe;

[FZ]=divergence(X,Y,Z,permute(PXZ,[2 1 3]),permute(PYZ,[2 1 3]),permute(PZZ,[2 1 3])); % Valid only because dx=dy=dz
FZ=permute(FZ,[2 1 3])./rhoe;


end

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);

w=1.0./sqrt(Bx(i,:,k).^2+1e-10);

ymax(i,k)=round(sum(w.*y)./sum(w));

nav=0;
jrange=round(ymax(i,k))-nav:round(ymax(i,k))+nav;

Tepar_max(i,k)=Tepar(i,round(ymax(i,k)),k);
Teper1_max(i,k)=Teper1(i,round(ymax(i,k)),k);
Teper2_max(i,k)=Teper2(i,round(ymax(i,k)),k);
E_max(i,k)=E(i,round(ymax(i,k)),k);
Ei_max(i,k)=Ei(i,round(ymax(i,k)),k);
Vx_max(i,k)=Vx(i,round(ymax(i,k)),k);
Vy_max(i,k)=Vy(i,round(ymax(i,k)),k);
Vz_max(i,k)=Vz(i,round(ymax(i,k)),k);

Fx_max(i,k)=mean(FX(i,jrange,k));
Fy_max(i,k)=mean(FY(i,jrange,k));
Fz_max(i,k)=mean(FZ(i,jrange,k));

Vix_max(i,k)=Vix(i,round(ymax(i,k)),k);
Viy_max(i,k)=Viy(i,round(ymax(i,k)),k);
Viz_max(i,k)=Viz(i,round(ymax(i,k)),k);
PdVi_par_max(i,k) = PdVi_par(i,round(ymax(i,k)),k);
PdVi_per1_max(i,k) = PdVi_per1(i,round(ymax(i,k)),k);
PdVi_per2_max(i,k) = PdVi_per2(i,round(ymax(i,k)),k);
end
end

maxE= max(E_max(:))*code_T


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tepar_max*code_T,['Tepar' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E','Tepar [keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper1_max*code_T,['Teper1' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E', 'Teperp1 [keV]' )
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper2_max*code_T,['Teper2' ncycle1],[5 15],3,ncycle1,[],3,'x/R_E','y/R_E', 'Teperp2 [keV]')

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),2.*(Teper2_max-Teper1_max)./(Teper2_max+Teper1_max),['Teagyr' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E', 'TeAgyr')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tepar_max./Teper1_max,['Teanis' ncycle1],[0.1 2],3,ncycle1,[],3,'x/R_E','y/R_E', 'TeAnis')


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(E_max*code_T),['E' ncycle1],[-1 0],3,ncycle1,[],3,'x/R_E','y/R_E','log10(Tebulk) [keV]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vx_max*code_V,['Vex' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vex(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vy_max*code_V,['Vez' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vez(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vz_max*code_V,['Vey' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vey(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vix_max*code_V,['Vix' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vix(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viy_max*code_V,['Viz' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Viz(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viz_max*code_V,['Viy' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Viy(gsm) [Km/s]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(Ei_max*code_T),['Ei' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','log10(Tibulk) [keV]')

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_par_max,['PdVi_par' ncycle1],[-2 2]*1e-9,3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{||} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_per1_max,['PdVi_per1' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{\perp 1} [code]')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_per2_max,['PdVi_per1' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{\perp 2} [code]')

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Fx_max/B0^2,['Ex' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ex(gsm) /B_0/V_A')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Fy_max/B0^2,['Ez' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ez(gsm) /B_0/V_A')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Fz_max/B0^2,['Ey' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ey(gsm) /B_0/V_A')

end
        

