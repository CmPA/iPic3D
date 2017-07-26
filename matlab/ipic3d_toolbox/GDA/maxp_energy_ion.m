global color_choice symmetric_color


for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);


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

file=[dir 'Pi_xy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXY=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_xz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXZ=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_yz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PYZ=reshape(V,Nx,Ny,Nz);


file=[dir 'Pi_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PXX=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PYY=reshape(V,Nx,Ny,Nz);

file=[dir 'Pi_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
PZZ=reshape(V,Nx,Ny,Nz);

P=(PXX+PYY+PZZ)/3;

[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

[FX]=divergence(X,Y,Z,permute(PXX,[2 1 3]),permute(PXY,[2 1 3]),permute(PXZ,[2 1 3])); % Valid only because dx=dy=dz
FX=permute(FX,[2 1 3])/4/pi;

[FZ]=divergence(X,Y,Z,permute(PXZ,[2 1 3]),permute(PYZ,[2 1 3]),permute(PZZ,[2 1 3])); % Valid only because dx=dy=dz
FZ=permute(FZ,[2 1 3])/4/pi;


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

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper1=fread(fid,'real*8');
fclose(fid);
Peper1=reshape(Peper1,Nx,Ny,Nz);
Teper1=Peper1./(rhoi);

file=[dir 'Pi_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper2=fread(fid,'real*8');
fclose(fid);
Peper2=reshape(Peper2,Nx,Ny,Nz);
Teper2=Peper2./(rhoi);

file=[dir 'Pi_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pepar=fread(fid,'real*8');
fclose(fid);
Pepar=reshape(Pepar,Nx,Ny,Nz);
Tepar=Pepar./(rhoi);

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

Vx=Vx./rhoi;
Vy=Vy./rhoi;
Vz=Vz./rhoi;

E= 4*pi*rhoi.*.5.*(Vx.^2+Vy.^2+Vz.^2)./qom;
E= E./(4*pi*rhoi);
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

end

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
w=1.0./sqrt(Bx(i,:,k).^2+1e-10);
ymax(i,k)=round(sum(w.*y)./sum(w));

%[dummy ymax(i,k)]=max(w);

navg=10;
iy=ymax(i,k)-navg:ymax(i,k)+navg;
Tepar_max(i,k)=Tepar(i,round(ymax(i,k)),k);
Teper1_max(i,k)=Teper1(i,round(ymax(i,k)),k);
Teper2_max(i,k)=Teper2(i,round(ymax(i,k)),k);
TT=[Tepar(i,round(ymax(i,k)),k) Teper1(i,round(ymax(i,k)),k) Teper2(i,round(ymax(i,k)),k)];
T_std_max(i,k)=std(TT)./mean(TT);
E_max(i,k)=E(i,round(ymax(i,k)),k);
Ei_max(i,k)=Ei(i,round(ymax(i,k)),k);
Vx_max(i,k)=Vx(i,round(ymax(i,k)),k);
Vy_max(i,k)=Vy(i,round(ymax(i,k)),k);
Vz_max(i,k)=Vz(i,round(ymax(i,k)),k);
Vix_max(i,k)=Vix(i,round(ymax(i,k)),k);
Viy_max(i,k)=Viy(i,round(ymax(i,k)),k);
Viz_max(i,k)=Viz(i,round(ymax(i,k)),k);

PXY_max(i,k) = mean(PXY(i,iy,k)./P(i,iy,k));
PXZ_max(i,k) = mean(PXZ(i,iy,k)./P(i,iy,k));
PYZ_max(i,k) = mean(PYZ(i,iy,k)./P(i,iy,k));
FX_max(i,k) = mean(FX(i,iy,k));
FZ_max(i,k) = mean(FZ(i,iy,k));

PdVi_par_max(i,k) = PdVi_par(i,round(ymax(i,k)),k);
PdVi_per1_max(i,k) = PdVi_per1(i,round(ymax(i,k)),k);
PdVi_per2_max(i,k) = PdVi_per2(i,round(ymax(i,k)),k);
end
end

maxE= max(E_max(:))*code_T



immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Tepar_max*code_T,['Tipar' ncycle1],5*[5 15],3,ncycle1,[],3,'x/R_E','y/R_E','Tipar [keV]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper1_max*code_T,['Tiper1' ncycle1],5*[5 15],3,ncycle1,[],3,'x/R_E','y/R_E', 'Tiperp1 [keV]' );
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Teper2_max*code_T,['Tiper2' ncycle1],5*[5 15],3,ncycle1,[],3,'x/R_E','y/R_E', 'Tiperp2 [keV]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),(Teper2_max+Teper1_max)/2*code_T,['Tiper' ncycle1],5*[5 15],3,ncycle1,[],3,'x/R_E','y/R_E', 'Tiperp [keV]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),T_std_max,['Tistd' ncycle1],[0 0],3,ncycle1,[],6,'x/R_E','y/R_E', 'TiSTD');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),2*abs(Teper2_max-Teper1_max)./(Teper2_max+Teper1_max),['Tiagyr' ncycle1],[0 0],3,ncycle1,[],6,'x/R_E','y/R_E', 'TiAgyr');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),2*Tepar_max./(Teper2_max+Teper1_max),['Tianis' ncycle1],[0.5 1.5],6,ncycle1,[],3,'x/R_E','y/R_E', 'TiAnis');

color_choice=1

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PXY_max,['PXZ' ncycle1],[-1 1]*1e-1,6,ncycle1,[],3,'x/R_E','y/R_E', 'P_{iXZ(gsm)}/Tr(P)');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PXZ_max,['PXY' ncycle1],[-1 1]*1e-1,6,ncycle1,[],3,'x/R_E','y/R_E', 'P_{iXY(gsm)}/Tr(P)');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PYZ_max,['PYZ' ncycle1],[-1 1]*1e-1,6,ncycle1,[],3,'x/R_E','y/R_E', 'P_{iYZ(gsm)}/Tr(P)');
                   
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),FX_max,['FX' ncycle1],[-1 1]*.5e-8,10,ncycle1,[],3,'x/R_E','y/R_E', 'DivP_{iX(gsm)}');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),FZ_max,['FY' ncycle1],[-1 1]*.5e-8,10,ncycle1,[],3,'x/R_E','y/R_E', 'DivP_{iY(gsm)}');


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(E_max*code_T),['E' ncycle1],[-1 0],3,ncycle1,[],3,'x/R_E','y/R_E','log10(Tebulk) [keV]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vx_max*code_V,['Vex' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vex(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vy_max*code_V,['Vez' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vez(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vz_max*code_V,['Vey' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vey(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vix_max*code_V,['Vix' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Vix(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viy_max*code_V,['Viz' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Viz(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viz_max*code_V,['Viy' ncycle1],[0 0]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','Viy(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(Ei_max*code_T),['Ei' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','log10(Tibulk) [keV]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_par_max,['PdVi_par' ncycle1],[-2 2]*1e-9,3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{||} [code]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_per1_max,['PdVi_per1' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{\perp 1} [code]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),PdVi_per2_max,['PdVi_per1' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','PdVi_{\perp 2} [code]');

end
        

