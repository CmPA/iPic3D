%maxp_common

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


file=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny,Nz);

file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny,Nz);

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);

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

Vexbx = (Ey.* Bz - Ez.* By)./B2;
Vexby = (Ez.* Bx - Ex.* Bz)./B2;
Vexbz = (Ex.* By - Ey.* Bx)./B2;

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);

close all
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


JdotE = Vx.*Ex+ Vy.*Ey + Vz.*Ez;

Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;

VparB= (Vx.*Bx + Vy.*By + Vz.*Bz);

Vparx=VparB.*Bx./B2;
Vpary=VparB.*By./B2;
Vparz=VparB.*Bz./B2;

Vslipx = Vx - Vexbx;
Vslipy = Vy - Vexby;
Vslipz = Vz - Vexbz;

Vslipperpx = Vx - Vexbx- Vparx;
Vslipperpy = Vy - Vexby- Vpary;
Vslipperpz = Vz - Vexbz- Vparz;
end

[nx ny nz]= size(Vx);
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
%w=1.0./sqrt(Bx(i,:,k).^2+1e-10);
ymax(i,k)=round(sum(w.*y)./sum(w));
JdotE_max(i,k)=JdotE(i,round(ymax(i,k)),k);


Vslipx_max(i,k)=Vslipx(i,round(ymax(i,k)),k);
Vslipy_max(i,k)=Vslipy(i,round(ymax(i,k)),k);
Vslipz_max(i,k)=Vslipz(i,round(ymax(i,k)),k);

Vslipperpx_max(i,k)=Vslipperpx(i,round(ymax(i,k)),k);
Vslipperpy_max(i,k)=Vslipperpy(i,7+round(ymax(i,k)),k);
Vslipperpz_max(i,k)=Vslipperpz(i,round(ymax(i,k)),k);

end
end

Vslipperp_max = sqrt(Vslipperpx_max.^2+Vslipperpy_max.^2+Vslipperpz_max.^2);

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),log10(abs(JdotE)),['logJdotEI' ncycle1],[0 0],3,ncycle1,[],3,'x/R_E','y/R_E','log(JdotE [code units]');

VSMX=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vslipx_max*code_V,['VslipIx' ncycle1],[-.5 .5]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vx(gsm) [Km/s]');
VSMY=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipy_max*code_V,['VslipIz' ncycle1],[-.5 .5]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vz(gsm) [Km/s]');
VSMZ=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipz_max*code_V,['VslipIy' ncycle1],[-.5 .5]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vy(gsm) [Km/s]');

VPSMX=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vslipperpx_max*code_V,['VslipIperpx' ncycle1],[-5 5]*1e2,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp Vx(gsm) [Km/s]');
VPSMY=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperpy_max*code_V,['VslipIperpz' ncycle1],[-5 5]*1e2,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_perp Vz(gsm) [Km/s]');
VPSMZ=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperpz_max*code_V,['VslipIperpy' ncycle1],[-5 5 ]*1e2,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp Vy(gsm) [Km/s]');

VPSM=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperp_max*code_V,['VslipIperp' ncycle1],[0 5 ]*1e2,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp V_i [Km/s]');



end
        

