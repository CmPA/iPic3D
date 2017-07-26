for cycle=Ncyc_ini:1000:Ncyc_max

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

biskamp=0
if(biskamp)
file=[dir 'Biskamp_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bisx=fread(fid,'real*8');
fclose(fid);
Bisx=reshape(Bisx,Nx,Ny,Nz);

file=[dir 'Biskamp_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bisy=fread(fid,'real*8');
fclose(fid);
Bisy=reshape(Bisy,Nx,Ny,Nz);

file=[dir 'Biskamp_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bisz=fread(fid,'real*8');
fclose(fid);
Bisz=reshape(Bisz,Nx,Ny,Nz);
end

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

[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

[Jx, Jy, Jz]=curl(X,Y,Z,permute(Bx,[2 1 3]),permute(By,[2 1 3]),permute(Bz,[2 1 3])); % Valid only because dx=dy=dz
Jx=permute(Jx,[2 1 3])/4/pi;
Jy=permute(Jy,[2 1 3])/4/pi;
Jz=permute(Jz,[2 1 3])/4/pi;


B2=Bx.^2+By.^2+Bz.^2;

[Vexbx, Vexby, Vexbz] = cross_overmod2(Ex, Ey, Ez, Bx, By, Bz, B2);


file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoi=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(rhoi,Nx,Ny,Nz);
rhoc=rhoi+rho;

close all
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

close all
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

Vx=Jx-Vix;
Vy=Jy-Viy;
Vz=Jz-Viz;

JedotE = Vx.*Ex+ Vy.*Ey + Vz.*Ez;
JdotE = Jx.*Ex+ Jy.*Ey + Jz.*Ez;

Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;

mi=1;
me=-1/qom;
Vmhdx=(mi*Vix./rhoi+me*Vx)/(me+mi);
Vmhdy=(mi*Viy./rhoi+me*Vy)/(me+mi);
Vmhdz=(mi*Viz./rhoi+me*Vz)/(me+mi);

[Epx, Epy, Epz] = cross_prod(Vmhdx, Vmhdy, Vmhdz, Bx, By, Bz);
Epx = Ex + Epx;
Epy = Ey + Epy;
Epz = Ez + Epz;

JdotEp= Epx.*Jx +Epy.*Jy +Epz.*Jz;

[Zx, Zy, Zz] = cross_prod(Vx, Vy, Vz, Bx, By, Bz);

Zx = Ex + Zx;
Zy = Ey + Zy;
Zz = Ez + Zz;

EdotB = (Ex.*Bx+Ey.*By+Ez.*Bz);
Eparx= EdotB.*Bx./B2;
Epary= EdotB.*By./B2;
Eparz= EdotB.*Bz./B2;

%[curlZx, curlZy, curlZz]=curl(X,Y,Z,permute(Zx,[2 1 3]),permute(Zy,[2 1 3]),permute(Zz,[2 1 3])); % Valid only because dx=dy=dz

[tmpx, tmpy, tmpz]=curl(X,Y,Z,permute(Eparx,[2 1 3]),permute(Epary,[2 1 3]),permute(Eparz,[2 1 3]));
tmpx=permute(tmpx,[2 1 3]);
tmpy=permute(tmpy,[2 1 3]);
tmpz=permute(tmpz,[2 1 3]);

[DBx, DBy, DBz] = cross_overmod2(Bx, By, Bz, tmpx, tmpy, tmpz,B2);


Zenitani_old = rho.*(Vx.*Zx+Vy.*Zy+Vz.*Zz);

Zenitani = (Jx.*Zx+Jy.*Zy+Jz.*Zz) -rhoc.*(Vx.*Ex+Vy.*Ey+Vz.*Ez);

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

Va=sqrt(B2./4./pi./rhoi);
end

[nx ny nz]= size(Vx)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
%w=1.0./sqrt(Bx(i,:,k).^2+1e-10);
%w=1.0./sqrt(B2(i,:,k)+1e-10);
ymax(i,k)=round(sum(w.*y)./sum(w));

nav=0;
jrange=round(ymax(i,k))-nav:round(ymax(i,k))+nav;
JdotE_max(i,k)=mean(JdotE(i,jrange,k));
JedotE_max(i,k)=mean(JedotE(i,jrange,k));
JdotEp_max(i,k)=mean(JdotEp(i,jrange,k));

Va_max=Va(i,round(ymax(i,k)),k);

Bz_max(i,k)=mean(Bz(i,jrange,k));

Bisx_max(i,k)=mean(DBx(i,jrange,k));%./sqrt(B2(i,iy,k)));
Bisy_max(i,k)=mean(DBy(i,jrange,k));%./sqrt(B2(i,iy,k)));
Bisz_max(i,k)=mean(DBz(i,jrange,k));%./sqrt(B2(i,iy,k)));

if(biskamp)
Bisx_max(i,k)=mean(Bisx(i,jrange,k))./(B2(i,jrange,k));
Bisy_max(i,k)=mean(Bisy(i,jrange,k))./(B2(i,jrange,k));
Bisz_max(i,k)=mean(Bisz(i,jrange,k))./(B2(i,jrange,k));
end

Jx_max(i,k)=Jx(i,round(ymax(i,k)),k);
Jy_max(i,k)=Jy(i,round(ymax(i,k)),k);
Jz_max(i,k)=Jz(i,round(ymax(i,k)),k);

Zx_max(i,k)=Zx(i,round(ymax(i,k)),k);
Zy_max(i,k)=Zy(i,round(ymax(i,k)),k);
Zz_max(i,k)=Zz(i,round(ymax(i,k)),k);

%curlZx_max(i,k)=curlZx(i,round(ymax(i,k)),k);
%curlZy_max(i,k)=curlZy(i,round(ymax(i,k)),k);
%curlZz_max(i,k)=curlZz(i,round(ymax(i,k)),k);


Vslipx_max(i,k)=Vslipx(i,round(ymax(i,k)),k);
Vslipy_max(i,k)=Vslipy(i,round(ymax(i,k)),k);
Vslipz_max(i,k)=Vslipz(i,round(ymax(i,k)),k);

Vslipperpx_max(i,k)=Vslipperpx(i,round(ymax(i,k)),k);
Vslipperpy_max(i,k)=Vslipperpy(i,round(ymax(i,k)),k);
Vslipperpz_max(i,k)=Vslipperpz(i,round(ymax(i,k)),k);

%[Zenitani_max(i,k) jjmax]=max(Zenitani(i,jrange,k));
%JdotE_max(i,k)=(JdotE(i,jrange(jjmax),k));
Zenitani_max(i,k)=mean(Zenitani(i,jrange,k));

end
end

J_max = sqrt(Jx_max.^2+0*Jy_max.^2+Jz_max.^2);


Vslip_max = sqrt(Vslipx_max.^2+Vslipy_max.^2+Vslipz_max.^2);
Vslipperp_max = sqrt(Vslipperpx_max.^2+Vslipperpy_max.^2+Vslipperpz_max.^2);

immagine_xyb(gsmx([0 Lx]),gsmz2y([0 Lz]),J_max*code_J,gsmx(xc),gsmz2y(zc),Bz_max,['J' ncycle1],[-2 2]*0e-10,0,ncycle1,[],6,'x/R_E','y/R_E','J [code units]');

end

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bz_max*code_B,['BYgsm' ncycle1],[-1 1]*1e1,0,ncycle1,[],1,'x/R_E','y/R_E','BY(gsm) [nT]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JdotE_max,['JdotE' ncycle1],[-2 2]*1e-10,6,ncycle1,[],6,'x/R_E','y/R_E','JdotE [code units]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JdotEp_max,['JdotEp' ncycle1],[-3 3]*1e-10,6,ncycle1,[],6,'x/R_E','y/R_E','JdotEp [code units]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JedotE_max,['JedotE' ncycle1],[-3 3]*1e-10,6,ncycle1,[],6,'x/R_E','y/R_E','JedotE [code units]');


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zenitani_max,['Zenitani' ncycle1],[-3 3]*1e-10,6,ncycle1,[],12,'x/R_E','y/R_E','Zenitani [code units]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vslipx_max*code_V,['VslipEx' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vx(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipy_max*code_V,['VslipEz' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vz(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipz_max*code_V,['VslipEy' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta Vy(gsm) [Km/s]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslip_max*code_V,['VslipE' ncycle1],[0 2]*1e3,3,ncycle1,[],6,'x/R_E','y/R_E','\delta V [Km/s]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslip_max./Va_max,['VslipEoVa' ncycle1],[0 3],3,ncycle1,[],6,'x/R_E','y/R_E','\delta V/V_A');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vslipperpx_max*code_V,['VslipEperpx' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp Vx(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperpy_max*code_V,['VslipEperpz' ncycle1],[-2 2]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_perp Vz(gsm) [Km/s]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperpz_max*code_V,['VslipEperpy' ncycle1],[-2 2 ]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp Vy(gsm) [Km/s]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperp_max*code_V,['VslipEperp' ncycle1],[0 2 ]*1e3,3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp V_e [Km/s]');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vslipperp_max./Va_max,['VslipEperpoVa' ncycle1],[0 3],3,ncycle1,[],3,'x/R_E','y/R_E','\delta_\perp V_e /V_A');


Zxavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zx_max./B0./Va_max,['OHM_x' ncycle1],[-1 1]*3e-1,6,ncycle1,[],3,'x/R_E','y/R_E',' \Omega_x(gsm)/B_0V_{A}');
Zyavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zy_max./B0./Va_max,['OHM_z' ncycle1],[-1 1]*3e-1,6,ncycle1,[],3,'x/R_E','y/R_E',' \Omega_z(gsm)/B_0V_{A}');

Zzavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zz_max./B0./Va_max,['OHM_y' ncycle1],[-1 1]*3e-1,6,ncycle1,[],3,'x/R_E','y/R_E',' \Omega_y(gsm)/B_0V_{A}');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),sqrt(Zx_max.^2+Zy_max.^2+Zz_max.^2)./B0./Va_max,['OHM_mod' ncycle1],[-1 1]*3e-1,6,ncycle1,[],3,'x/R_E','y/R_E','\Omega_{mod}(gsm)/B_0V_{A}');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),sqrt(Zxavg.^2+Zyavg.^2+Zzavg.^2),['OHM_avg' ncycle1],[0 1]*6e-1,6,ncycle1,[],3,'x/R_E','y/R_E','\Omega_{avg}/B_0V_{A}');


%curlZxavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),curlZx_max*code_E,['curlOHM_x' ncycle1],[0 0]*1e-2,6,ncycle1,[],3,'x/R_E','y/R_E',' \nabla \times\Omega_x(gsm) ');
%curlZyavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),curlZy_max*code_E,['culrOHM_z' ncycle1],[0 0]*1e-2,6,ncycle1,[],3,'x/R_E','y/R_E',' \nabla \times\Omega_z(gsm) ');

%curlZzavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),curlZz_max*code_E,['curlOHM_y' ncycle1],[0 0]*1e-2,6,ncycle1,[],3,'x/R_E','y/R_E',' \nabla \times\Omega_y(gsm) ');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),sqrt(curlZx_max.^2+curlZy_max.^2+curlZz_max.^2)./Va_max.^2,['curlOHM_mod' ncycle1],[0 0]*1e-0,6,ncycle1,[],3,'x/R_E','y/R_E','\nabla \times\Omega_{mod}(gsm) /B_0^2');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),sqrt(curlZxavg.^2+curlZyavg.^2+curlZzavg.^2),['curlOHM_avg' ncycle1],[0 1.0]*1e-2,6,ncycle1,[],3,'x/R_E','y/R_E','\nabla \times\Omega_{avg}(gsm) /B_0^2');


Bisxavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bisx_max,['Biskamp_x' ncycle1],[-1 1]*0e-9,6,ncycle1,[],5,'x/R_E','y/R_E','Biskamp_x');

Bisyavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bisy_max,['Biskamp_y' ncycle1],[-1 1]*0e-9,6,ncycle1,[],5,'x/R_E','y/R_E','Biskamp_y');

Biszavg=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bisz_max,['Biskamp_z' ncycle1],[-1 1]*0e-9,6,ncycle1,[],5,'x/R_E','y/R_E','Biskamp_z');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),sqrt(Bisxavg.^2+Bisyavg.^2+Biszavg.^2),['Biskamp_avg' ncycle1],[0 1]*4e-3,6,ncycle1,[],3,'x/R_E','y/R_E','Biskamp_avg');

%end
        

