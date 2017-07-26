%maxp_common

cuttone_Ey=[];
cuttone_Vex=[];
cuttone_Vix=[];

for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)

file=[dir 'Pi_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pxx=fread(fid,'real*8');
fclose(fid);
Pxx=reshape(Pxx,Nx,Ny,Nz);

file=[dir 'Pi_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pyy=fread(fid,'real*8');
fclose(fid);
Pyy=reshape(Pyy,Nx,Ny,Nz);

file=[dir 'Pi_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pzz=fread(fid,'real*8');
fclose(fid);
Pzz=reshape(Pzz,Nx,Ny,Nz);

P=Pxx+Pyy+Pzz;

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

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);

file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vex=fread(fid,'real*8');
fclose(fid);
Vex=reshape(Vex,Nx,Ny,Nz);

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vey=fread(fid,'real*8');
fclose(fid);
Vey=reshape(Vey,Nx,Ny,Nz);

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vez=fread(fid,'real*8');
fclose(fid);
Vez=reshape(Vez,Nx,Ny,Nz);

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny,Nz);
Vex=Vex./rhoe;
Vey=Vey./rhoe;
Vez=Vez./rhoe;

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

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoi=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(rhoi,Nx,Ny,Nz);
Vix=Vix./rhoi;
Viy=Viy./rhoi;
Viz=Viz./rhoi;

end

[nx ny nz]= size(P)
y=1:ny;
for i=1:nx
for k=1:nz
w=P(i,:,k);
ymax(i,k)=round(sum(w.*y)./sum(w));
%ymax(i,k)=round(sum(y./V(i,:,k).^2)./sum(1./V(i,:,k).^2));
%[dum j] = min(Vx(i,:,k).^2);
%ymax(i,k)=j;
Bxmax(i,k)=Bx(i,round(ymax(i,k)),k);
%Vxmax(i,k)=Vx(i,end/2,k);
Bymax(i,k)=By(i,round(ymax(i,k)),k);
Bzmax(i,k)=Bz(i,round(ymax(i,k)),k);
Ezmax(i,k)=Ez(i,round(ymax(i,k)),k);
Vexmax(i,k)=Vex(i,round(ymax(i,k)),k);
Veymax(i,k)=Vey(i,14+round(ymax(i,k)),k); %14 is oine RE above
Vezmax(i,k)=Vez(i,round(ymax(i,k)),k);
Vixmax(i,k)=Vix(i,round(ymax(i,k)),k);
Viymax(i,k)=Viy(i,14+round(ymax(i,k)),k); %14 is oine RE above
Vizmax(i,k)=Viz(i,round(ymax(i,k)),k);
end
end

ixcut=180
ixcut=round(Nx*4/6)
ixcut=round(Nx/2)
ixcut2=round(Nx/4)
ixcut2=ixcut

tmp=smooth(Ezmax,3);
cuttone_Ey=[cuttone_Ey; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

tmp=smooth(Vexmax,3);
cuttone_Vex=[cuttone_Vex; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

tmp=smooth(Vixmax,3);
cuttone_Vix=[cuttone_Vix; mean(tmp(round(ixcut2-10:ixcut2+10),:),1)];


%subplot(3,1,1)
%surf(ymax)
figure(1)
%pcolor(ymax')
%shading interp
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),gsmy2z(ymax*dy),['ZGSMmax' ncycle1],[-3.2 -2.7],5,ncycle1,[],1,'x/R_E','y/R_E','Zgsm [RE]')
colorbar
title(['time (s) = ' ntime])

load('ymax0_HRmaha3D1.mat')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),((ymax-ymax0)*dy),['deltaycode' ncycle1],[-1 1],5,ncycle1,[],99,'x/R_E','y/R_E','\deltaZgsm [RE]')

i1=10
i2=Nx-10
j1=10
j2=Nz-10
%subplot(3,1,2)
%pcolor(Ezmax')
immagine_xy(gsmx([xc(i1) xc(i2)]),gsmz2y([zc(j1) zc(j2)]),Ezmax(i1:i2,j1:j2)/B0^2,['EYgsm' ncycle1] ,[-0.4 0.4],0,ncycle1,[],2,'x/R_E','y/R_E','EYgsm /B_0/V_A')
%shading interp
%caxis([-2 2]*1e-6)
colorbar
%subplot(3,1,3)
%pcolor(Bymax')
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bymax*code_B,['BZgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bx[nT]')
colorbar
%caxis([-4 4]*1e-4)
%shading interp
%set(gcf,'Renderer','zbuffer');
%print -dpng caz.png
       
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vizmax*code_V,['ViYgsm' ncycle1],[-600 600],0,ncycle1,[],3,'x/R_E','y/R_E','Vy[km/s]')
colorbar
        
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vixmax*code_V,['ViXgsm' ncycle1],[-600 600],0,ncycle1,[],3,'x/R_E','y/R_E','Vx [km/s]')
colorbar
 
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viymax*code_V,['ViZgsm' ncycle1],[-600 600],0,ncycle1,[],3,'x/R_E','y/R_E','Vz [km/s]')
colorbar
        
        val_range = [-2000 2000];
        immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vezmax*code_V,['VeYgsm' ncycle1],val_range,0,ncycle1,[],3,'x/R_E','y/R_E','Vy[km/s]')
        colorbar
        
        immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vexmax*code_V,['VeXgsm' ncycle1],val_range,0,ncycle1,[],3,'x/R_E','y/R_E','Vx [km/s]')
        colorbar
        
        immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Veymax*code_V,['VeZgsm' ncycle1],val_range,0,ncycle1,[],3,'x/R_E','y/R_E','Vz [km/s]')
        colorbar

        
figure(4)
subplot(2,1,1)
plot(mean(Ezmax,2)/B0^2)
subplot(2,1,2)
plot(mean(Bymax,2))
grid on
print('-dpng',['AVGalongcodez' ncycle1]) 
end

figure
[nt,ndummy]=size(cuttone_Ey);

tc=linspace(0, 154,nt);
pcolor(tc,gsmz2y(yc(j1:j2)),flipud(cuttone_Ey(:,j1:j2)')./B0^2)
colorbar
title(num2str(gsmx(xc(ixcut))),'fontsize',14)
xlabel('t (s)','fontsize',14)
ylabel('Y(gsm)','fontsize',14)
set(gca,'ydir','Reverse')
shading interp
set(gcf,'Renderer','zbuffer');
print('-dpng', ['cuttoneEy' num2str(ixcut) '.png'])

figure


pcolor(tc,gsmz2y(yc(j1:j2)),-flipud(cuttone_Vex(:,j1:j2)')*code_V)
colorbar
title(num2str(gsmx(xc(ixcut))),'fontsize',14)
xlabel('t (s)','fontsize',14)
ylabel('Y(gsm)','fontsize',14)
set(gca,'ydir','Reverse')
shading interp
set(gcf,'Renderer','zbuffer');
print('-dpng', ['cuttoneVex' num2str(ixcut) '.png'])

pcolor(tc,gsmz2y(yc(j1:j2)),-flipud(cuttone_Vix(:,j1:j2)')*code_V)
colorbar
title(num2str(gsmx(xc(ixcut2))),'fontsize',14)
xlabel('t (s)','fontsize',14)
ylabel('Y(gsm)','fontsize',14)
set(gca,'ydir','Reverse')
shading interp
set(gcf,'Renderer','zbuffer');
print('-dpng', ['cuttoneVix' num2str(ixcut2) '.png'])


%!mkdir maxp_frame
%!mv *.png maxp_frame
