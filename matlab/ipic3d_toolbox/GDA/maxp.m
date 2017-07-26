
cuttone=[];
cuttone2=[];

for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
if(exist(file)==2)
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);
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

Ej=(Ex.*Bz-Ez.*Bx)./sqrt(Bx.^2+Bz.^2);

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

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;

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
Teperp1=Peper1./rhoe;

file=[dir 'Pe_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper2=fread(fid,'real*8');
fclose(fid);
Peper2=reshape(Peper2,Nx,Ny,Nz);
Teperp2=Peper2./rhoe;

file=[dir 'Pe_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pepar=fread(fid,'real*8');
fclose(fid);
Pepar=reshape(Pepar,Nx,Ny,Nz);
Tepar=Pepar./rhoe;

end

[nx ny nz]= size(Bx)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
%w=1.0./sqrt(Bx(i,:,k).^2+1e-10);

ymax(i,k)=round(sum(w.*y)./sum(w));
%ymax(i,k)=round(sum(y./V(i,:,k).^2)./sum(1./V(i,:,k).^2));
%[dum j] = min(Vx(i,:,k).^2);
%ymax(i,k)=j;
Bxmax(i,k)=Bx(i,round(ymax(i,k)),k);
%Vxmax(i,k)=Vx(i,end/2,k);
Bymax(i,k)=By(i,round(ymax(i,k)),k);
Bzmax(i,k)=Bz(i,round(ymax(i,k)),k);
Exmax(i,k)=Ex(i,round(ymax(i,k)),k);
Eymax(i,k)=Ey(i,round(ymax(i,k)),k);
Ezmax(i,k)=Ez(i,round(ymax(i,k)),k);
Ejmax(i,k)=Ej(i,round(ymax(i,k)),k);
Vxmax(i,k)=Vx(i,round(ymax(i,k)),k);
Vymax(i,k)=Vy(i,round(ymax(i,k)),k); 
Vzmax(i,k)=Vz(i,round(ymax(i,k)),k);
%Pparmax(i,k)=Pxx(i,end/2,k);
end
end

Emax=sqrt(Exmax.^2+Eymax.^2+Ezmax.^2);
Bmax=sqrt(Bxmax.^2+Bymax.^2+Bzmax.^2);
Vmax=sqrt(Vxmax.^2+Vymax.^2+Vzmax.^2);

Aymax=vecpot_uniform(xc,zc,-Vxmax,Vzmax);

for i=0:19
fact = ((i+1)/20).^2;
Ezmax(1+i,:)=Ezmax(1+i,:)*fact;
Ezmax(end-i,:)=Ezmax(end-i,:)*fact;
Ezmax(:,1+i)=Ezmax(:,1+i)*fact;
Ezmax(:,end-i)=Ezmax(:,end-i)*fact;
end

ixcut=180

%tmp=smooth(Ezmax,3);
tmp=Ezmax;
cuttone=[cuttone; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

tmp=smooth(Vxmax,3);
cuttone2=[cuttone2; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),gsmy2z(ymax*dy),['ZGSMmax' ncycle1],[0 0],5,ncycle1,[],1,'x/R_E','y/R_E','Zgsm [RE]');

%load('ymax0_HRmaha3D1.mat')
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),((ymax-ymax0)*dy),['deltaycode' ncycle1],[-1 1],5,ncycle1,[],99,'x/R_E','y/R_E','\deltaZgsm [RE]');

i1=10
i2=Nx-10
j1=10
j2=Nz-10
        
        i1=1
        i2=Nx
        j1=1
        j2=Nz
        

global color_choice symmetric_color initial_time Dt
symmetric_color=0;
color_choice=0;


Vzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vzmax*code_V,['VeYgsm' ncycle1],[-2000 2000],3,ncycle1,[],3,'x/R_E','y/R_E','Vy[km/s]');


Vxsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vxmax*code_V,['VeXgsm' ncycle1],[-2000 2000],3,ncycle1,[],3,'x/R_E','y/R_E','Vx [km/s]');


Vysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vymax*code_V,['VeZgsm' ncycle1],[-2000 2000],3,ncycle1,[],3,'x/R_E','y/R_E','Vz [km/s]');

ox=max(Xgsmrange);
oy=min(Zgsmrange);
oz=min(Ygsmrange);
ddx=dx/Lx*(Xgsmrange(2)-Xgsmrange(1));
ddy=dy/Ly*(Zgsmrange(2)-Zgsmrange(1));
ddz=dz/Lz*(Ygsmrange(2)-Ygsmrange(1));
savevtkvector(-Vxmax,Vymax,  Vzmax, ['Veonmaxp' ncycle1 '.vtk'],'Ve',ddx,ddy,ddz,ox,oy,oz)

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vmax*code_V,['Ve' ncycle1],[0 2000],3,ncycle1,[],3,'x/R_E','y/R_E','Ve [km/s]');



immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Ejmax/B0^2,xc,zc,Vxsm,Vzsm,['ErecELE' ncycle1] ,[-0.2 0.2],0,ncycle1,[],2,'x/R_E','y/R_E','Erec /B_0/V_A');


immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Bymax*code_B,xc,zc,Vxsm,Vzsm,['BZELEgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bx[nT]');


%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Emax*code_E,['E' ncycle1],[0 0],0,ncycle1,[],3,'x/R_E','y/R_E','E[mV/m]');
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bmax*code_B,['B' ncycle1],[0 0],0,ncycle1,[],3,'x/R_E','y/R_E','B[nT]');
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Aymax,['Az' ncycle1],[0 0],0,ncycle1,[],3,'x/R_E','y/R_E','Az[nT]');

close all
%contour(Aymax,30)
%print('-dpng',['Aycontrou' ncycle1])
  

figure(4)
subplot(2,1,1)
plot(mean(Ezmax,2)/B0^2)
subplot(2,1,2)
plot(mean(Bymax,2))
grid on
print('-dpng',['AVGalongcodez' ncycle1]) 
end

return
figure
[nt,ndummy]=size(cuttone);

tc=linspace(0, 154,nt);
pcolor(tc,gsmz2y(yc(j1:j2)),flipud(cuttone(:,j1:j2)')./B0^2)
colorbar
caxis([-.4 .4])
title(num2str(gsmx(xc(ixcut))),'fontsize',18)
xlabel('t (s)','fontsize',18)
ylabel('Y(gsm)','fontsize',18)
set(gca,'ydir','Reverse')
shading interp
                                   set(gca,'fontsize',[18])
set(gcf,'Renderer','zbuffer');
print('-dpng', ['cuttoneEY' num2str(ixcut) '.png'])

figure
[nt,ndummy]=size(cuttone);

tc=linspace(0, 154,nt);
pcolor(tc,gsmz2y(yc(j1:j2)),-flipud(cuttone2(:,j1:j2)')*code_V)
colorbar
title(num2str(gsmx(xc(ixcut))),'fontsize',14)
xlabel('t (s)','fontsize',14)
ylabel('Y(gsm)','fontsize',14)
set(gca,'ydir','Reverse')
shading interp
set(gcf,'Renderer','zbuffer');
print('-dpng', ['cuttoneVx' num2str(ixcut) '.png'])
%!mkdir maxp_frame
%!mv *.png maxp_frame
