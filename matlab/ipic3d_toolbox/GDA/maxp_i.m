%maxp_common

cuttone=[];
cuttone2=[];
sqrt(2)
for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

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

B2=(Bx.*Bx+By.*By+Bz.*Bz);


[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

%[Jx, Jy, Jz]=curl(X,Y,Z,permute(Bx,[2 1 3]),permute(By,[2 1 3]),permute(Bz,[2 1 3])); % Valid only because dx=dy=dz
%Jx=permute(Jx,[2 1 3])/4/pi;
%Jy=permute(Jy,[2 1 3])/4/pi;
%Jz=permute(Jz,[2 1 3])/4/pi;

%J2=(Jx.*Jx+Jy.*Jy+Jz.*Jz);


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

Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./sqrt(B2);
%Ej=(Ex.*Jx+Ey.*Jy+Ez.*Jz)./sqrt(J2);
Ej=(Ex.*Bz-Ez.*Bx)./sqrt(Bx.^2+Bz.^2);

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

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;


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
Ezmax(i,k)=Ez(i,round(ymax(i,k)),k);
Exmax(i,k)=Ex(i,round(ymax(i,k)),k);
Ejmax(i,k)=Ej(i,round(ymax(i,k)),k);
Vxmax(i,k)=Vx(i,round(ymax(i,k)),k);
Vymax(i,k)=Vy(i,round(ymax(i,k)),k); %14 is oine RE above
Vzmax(i,k)=Vz(i,round(ymax(i,k)),k);
%Pparmax(i,k)=Pxx(i,end/2,k);
end
end

for i=0:19
fact = ((i+1)/20).^2;
Ezmax(1+i,:)=Ezmax(1+i,:)*fact;
Ezmax(end-i,:)=Ezmax(end-i,:)*fact;
Ezmax(:,1+i)=Ezmax(:,1+i)*fact;
Ezmax(:,end-i)=Ezmax(:,end-i)*fact;
end


ixcut=180

tmp=smooth(Ezmax,3);
cuttone=[cuttone; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

tmp=smooth(Vxmax,3);
cuttone2=[cuttone2; mean(tmp(round(ixcut-10:ixcut+10),:),1)];

global color_choice symmetric_color initial_time Dt
symmetric_color=0;
color_choice=0;

%subplot(3,1,1)
%surf(ymax)
figure(1)
%pcolor(ymax')
%shading interp
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),gsmy2z(ymax*dy),['ZGSMmax' ncycle1],[-3.2 -2.7],5,ncycle1,[],1,'x/R_E','y/R_E','Zgsm [RE]')
%colorbar
%title(['time (s) = ' ntime])

%load('ymax0_HRmaha3D1.mat')
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),((ymax-ymax0)*dy),['deltaycode' ncycle1],[-1 1],5,ncycle1,[],99,'x/R_E','y/R_E','\deltaZgsm [RE]')

i1=10
i2=Nx-10
j1=10
j2=Nz-10
        
        
        i1=1
        i2=Nx
        j1=1
        j2=Nz
        
%subplot(3,1,2)
%pcolor(Ezmax')

%immagine_xy(gsmx([xc(i1) xc(i2)]),gsmz2y([zc(j1) zc(j2)]),Ezmax(i1:i2,j1:j2)/B0^2,['EYgsm' ncycle1] ,[-0.4 0.4],0,ncycle1,[],2,'x/R_E','y/R_E','EYgsm /B_0/V_A')
%colorbar

Vzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vzmax*code_V,['ViYgsm' ncycle1],[-600 600],3,ncycle1,[],1,'x/R_E','y/R_E','Vy[km/s]');

Vxsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vxmax*code_V,['ViXgsm' ncycle1],[-600 600],3,ncycle1,[],1,'x/R_E','y/R_E','Vx [km/s]');

        ox=max(Xgsmrange);
        oy=min(Zgsmrange);
        oz=min(Ygsmrange);
        ddx=dx/Lx*(Xgsmrange(2)-Xgsmrange(1));
        ddy=dy/Ly*(Zgsmrange(2)-Zgsmrange(1));
        ddz=dz/Lz*(Ygsmrange(2)-Ygsmrange(1));
savevtkvector(-Vxmax,  Vymax, Vzmax, ['Vionmaxp' ncycle1 '.vtk'],'Vi',ddx,ddy,ddz,ox,oy,oz)
savevtkvector(-Bxmax,  Bymax, Bzmax, ['Bonmaxp' ncycle1 '.vtk'],'B',ddx,ddy,ddz,ox,oy,oz)


dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Ejmax/B0^2,xc,zc,Vxsm,Vzsm,['Erec' ncycle1] ,[-0.2 0.2],3,ncycle1,[],2,'x/R_E','y/R_E','Erec /B_0/V_A');


dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Ezmax/B0^2,xc,zc,Vxsm,Vzsm,['EYgsm' ncycle1] ,[-0.2 0.2],3,ncycle1,[],2,'x/R_E','y/R_E','EYgsm /B_0/V_A');

dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Exmax/B0^2,xc,zc,Vxsm,Vzsm,['Exgsm' ncycle1] ,[-0.2 0.2],3,ncycle1,[],2,'x/R_E','y/R_E','EXgsm /B_0/V_A');

dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Bymax*code_B,xc,zc,Vxsm,Vzsm,['BZgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bz[nT]');

dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),-Bxmax*code_B,xc,zc,Vxsm,Vzsm,['BXgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bz[nT]');
dummy=immagine_xy_vel(gsmx([0 Lx]),gsmz2y([0 Lz]),Bzmax*code_B,xc,zc,Vxsm,Vzsm,['BYgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bz[nT]');
%caxis([-4 4]*1e-4)
%shading interp
%set(gcf,'Renderer','zbuffer');
%print -dpng caz.png
        overplot=0;
        
        if(overplot)
%
% add velocity streamfunction
%
        
hlines=streamslice(gsmx(xc),gsmz2y(zc),-(Vxmax'),(Vzmax'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
name=['Bn_XYIONcombo' ncycle1];
print('-dpng','-r300',[name '.png'])

        close all
hold on
        skip=10
        quiver(gsmx(xc(1:skip:end)),gsmz2y(zc(1:skip:end)),-(Vxmax(1:skip:end,1:skip:end)'),(Vzmax(1:skip:end,1:skip:end)'),'color',[1 1 1]);
        name=['Bn_XYquiver' ncycle1];
        print('-dpng','-r300',[name '.png'])
return
end
        
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vzmax*code_V,['ViYgsm' ncycle1],[-600 600],0,ncycle1,[],1,'x/R_E','y/R_E','Vy[km/s]')
%colorbar
        
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vxmax*code_V,['ViXgsm' ncycle1],[-600 600],0,ncycle1,[],1,'x/R_E','y/R_E','Vx [km/s]')

        colorbar
 
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vymax*code_V,['ViZgsm' ncycle1],[-600 600],0,ncycle1,[],1,'x/R_E','y/R_E','Vz [km/s]')
%colorbar
  
%immagine_xy(gsmx([xc(i1) xc(i2)]),gsmz2y([zc(j1) zc(j2)]),Vzmax*code_V,['VYgsm' ncycle1],[-600 600],0,ncycle1,[],3,'x/R_E','y/R_E','Vy[km/s]')
%colorbar
        
%immagine_xy(gsmx([xc(i1) xc(i2)]),gsmz2y([zc(j1) zc(j2)]),-Vxmax*code_V,['VXgsm' ncycle1],[-500 500],3,ncycle1,[],3,'x/R_E','y/R_E','Vx [km/s]')
%colorbar
 
immagine_xy(gsmx([xc(i1) xc(i2)]),gsmz2y([zc(j1) zc(j2)]),Vymax*code_V,['ViZgsm' ncycle1],[-600 600],0,ncycle1,[],3,'x/R_E','y/R_E','Vz [km/s]')
colorbar
      
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
title(num2str(gsmx(xc(ixcut))),'fontsize',14)
xlabel('t (s)','fontsize',14)
ylabel('Y(gsm)','fontsize',14)
set(gca,'ydir','Reverse')
shading interp
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
