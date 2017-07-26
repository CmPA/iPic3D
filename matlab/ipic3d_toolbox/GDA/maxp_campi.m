%maxp_common

cuttone=[];
cuttone2=[];

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


file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vix=fread(fid,'real*8');
fclose(fid);
Vix=reshape(Vix,Nx,Ny,Nz);

%file=[dir 'Ji_y_cycle' ncycle '.gda'];
%fid= fopen(file,'rb');
%Viy=fread(fid,'real*8');
%fclose(fid);
%Viy=reshape(Viy,Nx,Ny,Nz);

%file=[dir 'Ji_z_cycle' ncycle '.gda'];
%fid= fopen(file,'rb');
%Viz=fread(fid,'real*8');
%fclose(fid);
%Viz=reshape(Viz,Nx,Ny,Nz);


end


[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

[Jx, Jy, Jz]=curl(X,Y,Z,permute(Bx,[2 1 3]),permute(By,[2 1 3]),permute(Bz,[2 1 3])); % Valid only because dx=dy=dz
Jx=permute(Jx,[2 1 3])/4/pi;
Jy=permute(Jy,[2 1 3])/4/pi;
Jz=permute(Jz,[2 1 3])/4/pi;

Vx=Jx-Vix;
%Vy=Jy-Viy;
%Vz=Jz-Viz;

[nx ny nz]= size(Bx)
y=1:ny;
for i=1:nx
for k=1:nz
%w=Pper1(i,:,k);
w=1.0./sqrt(Bx(i,:,k).^2+1e-10);

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
Vxmax(i,k)=Vx(i,round(ymax(i,k)),k);
%Vymax(i,k)=Vy(i,round(ymax(i,k)),k); 
%Vzmax(i,k)=Vz(i,round(ymax(i,k)),k);
Vixmax(i,k)=Vix(i,round(ymax(i,k)),k);
%Viymax(i,k)=Viy(i,round(ymax(i,k)),k);
%Vizmax(i,k)=Viz(i,round(ymax(i,k)),k);
%Pparmax(i,k)=Pxx(i,end/2,k);
end
end

Emax=sqrt(Exmax.^2+Eymax.^2+Ezmax.^2);
Bmax=sqrt(Bxmax.^2+Bymax.^2+Bzmax.^2);
%Vmax=sqrt(Vxmax.^2+Vymax.^2+Vzmax.^2);
%Vimax=sqrt(Vixmax.^2+Viymax.^2+Vizmax.^2);

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

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),gsmy2z(ymax*dy),['ZGSMmax' ncycle1],[0 0],5,ncycle1,[],1,'x/R_E','y/R_E','Zgsm [RE]');
colorbar
title(['time (s) = ' ntime])

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
        


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Ezmax/B0^2,['EYgsm' ncycle1] ,[-0.1 0.1],0,ncycle1,[],2,'x/R_E','y/R_E','EYgsm /B_0/V_A');



immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bymax*code_B,['BZgsm' ncycle1],[-10 10],0,ncycle1,[],3,'x/R_E','y/R_E','Bx[nT]');


immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Emax*code_E,['E' ncycle1],[0 0],0,ncycle1,[],3,'x/R_E','y/R_E','E[mV/m]');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bmax*code_B,['B' ncycle1],[0 0],0,ncycle1,[],3,'x/R_E','y/R_E','B[nT]');

close all
  
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vzmax*code_J,['JeYgsm' ncycle1],[-2000 2000]*0,6,ncycle1,[],3,'x/R_E','y/R_E','Jy');

        
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vxmax*code_J,['JeXgsm' ncycle1],[-2000 2000]*0,6,ncycle1,[],3,'x/R_E','y/R_E','Jx ');

 
%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vymax*code_J,['JeZgsm' ncycle1],[-2000 2000]*0,6,ncycle1,[],3,'x/R_E','y/R_E','Jz ');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vmax*code_J,['Je' ncycle1],[0 2000]*0,6,ncycle1,[],3,'x/R_E','y/R_E','Je ');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vizmax*code_J,['JiYgsm' ncycle1],[-1 1]*3,6,ncycle1,[],3,'x/R_E','y/R_E','Jy');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Vixmax*code_J,['JiXgsm' ncycle1],[-1 1]*3,6,ncycle1,[],3,'x/R_E','y/R_E','Jx ');


%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Viymax*code_J,['JiZgsm' ncycle1],[-1 1]*3,6,ncycle1,[],3,'x/R_E','y/R_E','Jz ');

%immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Vimax*code_J,['Ji' ncycle1],[0 1]*3,6,ncycle1,[],3,'x/R_E','y/R_E','Je ');


end
