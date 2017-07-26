%maxp_common

for cycle=Ncyc_ini:1000:Ncyc_max

time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
ntime=num2str(time);

qom=-256;

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)



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

B2=sqrt(Bx.^2+By.^2+Bz.^2);



[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

[Jx, Jy, Jz]=curl(X,Y,Z,permute(Bx,[2 1 3]),permute(By,[2 1 3]),permute(Bz,[2 1 3])); % Valid only because dx=dy=dz
Jx=permute(Jx,[2 1 3])/4/pi;
Jy=permute(Jy,[2 1 3])/4/pi;
Jz=permute(Jz,[2 1 3])/4/pi;


file=[dir 'Pi_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);


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

Vx=Jx-Vix;
Vy=Jy-Viy;
Vz=Jz-Viz;

Vx=Vx./rhoe;
Vy=Vy./rhoe;
Vz=Vz./rhoe;


Vix=Vix./rhoi;
Viy=Viy./rhoi;
Viz=Viz./rhoi;

%Jx=Vex.*rhoe+Vix.*rhoi;
%Jy=Vey.*rhoe+Viy.*rhoi;
%Jz=Vez.*rhoe+Viz.*rhoi;


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


%
% Compute divergence of pressure tensor
%

P=(PXX+PYY+PZZ)/3;

[FX]=divergence(X,Y,Z,permute(PXX,[2 1 3]),permute(PXY,[2 1 3]),permute(PXZ,[2 1 3])); % Valid only because dx=dy=dz
FX=permute(FX,[2 1 3])./rhoe;

[FY]=divergence(X,Y,Z,permute(PXY,[2 1 3]),permute(PYY,[2 1 3]),permute(PYZ,[2 1 3])); % Valid only because dx=dy=dz
FY=permute(FY,[2 1 3])./rhoe;

[FZ]=divergence(X,Y,Z,permute(PXZ,[2 1 3]),permute(PYZ,[2 1 3]),permute(PZZ,[2 1 3])); % Valid only because dx=dy=dz
FZ=permute(FZ,[2 1 3])./rhoe;

%%% Attention, the -1/ne in the defintiin is implemented as 1/rhoe becasue rhoe is negative.
%
% Compute Ve x B
%

[Zx, Zy, Zz] = cross_prod(-Vx, -Vy, -Vz, Bx, By, Bz);

[Hx, Hy, Hz] = cross_prod(Jx, Jy, Jz, Bx, By, Bz);

Hx=Hx./abs(rhoe);
Hy=Hy./abs(rhoe);
Hz=Hz./abs(rhoe);

end

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);

w=1.0./(B2(i,:,k).^2+1e-10);

%w=(Jx(i,:,k).^2+Jy(i,:,k).^2+Jz(i,:,k).^2)

ymax(i,k)=round(sum(w.*y)./sum(w));

%[dummy ymax(i,k)]=max(w);
nav=20;
jrange=max(1,round(ymax(i,k))-nav):min(round(ymax(i,k))+nav,ny);


Fx_max(i,k)=mean(FX(i,jrange,k),2);
Fy_max(i,k)=mean(FY(i,jrange,k),2);
Fz_max(i,k)=mean(FZ(i,jrange,k),2);

Zx_max(i,k)=mean(Zx(i,jrange,k),2);
Zy_max(i,k)=mean(Zy(i,jrange,k),2);
Zz_max(i,k)=mean(Zz(i,jrange,k),2);

Hx_max(i,k)=mean(Hx(i,jrange,k),2);
Hy_max(i,k)=mean(Hy(i,jrange,k),2);
Hz_max(i,k)=mean(Hz(i,jrange,k),2);

Ex_max(i,k)=mean(Ex(i,jrange,k),2);
Ey_max(i,k)=mean(Ey(i,jrange,k),2);
Ez_max(i,k)=mean(Ez(i,jrange,k),2);

Bx_max(i,k)=mean(Bx(i,jrange,k),2);
By_max(i,k)=mean(By(i,jrange,k),2);
Bz_max(i,k)=mean(Bz(i,jrange,k),2);

end
end

% For 1mar08C
%Xgsm=-10;
% for 15feb08 (HRmaha3D3)
Xgsm=-28;

x= -Lx*(Xgsm-Xgsmrange(2))/(Xgsmrange(2)-Xgsmrange(1))
ix=round(x/dx)

plot_maxp=0;
if (plot_maxp)
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Zx_max/B0^2,['Zx' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[Vi \times B]x(gsm) /B_0/V_A');
Zysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zy_max/B0^2,['Zz' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[Vi \times B]z(gsm) /B_0/V_A');
Zzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Zz_max/B0^2,['Zy' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[Vi \times B]y(gsm) /B_0/V_A');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Hx_max/B0^2,['Hx' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','[J \times B]x(gsm) /B_0/V_A/\rho');
Hysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Hy_max/B0^2,['Hz' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','[J \times B]z(gsm) /B_0/V_A/\rho');
Hzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Hz_max/B0^2,['Hy' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','[J \times B]y(gsm) /B_0/V_A/\rho');

Bxsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Bx_max*code_B,['Bx' ncycle1],[-0.1 0.1]*0,3,ncycle1,[],3,'x/R_E','y/R_E','Bx(gsm) ');
Bysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),By_max*code_B,['Bz' ncycle1],[-0.1 0.1]*0,3,ncycle1,[],3,'x/R_E','y/R_E','Bz(gsm) ');
Bzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Bz_max*code_B,['By' ncycle1],[-0.1 0.1]*0,3,ncycle1,[],3,'x/R_E','y/R_E','By(gsm) ');


Exsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Ex_max/B0^2,['Ex' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ex(gsm) /B_0/V_A');
Eysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Ey_max/B0^2,['Ez' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ez(gsm) /B_0/V_A');
Ezsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Ez_max/B0^2,['Ey' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','Ey(gsm) /B_0/V_A');

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-Fx_max/B0^2,['Fx' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[\nabla \cdot P]x(gsm) /B_0/V_A/\rho');
Fysm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Fy_max/B0^2,['Fz' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[\nabla \cdot P]z(gsm) /B_0/V_A/\rho');
Fzsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Fz_max/B0^2,['Fy' ncycle1],[-0.1 0.1]*4,3,ncycle1,[],3,'x/R_E','y/R_E','-[\nabla \cdot P]y(gsm) /B_0/V_A/\rho');


Bx=Bx./B2;
By=By./B2;
Bz=Bz./B2;

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
w=1.0./sqrt(B2(i,:,k).^2+1e-10);

ymax(i,k)=round(sum(w.*y)./sum(w));
%[dummy ymax(i,k)]=max(w);

iy=round(ymax(i,k))-5:round(ymax(i,k))+5;
p(1,1)=mean(PXX(i,iy,k));
p(1,2)=mean(PXY(i,iy,k));
p(1,3)=mean(PXZ(i,iy,k));
p(2,2)=mean(PYY(i,iy,k));
p(2,3)=mean(PYZ(i,iy,k));
p(3,3)=mean(PZZ(i,iy,k));
p(2,1)=p(1,2);
p(3,1)=p(1,3);
p(3,2)=p(2,3);

b(1)=mean(Bx(i,iy,k));
b(2)=mean(By(i,iy,k));
b(3)=mean(Bz(i,iy,k));
b=reshape(b,3,1);
b=b./sqrt(sum(b.^2));

%%%%%%%%%%
% Scudder
%%%%%%%%%%

for l=1:3
N1(l,:)=cross(b,p(l,:));
end

for l=1:3
N(:,l)=cross(b,N1(:,l));
end

lamb=sort(eig(N));
lambda1(i,k)=lamb(1);
lambda2(i,k)=lamb(2);
lambda3(i,k)=lamb(3);

Agyro(i,k)= 2*(lamb(3)-lamb(2))/(lamb(3)+lamb(2));


%%%%%%%%%%%
% Aunai
%%%%%%%%%%%

Tr=(p(1,1)+p(2,2)+p(3,3));
Ppar=b'*p*b;
Pper=(Tr-Ppar)/2;
G=eye(3,3)*Pper+(Ppar-Pper)*kron(b,b');
N=p-G;
Agyro_aunai(i,k)=sqrt(sum(N(:).^2))./Tr;

%%%%%%%%%%%
% Swisdak
%%%%%%%%%%%

I2=p(1,1)*p(2,2)+p(1,1)*p(3,3)+p(2,2)*p(3,3);
I2=I2-(p(1,2).^2+p(1,3).^2+p(2,3).^2);
Q=1-4*I2./((Tr-Ppar).*(Tr+3*Ppar));
Nongyro_swisdak(i,k)=sqrt(Q);
% The following formula form Swisdak paper is actually wrong
%Nongyro_aunai(i,k)=sqrt(8*(p(1,2).^2+p(1,3).^2+p(2,3).^2))./(Ppar+2*Pper);
                                 
end
end

Agsm=immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Agyro,['Agyro' ncycle1],[0 .4],3,ncycle1,[],3,'x/R_E','y/R_E','Agyro');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Agyro_aunai,['Agyro_Aunai' ncycle1],[0 .2],3,ncycle1,[],3,'x/R_E','y/R_E','Agyro');
immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),Nongyro_swisdak,['Nongyro_Swisdak' ncycle1],[0 .2],3,ncycle1,[],3,'x/R_E','y/R_E','NonGyro');

close all

end

limit_plot=[Ygsmrange(1)+1 Ygsmrange(2)-1]
                                 
figure(1)
%fact=(Agsm(ix,:)-mean(Agsm(ix,:)))./max(Agsm(ix,:)-mean(Agsm(ix,:)));
fact=1;

if(plot_maxp)
plot(gsmz2y(zc),Eysm(ix,:)./fact,'k',gsmz2y(zc),Zysm(ix,:)./fact,'b',gsmz2y(zc),Fysm(ix,:).*fact,'r',gsmz2y(zc),Hysm(ix,:).*fact,'m',gsmz2y(zc),Agsm(ix,:),'g')
legend('Ez','[Vi\times B]z','[\nabla \cdot P]z','Hall','location','SouthWest')
grid on
xlim(limit_plot)
print -dpng flythroughZgsm.png 
figure(1)
%plot(gsmz2y(zc),Ezsm(ix,:).*fact,'k',gsmz2y(zc),Zzsm(ix,:).*fact,'b',gsmz2y(zc),Fzsm(ix,:).*fact,'r',gsmz2y(zc),Hzsm(ix,:).*fact,'m',gsmz2y(zc),Agsm(ix,:),'g')
end
                                 
%ixr=ix-5:ix+5;
ixr=ix;
plot(gsmz2y(zc),code_E*mean(Ex_max(ixr,:),1),'k',gsmz2y(zc),code_E*mean(Ey_max(ixr,:),1),'r',gsmz2y(zc),code_E*mean(Ez_max(ixr,:),1),'g')
legend('Ex','Ez','Ey','location','SouthWest')
grid on 
pbaspect([5 1 1])
xlim(limit_plot)
print -dpng flythroughYEgsm.png
if(plot_maxp)
plot(gsmz2y(zc),mean(Bxsm(ixr,:),1),'k',gsmz2y(zc),mean(Bysm(ixr,:),1),'r',gsmz2y(zc),mean(Bzsm(ixr,:),1),'g')
legend('Bx','Bz','By','location','SouthWest')
grid on 
pbaspect([5 1 1])
xlim(limit_plot)
print -dpng flythroughYBgsm.png 
plot(gsmz2y(zc),mean(Agsm(ixr,:),1),'g')
grid on
pbaspect([5 1 1])
xlim(limit_plot)
print -dpng flythroughAgsm.png 
plot(gsmz2y(zc),mean(Fzsm(ixr,:),1),'r',gsmz2y(zc),mean(Ezsm(ixr,:)-Zzsm(ixr,:),1),'y')
legend('[\nabla \cdot P]y','[E+Ve\times B]y','location','SouthWest')
grid on
pbaspect([5 1 1])
xlim(limit_plot)
print -dpng flythroughYOHMgsm.png 
plot(gsmz2y(zc),mean(Fysm(ixr,:),1),'r',gsmz2y(zc),mean(Eysm(ixr,:)-Zysm(ixr,:),1),'y')
legend('[\nabla \cdot P]z','[E+Ve\times B]z','location','SouthWest')
grid on
pbaspect([5 1 1])
xlim(limit_plot)
print -dpng flythroughZOHMgsm.png
end
                                 
end
        

