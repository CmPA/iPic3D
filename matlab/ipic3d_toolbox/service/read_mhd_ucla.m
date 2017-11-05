

%nome='psy400.iPic_mp_box_ttcor.0015.dat'; %Jean7 a


%nome='psy400.iPic_mp_box_ttcor.0015.dat'; %Jean13
%nome='psy501.iPic_mp_gbox.0020.dat'; %Jean17
nome='~/Dropbox/Science/ucla/7feb09/feb0709iPICBox.035800UT.dat'; %Jean12

[code_n, code_J, code_V, code_T, code_E, code_B, momentum_corrector] =   code_units();
e= 1.6022e-19;

Tratio=[1/5,1];
 
readdo=0;
if(readdo)
% Data structure:
%x(RE) y(RE) z(RE) bx(nT) by(nT) bz(nT) vx(km/s) vy(km/s) vz(km/s)
% density(cm-3) pressure(pPa) jx(nA/m2) jy(nA/m2) jz(nA/m2)
fid=fopen(nome);

s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);
xmin=x(1); xmax=x(2); dx=x(3); Nx= x(4);
s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);
ymin=x(1); ymax=x(2); dy=x(3); Ny= x(4);
s=fscanf(fid,'%s',4);
x =  fscanf(fid,'%f',4);;
zmin=x(1); zmax=x(2); dz=x(3); Nz= x(4);



a=fscanf(fid,'%f',[14 inf])';


fclose(fid)
end

x=reshape(a(:,1),Nz,Ny,Nx);xmax=max(x(:));xmin=min(x(:));
y=reshape(a(:,2),Nz,Ny,Nx);ymax=max(y(:));ymin=min(y(:));
z=reshape(a(:,3),Nz,Ny,Nx);zmax=max(z(:));zmin=min(z(:));

n=reshape(a(:,10),Nz,Ny,Nx)*1e6;

Bx=reshape(a(:,4),Nz,Ny,Nx)*1e-9/code_B;
By=reshape(a(:,5),Nz,Ny,Nx)*1e-9/code_B;
Bz=reshape(a(:,6),Nz,Ny,Nx)*1e-9/code_B;

Jx=reshape(a(:,12),Nz,Ny,Nx)*1e-9;
Jy=reshape(a(:,13),Nz,Ny,Nx)*1e-9;
Jz=reshape(a(:,14),Nz,Ny,Nx)*1e-9;

Vx=reshape(a(:,7),Nz,Ny,Nx)*1e3*momentum_corrector;
Vy=reshape(a(:,8),Nz,Ny,Nx)*1e3*momentum_corrector;
Vz=reshape(a(:,9),Nz,Ny,Nx)*1e3*momentum_corrector;

Ex = - (Vy.*Bz - Vz.*Bx)/code_V;
Ey = - (Vz.*Bx - Vx.*Bz)/code_V;
Ez = - (Vx.*By - Vy.*Bx)/code_V;

Vix = Vx / code_V;
Viy = Vy / code_V;
Viz = Vz / code_V;

Vex = (Vx - Jx ./n /e)/ code_V;
Vey = (Vy - Jy ./n /e)/ code_V;
Vez = (Vz - Jz ./n /e)/ code_V;

p=reshape(a(:,11),Nz,Ny,Nx)*1e-12;

T = (p ./ n  ) ./ code_T ./(1.0+Tratio(1)/Tratio(2));

n= n/ code_n;

xpic=linspace(xmin,xmax,Nxpic);
pcolor(squeeze(x(Nz/2,:,:)),squeeze(y(Nz/2,:,:)),squeeze(T(Nz/2,:,:))); colorbar; shading interp


Nxpic=501;Nypic=101; Nzpic=2;


dx=(xmax-xmin)/Nxpic;
dy=(ymax-ymin)/Nypic;
dz=(zmax-zmin)/Nzpic;

xpic=linspace(xmin,xmax,Nxpic);

ypic=linspace(ymin,ymax,Nypic);
if(Nzpic<=2) 
    zpic=[(zmin+zmax)/2 (zmin+zmax)/2] ;
else
    zpic=linspace(zmin,zmax,Nzpic);
end

[Xpic,Ypic,Zpic]=ndgrid(xpic,ypic,zpic);

Tpic=interpmio(x,y,z,T,Xpic,Ypic,Zpic); imagesc(squeeze(Tpic(:,:,1)))


