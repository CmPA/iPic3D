

%nome='psy400.iPic_mp_box_ttcor.0015.dat'; %Jean7 a


%nome='psy400.iPic_mp_box_ttcor.0015.dat'; %Jean13
%nome='psy501.iPic_mp_gbox.0020.dat'; %Jean17
nome='~/Dropbox/Science/ucla/7feb09/feb0709iPICBox.035800UT.dat'; %7feb09
nome='~/Documents/storage/ucla/ucla/HRmaha3D3/feb1508iPIC.034800UT.dat'; %HRmaha3D3
%nome='~/psy501.mp_EC_3DT01_box.0020.dat'

[code_n, code_J, code_V, code_T, code_E, code_B, momentum_corrector] =   code_units(25,1);
e= 1.6022e-19;

Tratio=[1/5,1];
qom=[-256,1];
concentration=[-1.0,1.0];

readdo=1;
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
x =  fscanf(fid,'%f',4);
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

n= n/ code_n/4/pi;

Lx=xmax-xmin
Ly=ymax-ymin

%Nxpic=129;Nypic=65; Nzpic=2; %2D 7feb09

%Nxpic=400+1;Nypic=160+1; Nzpic=160+1; %HRmaha3D3

Nxpic=144+1; Nypic= 252+1; Nzpic=108+1;

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

Tpic=interpmio(x,y,z,T,Xpic,Ypic,Zpic); 

imagesc(squeeze(sqrt(256*Tratio(1)*Tpic(:,:,1))));axis image;colorbar

npic=interpmio(x,y,z,n,Xpic,Ypic,Zpic);

Bxpic=interpmio(x,y,z,Bx,Xpic,Ypic,Zpic);
Bypic=interpmio(x,y,z,By,Xpic,Ypic,Zpic);
Bzpic=interpmio(x,y,z,Bz,Xpic,Ypic,Zpic);

Expic=interpmio(x,y,z,Ex,Xpic,Ypic,Zpic);
Eypic=interpmio(x,y,z,Ey,Xpic,Ypic,Zpic);
Ezpic=interpmio(x,y,z,Ez,Xpic,Ypic,Zpic);

Jxpic(1:Nxpic,1:Nypic,1:Nzpic,1)=interpmio(x,y,z,n.*Vex,Xpic,Ypic,Zpic);
Jypic(1:Nxpic,1:Nypic,1:Nzpic,1)=interpmio(x,y,z,n.*Vey,Xpic,Ypic,Zpic);
Jzpic(1:Nxpic,1:Nypic,1:Nzpic,1)=interpmio(x,y,z,n.*Vez,Xpic,Ypic,Zpic);

Jxpic(1:Nxpic,1:Nypic,1:Nzpic,2)=interpmio(x,y,z,n.*Vix,Xpic,Ypic,Zpic);
Jypic(1:Nxpic,1:Nypic,1:Nzpic,2)=interpmio(x,y,z,n.*Viy,Xpic,Ypic,Zpic);
Jzpic(1:Nxpic,1:Nypic,1:Nzpic,2)=interpmio(x,y,z,n.*Viz,Xpic,Ypic,Zpic);

ns=2;
!rm Initial-Fields_000000.h5
opath='Initial-Fields_000000.h5'
h5create(opath,'/Step#0/Block/Bx/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Bx/0',Bxpic)

h5create(opath,'/Step#0/Block/By/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/By/0',Bypic)

h5create(opath,'/Step#0/Block/Bz/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Bz/0',Bzpic)

h5create(opath,'/Step#0/Block/Bx_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Bx_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,'/Step#0/Block/By_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/By_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,'/Step#0/Block/Bz_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Bz_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,'/Step#0/Block/Ex/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ex/0',Expic)

h5create(opath,'/Step#0/Block/Ey/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ey/0',Eypic)

h5create(opath,'/Step#0/Block/Ez/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ez/0',Ezpic)


h5create(opath,'/Step#0/Block/Ex_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ex_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,'/Step#0/Block/Ey_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ey_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,'/Step#0/Block/Ez_ext/0',[Nxpic, Nypic, Nzpic]);
h5write(opath,'/Step#0/Block/Ez_ext/0',zeros(Nxpic, Nypic, Nzpic))

h5create(opath,['/Step#0/Block/rho_avg/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/rho_avg/0'],zeros(Nxpic, Nypic, Nzpic))

for is=1:ns
h5create(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],Jxpic(1:Nxpic,1:Nypic,1:Nzpic,is)*concentration(is))

h5create(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],Jypic(1:Nxpic,1:Nypic,1:Nzpic,is)*concentration(is))

h5create(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],Jzpic(1:Nxpic,1:Nypic,1:Nzpic,is)*concentration(is))


h5create(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],npic*concentration(is))

vth2= abs(qom(is)) * Tratio(is) * Tpic; 
v0= Jxpic(1:Nxpic,1:Nypic,1:Nzpic,is)./(npic+1e-10);
h5create(opath,['/Step#0/Block/Pxx_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Pxx_' num2str(is-1) '/0'], ...
    (v0.^2 + vth2).* npic.*concentration(is));

v0 = Jypic(1:Nxpic,1:Nypic,1:Nzpic,is)./npic;
h5create(opath,['/Step#0/Block/Pyy_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Pyy_' num2str(is-1) '/0'], ...
    (v0.^2 + vth2).* npic.*concentration(is));

v0 = Jzpic(1:Nxpic,1:Nypic,1:Nzpic,is)./npic;
h5create(opath,['/Step#0/Block/Pzz_' num2str(is-1) '/0'],[Nxpic, Nypic, Nzpic]);
h5write(opath,['/Step#0/Block/Pzz_' num2str(is-1) '/0'], ...
    (v0.^2 + vth2).* npic.*concentration(is));


end
h5writeatt(opath,'/Step#0','nspec',int32(ns));
