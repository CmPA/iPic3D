clear all;
close all

Lx=20;nxc=128;dx=Lx/(nxc+1);
Ly=20;nyc=128;dy=Ly/(nyc+1);
Lz=1;nzc=1;dz=Lz/(nzc+1);

x = linspace(0,Lx,nxc+1);
y = linspace(0,Ly,nyc+1);
z = linspace(0,Lz,nzc+1);

[X,Y,Z] = ndgrid(x,y,z);
[Nx,Ny,Nz]=size(X)
% 
% X=permute(X,[3 2 1]);
% Y=permute(Y,[3 2 1]);
% Z=permute(Z,[3 2 1]);


ns=2;
B0=0.0097; L=0.5; 

N0=[1;1]/4/pi;
qom=[-256;1];
sqom=sign(qom);

U0=[0;0];V0=[0;0];W0=[.00325;-.01624];
UTH=[0.045   0.0063];
VTH=UTH;
WTH=UTH;

Bx=B0*tanh((Y-Ly/2)./L);
By=zeros(Nx,Ny,Nz);
Bz=zeros(Nx,Ny,Nz);

Ex=zeros(Nx,Ny,Nz);
Ey=zeros(Nx,Ny,Nz);
Ez=zeros(Nx,Ny,Nz);

rho = (.01+sech((Y-Ly/2)./L).^2);


!rm Harris-Initial-Fields_000000.h5
opath='Harris-Initial-Fields_000000.h5'

h5create(opath,'/Step#0/Block/Bx/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bx/0',Bx)

h5create(opath,'/Step#0/Block/By/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/By/0',By)

h5create(opath,'/Step#0/Block/Bz/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bz/0',Bz) 

h5create(opath,'/Step#0/Block/Bx_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bx_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,'/Step#0/Block/By_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/By_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,'/Step#0/Block/Bz_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Bz_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,'/Step#0/Block/Ex/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ex/0',Ex)

h5create(opath,'/Step#0/Block/Ey/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ey/0',Ey)

h5create(opath,'/Step#0/Block/Ez/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ez/0',Ez)

h5create(opath,'/Step#0/Block/Ex_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ex_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,'/Step#0/Block/Ey_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ey_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,'/Step#0/Block/Ez_ext/0',[Nx, Ny, Nz]);
h5write(opath,'/Step#0/Block/Ez_ext/0',zeros(Nx,Ny,Nz))

h5create(opath,['/Step#0/Block/rho_avg/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/rho_avg/0'],zeros(Nx, Ny, Nz))

for is=1:ns
h5create(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],N0(is)*rho*U0(is)*sqom(is))

h5create(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],N0(is)*rho*V0(is)*sqom(is))

h5create(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],N0(is)*rho*W0(is)*sqom(is))

h5create(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],N0(is)*rho*sqom(is))

h5create(opath,['/Step#0/Block/Pxx_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Pxx_' num2str(is-1) '/0'],N0(is)*rho*(U0(is).^2+UTH(is).^2)*sqom(is))

h5create(opath,['/Step#0/Block/Pyy_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Pyy_' num2str(is-1) '/0'],N0(is)*rho*(V0(is).^2+VTH(is).^2)*sqom(is))

h5create(opath,['/Step#0/Block/Pzz_' num2str(is-1) '/0'],[Nx, Ny, Nz]);
h5write(opath,['/Step#0/Block/Pzz_' num2str(is-1) '/0'],N0(is)*rho*(W0(is).^2+WTH(is).^2)*sqom(is))
end
h5writeatt(opath,'/Step#0','nspec',int32(ns));