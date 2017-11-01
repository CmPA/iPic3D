Lx=1;dx=0.2;
Ly=1;dy=0.1;
Lz=1;dz=0.5;

x = 0:dx:Lx;
y = 0:dy:Ly;
z = 0:dz:Lz;
[X,Y,Z] = ndgrid(x,y,z);
[Nx,Ny,Nz]=size(X)

X=permute(X,[3 2 1]);
Y=permute(Y,[3 2 1]);
Z=permute(Z,[3 2 1]);


ns=2;
B0=1; L=0.1; 

N0=[1;1];
qom=[-256;1];

U0=[0;0];V0=[0;0];W0=[1;-1];


Bx=B0*tanh(Y./L);
By=zeros(Nz,Ny,Nx);
Bz=zeros(Nz,Ny,Nx);

Ex=zeros(Nz,Ny,Nx);
Ey=zeros(Nz,Ny,Nx);
Ez=zeros(Nz,Ny,Nx);

rho = sech(Y./L).^2;


!rm Initial-Fields_000000.h5
opath='Initial-Fields_000000.h5'
h5create(opath,'/Step#0/Block/Bx/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/Bx/0',Bx)

h5create(opath,'/Step#0/Block/By/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/By/0',By)

h5create(opath,'/Step#0/Block/Bz/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/Bz/0',Bz)

h5create(opath,'/Step#0/Block/Ex/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/Ex/0',Ex)

h5create(opath,'/Step#0/Block/Ey/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/Ey/0',Ey)

h5create(opath,'/Step#0/Block/Ez/0',[Nz, Ny, Nx]);
h5write(opath,'/Step#0/Block/Ez/0',Ez)

for is=1:ns
h5create(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],[Nz, Ny, Nx]);
h5write(opath,['/Step#0/Block/Jx_' num2str(is-1) '/0'],N0(is)*rho*U0(is)*qom(is))

h5create(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],[Nz, Ny, Nx]);
h5write(opath,['/Step#0/Block/Jy_' num2str(is-1) '/0'],N0(is)*rho*V0(is)*qom(is))

h5create(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],[Nz, Ny, Nx]);
h5write(opath,['/Step#0/Block/Jz_' num2str(is-1) '/0'],N0(is)*rho*W0(is)*qom(is))

h5create(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],[Nz, Ny, Nx]);
h5write(opath,['/Step#0/Block/rho_' num2str(is-1) '/0'],N0(is)*rho*qom(is))
end
h5writeatt(opath,'/Step#0','nspec',int32(ns));