close all
clear all
addpath(genpath('~/iPic3D-github/matlab/ipic3d_toolbox'))

%In nT:
Bx=-0.4189;
By=1.5791;
Bz=12.7405; 
BnT=sqrt(Bx^2+By^2+Bz^2);

n_ref=6; % particle per cc

load('X3d_FAC.mat')
Np=max(size(X2));

load('X3dHE_FAC.mat')
NpHE=max(size(XHE2));

mrcode=256;


 np=n_ref*1e6; % puts ni in m^-3
 ne=np;
 
% Physics Constants
 mu0=4*pi*1e-7;
 eps0=8.8542*1.e-12;
 cphys=1/sqrt(mu0*eps0);
 k=1.3807e-23;
 e= 1.6022e-19;
 mp=1.6726e-27;
 me = mp / mrcode;
 
 wpp=sqrt(np*e^2/mp/eps0);
 wcp=e*BnT*1e-9/mp;
 

code_B = mp*wpp/e;

X2=X2/250*560*2-560; %from grid indeces to Km/s
X2=X2*1e3; % to m/s
X2=X2/cphys; %to code units

XHE2=XHE2/250*560*2-560; %from grid indeces to Km/s
XHE2=XHE2*1e3; % to m/s
XHE2=XHE2/cphys; %to code units

vthi=max(std(X2));
vthe=vthi*sqrt(mrcode/5);
beta_i=vthi.^2*np*mp/(BnT*1e-9)^2*2*mu0*cphys^2
beta_e=vthe.^2*ne*me/(BnT*1e-9)^2*2*mu0*cphys^2


vthiHE=mean(std(XHE2));
vtheHE=vthiHE*sqrt(mrcode/5);

Lx=200;nxc=128;dx=Lx/(nxc+1);
Ly=200;nyc=128;dy=Ly/(nyc+1);
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


nspec=4;
B0=BnT*1e-9/code_B; 

N0=[1;1;1;1]/4/pi;
qom=[-mrcode;1;-mrcode;1];


sqom=sign(qom);

U0=[0;0;0;0];V0=[0;0;0;0];W0=[.0;-.0;0;0];
UTH=[vthe   vthi vtheHE vthiHE];
VTH=UTH;
WTH=UTH;

Bx=zeros(Nx,Ny,Nz);
By=B0*ones(Nx,Ny,Nz);
Bz=zeros(Nx,Ny,Nz);

Ex=zeros(Nx,Ny,Nz);
Ey=zeros(Nx,Ny,Nz);
Ez=zeros(Nx,Ny,Nz);

rho = ones(Nx,Ny,Nz);

!rm FPI-Initial-Fields_000000.h5
opath='FPI-Initial-Fields_000000.h5'

%
% For lambda damping to work, we need to have the intial fields as external
% fields so they are not damped in the damping region
%

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

for is=1:nspec
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
h5writeatt(opath,'/Step#0','nspec',int32(nspec));



!rm FPI-Initial-Partcl_000000.h5
opath='FPI-Initial-Partcl_000000.h5'



qp=n_ref*Lx*Ly*Lz/Np
for is=1:2
    qp=N0(is)*Lx*Ly/Np*sqom(is)
h5create(opath,['/Step#0/q_' num2str(is-1)],Np);
h5write(opath,['/Step#0/q_' num2str(is-1) ],qp*ones(Np,1)) 

if(is==2)
    fpi=1
    if(fpi)
h5create(opath,['/Step#0/u_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/u_' num2str(is-1) ],X2(:,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/v_' num2str(is-1) ],X2(:,2)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/w_' num2str(is-1) ],X2(:,3))  
    else
h5create(opath,['/Step#0/u_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/u_' num2str(is-1) ],vthi*randn(Np,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/v_' num2str(is-1) ],vthi*randn(Np,1)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/w_' num2str(is-1) ],vthi*randn(Np,1)) 
    end
else
h5create(opath,['/Step#0/u_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/u_' num2str(is-1) ],vthe*randn(Np,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/v_' num2str(is-1) ],vthe*randn(Np,1)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/w_' num2str(is-1) ],vthe*randn(Np,1))   
    
end

h5create(opath,['/Step#0/x_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/x_' num2str(is-1) ],rand(Np,1)*Lx) 
h5create(opath,['/Step#0/y_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/y_' num2str(is-1) ],rand(Np,1)*Ly) 
h5create(opath,['/Step#0/z_' num2str(is-1) ],Np);
h5write(opath,['/Step#0/z_' num2str(is-1) ],rand(Np,1)*Lz)  
h5writeatt(opath,'/Step#0',['npart_' num2str(is-1)] ,int32(Np));
end

for is=3:4
    qp=N0(is)*Lx*Ly/Np*sqom(is)
h5create(opath,['/Step#0/q_' num2str(is-1)],NpHE);
h5write(opath,['/Step#0/q_' num2str(is-1) ],qp*QHE) 
is
qp*QHE

if(is==4)
    fpi=1
    if(fpi)
h5create(opath,['/Step#0/u_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/u_' num2str(is-1) ],XHE2(:,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/v_' num2str(is-1) ],XHE2(:,2)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/w_' num2str(is-1) ],XHE2(:,3))  
    else
h5create(opath,['/Step#0/u_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/u_' num2str(is-1) ],vthiHE*randn(NpHE,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/v_' num2str(is-1) ],vthiHE*randn(NpHE,1)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/w_' num2str(is-1) ],vthiHE*randn(NpHE,1)) 
    end
else
h5create(opath,['/Step#0/u_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/u_' num2str(is-1) ],vtheHE*randn(NpHE,1)) 
h5create(opath,['/Step#0/v_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/v_' num2str(is-1) ],vtheHE*randn(NpHE,1)) 
h5create(opath,['/Step#0/w_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/w_' num2str(is-1) ],vtheHE*randn(NpHE,1))   
    
end

h5create(opath,['/Step#0/x_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/x_' num2str(is-1) ],rand(NpHE,1)*Lx) 
h5create(opath,['/Step#0/y_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/y_' num2str(is-1) ],rand(NpHE,1)*Ly) 
h5create(opath,['/Step#0/z_' num2str(is-1) ],NpHE);
h5write(opath,['/Step#0/z_' num2str(is-1) ],rand(NpHE,1)*Lz)  
h5writeatt(opath,'/Step#0',['npart_' num2str(is-1)] ,int32(NpHE));
end

h5writeatt(opath,'/Step#0','nspec',int32(nspec));
