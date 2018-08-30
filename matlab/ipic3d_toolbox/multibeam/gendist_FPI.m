load bimodal_ions_FAC.mat           %small box
%load bimodal_ions_FAC_bigbox.mat    %big box
[nx,ny,nz]=size(fcutoff);
vmax=560                             %for small box
%vmax=560*3.75;                       %for big box

dv=2*vmax/(nx-1);
v1D=-vmax:dv:vmax;
[vvx,vvy,vvz]=meshgrid(v1D,v1D,v1D);
f1D=fcutoff(:);
fmin=min(f1D(f1D>0));

%Distribution cannot equal zero for method to work
%so set zeros (if the exist) to small fraction of 
%smallest nonzero value.
f1D(f1D==0)=fmin/1000;      
vx1D=vvx(:)'; vy1D=vvy(:)'; vz1D=vvz(:)';

Np=1e6;          %Number of randomly distributed particles

Ng=size(f1D,1);
ranarr=rand(4,Np);
fcum=cumsum(f1D);
fcum=Ng*fcum/fcum(Ng);
Pg=interp1(fcum',1:Ng,Ng*ranarr(1,:));
Pg=1+floor(Pg);
xp=vx1D(Pg)+dv*ranarr(2,:)-dv/2;
yp=vy1D(Pg)+dv*ranarr(3,:)-dv/2;
zp=vz1D(Pg)+dv*ranarr(4,:)-dv/2;

means=[mean(xp) mean(yp) mean(zp)]
figure(1)
scatter3(xp(1:1000:Np),yp(1:1000:Np),zp(1:1000:Np))
daspect([1 1 1])