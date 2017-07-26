addpath ~/matlab/matlab-parsek
% need to use parsek3D to read it
load cm_new

[nz ny nx nt]=size(Bx);

% Jz for electron should mimic the Az

for it=1:nt
iy=ny/2+1;

time=ig*wci*Dt

bx = squeeze(Bx(:,iy,:,it));
by = squeeze(By(:,iy,:,it));
bz = squeeze(Bz(:,iy,:,it));
ex = squeeze(Ex(:,iy,:,it));
ey = squeeze(Ey(:,iy,:,it));
ez = squeeze(Ez(:,iy,:,it));
jsx0 = squeeze(Jxs0(:,iy,:,it));
jsy0 = squeeze(Jys0(:,iy,:,it));
jsz0 = squeeze(Jzs0(:,iy,:,it));
if(background)
jsxb = squeeze(Jxs2(:,iy,:,it));
jsyb = squeeze(Jys2(:,iy,:,it));
jszb = squeeze(Jzs2(:,iy,:,it));
rhob = squeeze(rhos2(:,iy,:,it));
else
jsxb=0;
jsyb=0;
jszb=0;
rhob=0;
end
rho0 = squeeze(rhos0(:,iy,:,it));


[axis1 axis2]=meshgrid(1:nx,1:nz);
axis1=axis1/nx*Lx;
axis2=axis2/nz*Lz;

label1='x/d_i'
label2='z/d_i'
film='filmXZ/'

ay=vecpot(axis1,axis2,bx,bz);

sec3D_common

end
