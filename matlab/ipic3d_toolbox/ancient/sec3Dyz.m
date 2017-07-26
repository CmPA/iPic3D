addpath ~/matlab/matlab-parsek
% need to use parsek3D to read it
load cm_new

[nz ny nx nt]=size(Bx);

% Jz for electron should mimic the Az

for it=1:nt
ix=nx/2;

time=ig*Dt*wci;

bx = squeeze(Bx(:,:,ix,it));
by = squeeze(By(:,:,ix,it));
bz = squeeze(Bz(:,:,ix,it));
ex = squeeze(Ex(:,:,ix,it));
ey = squeeze(Ey(:,:,ix,it));
ez = squeeze(Ez(:,:,ix,it));
jsx0 = squeeze(Jxs0(:,:,ix,it));
jsy0 = squeeze(Jys0(:,:,ix,it));
jsz0 = squeeze(Jzs0(:,:,ix,it));
if(background)
jsxb = squeeze(Jxs2(:,:,ix,it));
jsyb = squeeze(Jys2(:,:,ix,it));
jszb = squeeze(Jzs2(:,:,ix,it));
rhob = squeeze(rhos2(:,:,ix,it));
else
jsxb=0;
jsyb=0;
jszb=0;
rhob=0;
end
rho0 = squeeze(rhos0(:,:,ix,it));


[yy zz]=meshgrid(1:ny,1:nz);
yy=yy/ny*Ly;
zz=zz/nz*Lz;

[axis1 axis2]=meshgrid(1:ny,1:nz);
axis1=axis1/ny*Ly;
axis2=axis2/nz*Lz;

label1='y/d_i'
label2='z/d_i'
film='filmYZ/'

ay=vecpot(axis1,axis2,by,bz);

sec3D_common

end
