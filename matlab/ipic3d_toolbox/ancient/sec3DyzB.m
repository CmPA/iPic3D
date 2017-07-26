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

sec3D_common_B

end
