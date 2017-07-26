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

[axis1 axis2]=meshgrid(1:nx,1:nz);
axis1=axis1/nx*Lx;
axis2=axis2/nz*Lz;

label1='x/d_i'
label2='z/d_i'
film='filmXZ/'

ay=vecpot(axis1,axis2,bx,bz);

sec3D_common_B

end
