addpath ~/matlab/matlab-parsek
% need to use parsek3D to read it

[nz ny nx nt]=size(Bx);


% Jz for electron should mimic the Az

for it=1:nt
iz=round(nz/2);

bx = squeeze(Bx(iz,:,:,it));
by = squeeze(By(iz,:,:,it));
bz = squeeze(Bz(iz,:,:,it));

[axis1 axis2]=meshgrid(1:nx,1:ny);
axis1=axis1/nx*Lx;
axis2=axis2/ny*Ly;

label1='x/d_i'
label2='y/d_i'
film='filmXY/'

ay=vecpot(axis1,axis2,bx,by);

sec3D_common_B

end
