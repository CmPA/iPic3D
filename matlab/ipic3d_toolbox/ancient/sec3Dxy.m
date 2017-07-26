addpath ~/matlab/matlab-parsek
% need to use parsek3D to read it

[nz ny nx nt]=size(Bx);


% Jz for electron should mimic the Az

for it=1:nt
iz=round(nz/2);

bx = squeeze(Bx(iz,:,:,it));
by = squeeze(By(iz,:,:,it));
bz = squeeze(Bz(iz,:,:,it));
ex = squeeze(Ex(iz,:,:,it));
ey = squeeze(Ey(iz,:,:,it));
ez = squeeze(Ez(iz,:,:,it));
jsx0 = squeeze(Jxs0(iz,:,:,it));
jsy0 = squeeze(Jys0(iz,:,:,it));
jsz0 = squeeze(Jzs0(iz,:,:,it));
if(background)
jsxb = squeeze(Jxs2(iz,:,:,it));
jsyb = squeeze(Jys2(iz,:,:,it));
jszb = squeeze(Jzs2(iz,:,:,it));
rhob = squeeze(rhos2(iz,:,:,it));
else
jsxb=0;
jsyb=0;
jszb=0;
rhob=0;
end
rho0 = squeeze(rhos0(iz,:,:,it));

[axis1 axis2]=meshgrid(1:nx,1:ny);
axis1=axis1/nx*Lx;
axis2=axis2/ny*Ly;

label1='x/d_i'
label2='y/d_i'
film='filmXY/'

ay=vecpot(axis1,axis2,bx,by);

sec3D_common

end
