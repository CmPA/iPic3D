close all
addpath(genpath('../../ipic3d_toolbox')); % Point to the directory where the iPic3D toolbox is

%clear all
leggo=2; 

for cycle=70010

if(leggo==2)   
ncycle = num2str(cycle,'%06d');
dir = '~/Dropbox/Science/san_diego/energy-test/'
Lx=25;Ly=50;Lz=1;
qom=[-40,2.5];
    namefile = 'TwoCoils2D-Fields';
    fn=[dir,namefile,'_',ncycle,'.h5'];


    Bx = hdf5read(fn,'/Step#0/Block/Bx/0/')+hdf5read(fn,'/Step#0/Block/Bx_ext/0/');
    [Nx,Ny,Nz]=size(Bx)
    dx=Lx/(Nx-1);
    dy=Ly/(Ny-1);
    dz=Lz/(Nz-1);

    By = hdf5read(fn,'/Step#0/Block/By/0/')+hdf5read(fn,'/Step#0/Block/By_ext/0/');
    Bz = hdf5read(fn,'/Step#0/Block/Bz/0/')+hdf5read(fn,'/Step#0/Block/Bz_ext/0/');
    
    
    Ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    Ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    Ez = hdf5read(fn,'/Step#0/Block/Ez/0/');
    Jex = hdf5read(fn,'/Step#0/Block/Jx_0/0/');%+hdf5read(fn,'/Step#0/Block/Jx_2/0/');
    Jey = hdf5read(fn,'/Step#0/Block/Jy_0/0/');%+hdf5read(fn,'/Step#0/Block/Jy_2/0/');
    Jez = hdf5read(fn,'/Step#0/Block/Jz_0/0/');%+hdf5read(fn,'/Step#0/Block/Jz_2/0/');
    Jix = hdf5read(fn,'/Step#0/Block/Jx_1/0/');%+hdf5read(fn,'/Step#0/Block/Jx_3/0/');
    Jiy = hdf5read(fn,'/Step#0/Block/Jy_1/0/');%+hdf5read(fn,'/Step#0/Block/Jy_3/0/');
    Jiz = hdf5read(fn,'/Step#0/Block/Jz_1/0/');%+hdf5read(fn,'/Step#0/Block/Jz_3/0/');
    
    rhoe = hdf5read(fn,'/Step#0/Block/rho_0/0/');%+hdf5read(fn,'/Step#0/Block/rho_2/0/');
    rhoi = hdf5read(fn,'/Step#0/Block/rho_1/0/');%+hdf5read(fn,'/Step#0/Block/rho_3/0/');

    Pexx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxx_2/0/');
    Peyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/');%+hdf5read(fn,'/Step#0/Block/Pyy_2/0/');
    Pezz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pzz_2/0/');
    Pexy = hdf5read(fn,'/Step#0/Block/Pxy_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxy_2/0/');    
    Pexz = hdf5read(fn,'/Step#0/Block/Pxz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pxz_2/0/');
    Peyz = hdf5read(fn,'/Step#0/Block/Pyz_0/0/');%+hdf5read(fn,'/Step#0/Block/Pyz_2/0/');
    
    Pixx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxx_3/0/');
    Piyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/');%+hdf5read(fn,'/Step#0/Block/Pyy_3/0/');
    Pizz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pzz_3/0/');
    Pixy = hdf5read(fn,'/Step#0/Block/Pxy_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxy_3/0/');    
    Pixz = hdf5read(fn,'/Step#0/Block/Pxz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pxz_3/0/');
    Piyz = hdf5read(fn,'/Step#0/Block/Pyz_1/0/');%+hdf5read(fn,'/Step#0/Block/Pyz_3/0/');
    B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
    B2D=sqrt(Bx.^2+By.^2);
    perp2x=Bz.*Bx./(B.*B2D);
    perp2y=Bz.*By./(B.*B2D);
    perp2z=-B2D./B;
    Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

    [Pexx,Peyy,Pezz,Pexy,Pexz,Peyz]=compute_pressure(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Jex,Jey,Jez,rhoe, qom(1));
    [Pixx,Piyy,Pizz,Pixy,Pixz,Piyz]=compute_pressure(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Jix,Jiy,Jiz,rhoi, qom(2));
    
    Qex = hdf5read(fn,'/Step#0/Block/EFx_0/0/');%+hdf5read(fn,'/Step#0/Block/EFx_2/0/');
    Qey = hdf5read(fn,'/Step#0/Block/EFy_0/0/');%+hdf5read(fn,'/Step#0/Block/EFy_2/0/');
    Qez = hdf5read(fn,'/Step#0/Block/EFz_0/0/');%+hdf5read(fn,'/Step#0/Block/EFz_2/0/');    
    Qix = hdf5read(fn,'/Step#0/Block/EFx_1/0/');%+hdf5read(fn,'/Step#0/Block/EFx_3/0/');
    Qiy = hdf5read(fn,'/Step#0/Block/EFy_1/0/');%+hdf5read(fn,'/Step#0/Block/EFy_3/0/');
    Qiz = hdf5read(fn,'/Step#0/Block/EFz_1/0/');%+hdf5read(fn,'/Step#0/Block/EFz_3/0/'); 
  
    [Qenthex,Qenthey,Qenthez,Qbulkex,Qbulkey,Qbulkez,Qhfex,Qhfey,Qhfez,Ubulke,Uthe] = ...
    compute_energy_fluxes(Pexx,Peyy,Pezz,Pexy,Pexz,Peyz,Qex,Qey,Qez,Jex,Jey,Jez,rhoe, qom(1));

    [Qenthix,Qenthiy,Qenthiz,Qbulkix,Qbulkiy,Qbulkiz,Qhfix,Qhfiy,Qhfiz,Ubulki,Uthi] = ...
    compute_energy_fluxes(Pixx,Piyy,Pizz,Pixy,Pixz,Piyz,Qix,Qiy,Qiz,Jix,Jiy,Jiz,rhoi, qom(2));

end

Lx=25;Ly=50;Lz=1;
[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

[X Y] = meshgrid(0:dx:Lx,0:dy:Ly);

qom_ele = qom(1);

bufferX=round(Nx/20);
bufferY=round(Ny/20);
ir=3:Nx-bufferX;
jr=bufferY:Ny-bufferY;



%
% Electrons
%


global color_choice symmetric_color labelx labely labelc reversex Ncycle cyl Nsm_div

cyl=1;
reversex=0;signx=1;
symmetric_color=1;
color_choice =3;
labelx ='r/d_i';
labely ='\zeta/d_i';
labelc = 'mW/m^2';
lablec = 'c.u.';
Nsm=10
Nsm_div=0;

divfactor=permute(x,[2 1 3]);
%divfactor=ones(size(divfactor));



code_E=1;code_J=1;code_dp=1.;code_B=1;mu0=1;eps0=1;

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);



for iz=1%135
kr=0
kr=kr+round(iz);

AAz=vecpot(xc,yc,mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));


[Lx, Ly, Lz] = cross_prod(Jex, Jey, Jez, Bx, By, Bz);


labelc = 'c.u.';range=[-1 1]*1e-7
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(Lx(ir,jr,kr),3),AAz(ir,jr) ,['Lorentz Fx Y=' num2str((z(1,1,iz)))],'Fx',range, Nsm, 1,1,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(Ly(ir,jr,kr),3),AAz(ir,jr) ,['Lorentz Fy Y=' num2str((z(1,1,iz)))],'Fy',range, Nsm, 2,1,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(Lz(ir,jr,kr),3),AAz(ir,jr) ,['Lorentz Fz Y=' num2str((z(1,1,iz)))],'Fz',range, Nsm, 3,1,3);


[divPx, divPy, divPz, divPx1, divPx2] = compute_divtensor_cyl(x,y,z,Pexx,Pexy,Pexz,Peyy,Peyz,Pezz);

tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(divPx(ir,jr,kr),3),AAz(ir,jr) ,['Pressure Fx Y=' num2str((z(1,1,iz)))],'Fx',range, Nsm, 1,2,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(divPy(ir,jr,kr),3),AAz(ir,jr) ,['Pressure Fy Y=' num2str((z(1,1,iz)))],'Fy',range, Nsm, 2,2,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(divPz(ir,jr,kr),3),AAz(ir,jr) ,['Pressure Fz Y=' num2str((z(1,1,iz)))],'Fz',range, Nsm, 3,2,3);

tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(divPx1(ir,jr,kr),3),AAz(ir,jr) ,['divPx1 Fx Y=' num2str((z(1,1,iz)))],'Fz',range, Nsm, 10,1,2);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(divPx2(ir,jr,kr),3),AAz(ir,jr) ,['divPx2 Fx Y=' num2str((z(1,1,iz)))],'Fz',range, Nsm, 10,2,2);

rhoep=rhoe-1e-10;
[ugradux, ugraduy, ugraduz, ugradux1, ugradux2] = compute_AgradB_cyl(x,y,z,Jex./rhoep,Jey./rhoep,Jez./rhoep,Jex./rhoep,Jey./rhoep,Jez./rhoep);
ugradux=ugradux.*rhoe;
ugradux1=ugradux1.*rhoe;
ugradux2=ugradux2.*rhoe;
ugraduy=ugraduy.*rhoe;
ugraduz=ugraduz.*rhoe;
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(ugradux(ir,jr,kr),3),AAz(ir,jr) ,['UgradU Fx Y=' num2str((z(1,1,iz)))],'Fx',range, Nsm, 1,3,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(ugraduy(ir,jr,kr),3),AAz(ir,jr) ,['UgradU Fy Y=' num2str((z(1,1,iz)))],'Fy',range, Nsm, 2,3,3);
tmp=common_image_multi((X(jr,ir)),(Y(jr,ir)),mean(ugraduz(ir,jr,kr),3),AAz(ir,jr) ,['UgradU Fz Y=' num2str((z(1,1,iz)))],'Fz',range, Nsm, 3,3,3);


figure(1)
print('-dpng','-r300',['Forcex' Ncycle '.png'])
figure(2)
print('-dpng','-r300',['Forcey' Ncycle '.png'])
figure(3)
print('-dpng','-r300',['Forcez' Ncycle '.png'])
figure(10)
print('-dpng','-r300',['TwotermsdivPx' Ncycle '.png'])

figure(4)
subplot(2,1,1)
jplot=round(Ny/2-Ny/10);
plot(X(jplot,ir)',mean(Lx(ir,jplot,kr),3),X(jplot,ir)',mean(divPx(ir,jplot,kr),3),X(jplot,ir)',mean(ugradux(ir,jplot,kr),3))
legend('Lorentz','divP','\rho UgradU')
title(['direction r   Y=' num2str(Y(jplot,1))])
subplot(2,1,2)
plot(X(jplot,ir)',mean(divPx1(ir,jplot,kr),3),X(jplot,ir)',mean(divPx2(ir,jplot,kr),3))
hold on
plot(X(jplot,ir)',mean(ugradux1(ir,jplot,kr),3),X(jplot,ir)',mean(ugradux2(ir,jplot,kr),3))

legend('divP1','div2','\rho ugradu1','\rho ugradu2')
print('-dpng','-r300',['traces' Ncycle '.png'])

end



end