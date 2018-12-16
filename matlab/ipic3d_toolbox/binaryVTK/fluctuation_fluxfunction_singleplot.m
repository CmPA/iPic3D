close all
clear all
!rm *.png
addpath(genpath('../../ipic3d_toolbox'));

must_read=true;
leggo='h5';
if(must_read)

sim_name='tred82'
switch sim_name
case 'tred77'
TRED77;
case_name='GEM';
cycle = 15000;
zcode = Lz/2;
case 'tred81'
tred81;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;
    case 'tred82'
tred82;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;
case 'AH'
generic;
case_name='AH';
cycle =4000;
zcode = Lz/2;
case 'HRmaha3D3'
HRmaha3D3;
leggo='gda';
    case_name='GEM';
dir='/data1/gianni/HRmaha3D3/h5/'; cycle= 80002; ncycle = num2str(cycle,'%06d');
cycle = 80002;  % for h5
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
case '7feb09'
FEB09;
cycle=18000
case_name='MHDUCLA'
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
otherwise
print('no recognised case selected')
end

% Prepare string
ntime = num2str(cycle,'%06d');
ncycle = num2str(cycle,'%06d');


import_h5_binvtk
end


[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
S=sqrt(Sx.^2+Sy.^2+Sz.^2);
[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);
radius=0.01;
divS = compute_div(x,y,z,Sx,Sy,Sz,radius, 0.);
divQe = compute_div(x,y,z,Qex,Qey,Qez,radius, 0.);
divQi = compute_div(x,y,z,Qix,Qiy,Qiz,radius, 0.);

UPex = (Jex.*Pexx + Jey.*Pexy + Jez.*Pexz)./rhoe;
UPey = (Jex.*Pexy + Jey.*Peyy + Jez.*Peyz)./rhoe;
UPez = (Jex.*Pexz + Jey.*Peyz + Jez.*Pezz)./rhoe;
UPix = (Jix.*Pixx + Jiy.*Pixy + Jiz.*Pixz)./rhoi;
UPiy = (Jix.*Pixy + Jiy.*Piyy + Jiz.*Piyz)./rhoi;
UPiz = (Jix.*Pixz + Jiy.*Piyz + Jiz.*Pizz)./rhoi;
radius=0.001;
[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradVe, ugradp, pdivve, divUP] = compute_energy_balance( ...
    rhoe, Jex, Jey, Jez,... 
    Qbulkex, Qbulkey, Qbulkez, Qenthex, Qenthey, Qenthez, Qhfex, Qhfey, Qhfez, ...
    Pexx, Peyy, Pezz, Pexy, Pexz, Peyz, x, y, z, dx, dy, dz, qom, radius,0.);
[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradVi, ugradp, pdivvi, divUP] = compute_energy_balance( ...
    rhoi, Jix, Jiy, Jiz,... 
    Qbulkix, Qbulkiy, Qbulkiz, Qenthix, Qenthiy, Qenthiz, Qhfix, Qhfiy, Qhfiz, ...
    Pixx, Piyy, Pizz, Pixy, Pixz, Piyz, x, y, z, dx, dy, dz, 1.0, radius,0.);


switch sim_name
case 'tred82'
    list_value=[-.02, -.01, 0, .005, .01, .015, .02] %tred82
    otherwise
    list_value=-.02:.01:.04 %tred81
end

xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=zeros(size(Bx));
for kr=1:Nz
AAz(:,:,kr)=vecpot(xc,yc,Bx(:,:,kr),By(:,:,kr));
AAz(:,:,kr)=AAz(:,:,kr)-AAz(round(Nx/2),round(Ny/2),kr);
end
figure(1)
imagesc(xc,yc,mean(AAz,3)')
%contourf(xc,yc,mean(AAz,3)',20)
%hold on; contour(xc,yc,mean(AAz,3)',list_value,'m','linewidth',2)
load cm_new; colormap(cm_kbwrk)
%colormap default
colorbar
xlabel('x')
ylabel('y')
cmax=max(max(abs(mean(AAz,3))));
caxis([-cmax cmax])
axis equal 
axis tight
print('-dpng','-r300',[ncycle '_Phi'])


% 
% [X,Y,Z]=ndgrid(1:Nx,1:Ny,1:Nz);
% figure
% plot3(AAz(:),Y(:),S(:),'.')
% figure
% plot(AAz(:),S(:),'.')

colormap hsv

figura(AAz,S,2,'S',ncycle, list_value)
figura(AAz,divS,2,'divS',ncycle, list_value)
figura(AAz,divQe,2,'divQe',ncycle, list_value)
figura(AAz,divQi,2,'divQi',ncycle, list_value)
figura(AAz,Bx.*Bx+By.*By+Bz.*Bz,2,'B2',ncycle, list_value)
figura(AAz,Ex.*Ex+Ey.*Ey+Ez.*Ez,2,'E2',ncycle, list_value)
figura(AAz,0.5*(Jex.^2+Jey.^2+Jez.^2)./rhoe/qom,2,'Ubulke',ncycle, list_value)
figura(AAz,0.5*(Pexx+Peyy+Pezz),2,'Uthe',ncycle, list_value);
figura(AAz,0.5*(Jix.^2+Jiy.^2+Jiz.^2)./rhoi,2,'Ubulki',ncycle, list_value)
figura(AAz,0.5*(Pixx+Piyy+Pizz),2,'Uthi',ncycle, list_value);
%figura(AAz,log10(S),3,'Slog',ncycle, list_value)
figura(AAz,Sx,4,'Sx',ncycle, list_value)
figura(AAz,Sy,5,'Sy',ncycle, list_value)
figura(AAz,Sz,6,'Sz',ncycle, list_value)
figura(AAz,Qex,4,'Qex',ncycle, list_value)
figura(AAz,Qey,4,'Qey',ncycle, list_value)
figura(AAz,Qez,4,'Qez',ncycle, list_value)
figura(AAz,Qix,4,'Qix',ncycle, list_value)
figura(AAz,Qiy,4,'Qiy',ncycle, list_value)
figura(AAz,Qiz,4,'Qiz',ncycle, list_value)
%figura(AAz,UPex,4,'uPex',ncycle, list_value)
%figura(AAz,UPey,4,'uPey',ncycle, list_value)
%figura(AAz,UPez,4,'uPez',ncycle, list_value)
%figura(AAz,UPix,4,'uPix',ncycle, list_value)
%figura(AAz,UPiy,4,'uPiy',ncycle, list_value)
%figura(AAz,UPiz,4,'uPiz',ncycle, list_value)
figura(AAz,pdivve,4,'pthetae',ncycle, list_value)
figura(AAz,pdivvi,4,'pthetai',ncycle, list_value)
figura(AAz,PgradVe-pdivve,4,'PiDe',ncycle, list_value)
figura(AAz,PgradVi-pdivvi,4,'PiDi',ncycle, list_value)
figura(AAz,Jex.*Ex+Jey.*Ey+Jez.*Ez,4,'JeE',ncycle, list_value)
figura(AAz,Jix.*Ex+Jiy.*Ey+Jiz.*Ez,4,'JiE',ncycle, list_value)

        
        Epx = Ex + (Jey.*Bz - Jez.*By)./rhoe;
        Epy = Ey + (Jez.*Bx - Jex.*Bz)./rhoe;
        Epz = Ez + (Jex.*By - Jey.*Bx)./rhoe;
        
        JdotEp=(Jex+Jix).*Epx + (Jey+Jiy).*Epy + (Jez+Jiz).*Epz;
figura(AAz,JdotEp,4,'JEp',ncycle, list_value)
     
xflow='inflow'
   structure_function
xflow='outflow'
   structure_function

        
function [] = figura(a,p,n,name,prename, list_value)
% MYMEAN Example of a local function.
close all

ndiv=100;
Np=max(size(a(:)));
p_avg=mean(p,3);
[Nx Ny Nz]=size(p);
dp=p;
for k=1:Nz
    dp(:,:,k)=p(:,:,k)-p_avg;
end
% figure(n)
% [totnum,nbinu,xrange,urange]=spaziofasi2(a(:),p(:),ones(Np,1),0,min(a(:)),max(a(:)),min(p(:)),max(p(:)),ndiv);
% imagesc(xrange,urange,log10(nbinu))
% xlabel('\Phi')
% ylabel(name)
% colorbar
% colormap hsv
% print('-dpng','-r300',[prename '_' name])
% close(n)
figure(n)
[totnum,nbinu,xrange,urange]=spaziofasi2(a(:),dp(:),ones(Np,1),0,min(a(:)),max(a(:)),min(dp(:)),max(dp(:)),ndiv);
durange=urange(2)-urange(1);
imagesc(xrange,urange,log10(nbinu))
xlabel('\Phi','fontsize',[14])
ylabel(name,'fontsize',[14])
colorbar
%colormap hsv
load gist_ncar.mat
colormap(flipud(gist_ncar))

Ncuts=7;
figure(n+1)
        urka=-20:.1:20;
        semilogy(urka,exp(-urka.^2/2)/sqrt(2*pi),'k--')
hold on
labelle=["normal"];

Nxr=max(size(xrange));

for i=1:Ncuts
%ip=round(ndiv/Ncuts*i-ndiv/Ncuts/2);
ip=(list_value(i)*(Nxr-1)+max(xrange)-min(xrange)*Nxr)/(max(xrange)-min(xrange));
ip=round(ip);
lr=num2str(xrange(ip),'%10.3f\n');
labelle=[labelle;string(lr)];
figure(n)
hold on
plot([xrange(ip) xrange(ip)],[min(urange) max(urange)],'k')
figure(n+1)
sig=sqrt(urange.^2*nbinu(:,ip)/sum(nbinu(:,ip)));
        semilogy(urange/sig,nbinu(:,ip)./sum(durange/sig*nbinu(:,ip)))%,'linewidth',[4])
%ylim([min(nbinu(:,ip)) max(nbinu(:,ip))])
end
ylim([1e-6, 10])
xlim([-20 20])
xlabel(['\Delta' name],'fontsize',[14])
title(name,'fontsize',[14])
legend(labelle)
set(gca,'fontsize',[14])
print('-dpng','-r300',[prename '_mp_' name])
close(n+1)
figure(n)  
set(gca,'fontsize',[14])
print('-dpng','-r300',[prename '_d_' name])
close(n)

end


