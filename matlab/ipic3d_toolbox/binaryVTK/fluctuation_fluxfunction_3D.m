close all
addpath(genpath('../../ipic3d_toolbox'));

must_read=false;
if(must_read)

sim_name='tred81'
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
case 'AH'
generic;
case_name='AH';
cycle =4000;
zcode = Lz/2;
case 'HRmaha3D3'
HRmaha3D3;
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

xc=linspace(0, Lx, Nx+1);
yc=linspace(0, Ly, Ny+1);
AAz=zeros(size(Bx));
for kr=1:Nz
AAz(:,:,kr)=vecpot(xc,yc,Bx(:,:,kr),By(:,:,kr));
AAz(:,:,kr)=AAz(:,:,kr)-AAz(round(Nx/2),round(Ny/2),kr);
end
figure(1)
imagesc(xc,yc,mean(AAz,3)')
load cm_new
colormap(cm_kbwrk)
colorbar
xlabel('x')
ylabel('y')
cmax=max(max(abs(mean(AAz,3))));
caxis([-cmax cmax])
print('-dpng','-r300',[ncycle '_Phi'])
% 
% [X,Y,Z]=ndgrid(1:Nx,1:Ny,1:Nz);
% figure
% plot3(AAz(:),Y(:),S(:),'.')
% figure
% plot(AAz(:),S(:),'.')

colormap hsv

figura(AAz,S,2,'S',ncycle)
figura(AAz,log10(S),3,'Slog',ncycle)
figura(AAz,Sx,4,'Sx',ncycle)
[X,Y,Z]=ndgrid(1:Nx,1:Ny,1:Nz);
figura3D(AAz,Sx,Y,4,'Sx',ncycle)
figura(AAz,Sy,5,'Sy',ncycle)
figura(AAz,Sz,6,'Sz',ncycle)

function [] = figura(a,p,n,name,prename)
% MYMEAN Example of a local function.
close all
figure(n)
ndiv=100;
Np=max(size(a(:)));
p_avg=mean(p,3);
[Nx Ny Nz]=size(p);
dp=p;
for k=1:Nz
    dp(:,:,k)=p(:,:,k)-p_avg;
end 
[totnum,nbinu,xrange,urange]=spaziofasi2(a(:),p(:),ones(Np,1),0,min(a(:)),max(a(:)),min(p(:)),max(p(:)),ndiv);
imagesc(xrange,urange,log10(nbinu))
xlabel('\Phi')
ylabel(name)
colorbar
colormap hsv
print('-dpng','-r300',[prename '_' name])
close(n)
figure(n)
[totnum,nbinu,xrange,urange]=spaziofasi2(a(:),dp(:),ones(Np,1),0,min(a(:)),max(a(:)),min(dp(:)),max(dp(:)),ndiv);
imagesc(xrange,urange,log10(nbinu))
xlabel('\Phi')
ylabel(name)
colorbar
colormap hsv
print('-dpng','-r300',[prename '_d_' name])
close(n)
end

function [] = figura3D(a,p,x,n,name,prename)
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
[vdf,umin,vmin,wmin,du,dv,dw]=spaziofasi3D_different(a(:),x(:),p(:),ones(Np,1),ndiv);

savevtk_bin(log10(vdf+1e-20),[prename '_VTK_' name '.vtk'],'vdf',1,1,1,0,0,0);

end
