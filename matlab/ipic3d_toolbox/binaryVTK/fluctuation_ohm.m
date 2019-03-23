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

colormap hsv


nE(AAz,rhoe,Ex,Ey, Ez,1,'rhoE',ncycle, list_value)
JxB(AAz,Jex,Jey,Jez ,Ex,Ey, Ez,2,'JxB',ncycle, list_value)




function [] = nE(a,p,q1,q2,q3,n,name,prename, list_value)
% MYMEAN Example of a local function.
close all

ndiv=300;
Np=max(size(a(:)));
p_avg=mean(p,3);
[Nx Ny Nz]=size(p);
dp=p;
for k=1:Nz
    dp(:,:,k)=p(:,:,k)-p_avg;
end

q1_avg=mean(q1,3);
dq1=q1;
for k=1:Nz
    dq1(:,:,k)=q1(:,:,k)-q1_avg;
end
q2_avg=mean(q2,3);
dq2=q2;
for k=1:Nz
    dq2(:,:,k)=q2(:,:,k)-q2_avg;
end
q3_avg=mean(q3,3);
dq3=q3;
for k=1:Nz
    dq3(:,:,k)=q3(:,:,k)-q3_avg;
end
figure(n)

dpq=dp.*dq1;
%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),dpq(:),ones(Np,1),0,min(a(:)),max(a(:)),min(dpq(:)),max(dpq(:)),ndiv);
%int1=urange*nbinpq./sum(nbinpq,1);
[int1, xrange] = spaziofasi_int(a(:),dpq(:),ones(Np,1),ndiv);

dpq=dp.*dq2;
%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),dpq(:),ones(Np,1),0,min(a(:)),max(a(:)),min(dpq(:)),max(dpq(:)),ndiv);
%int2=urange*nbinpq./sum(nbinpq,1);
int2 = spaziofasi_int(a(:),dpq(:),ones(Np,1),ndiv);

dpq=dp.*dq3;
%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),dpq(:),ones(Np,1),0,min(a(:)),max(a(:)),min(dpq(:)),max(dpq(:)),ndiv);
%int3=urange*nbinpq./sum(nbinpq,1);
int3 = spaziofasi_int(a(:),dpq(:),ones(Np,1),ndiv);

figure(n) 
plot(xrange,-int1,xrange,-int2,xrange,-int3);
legend('x','y','z')
set(gca,'fontsize',[14])
print('-dpng','-r300',[prename '_sum_' name])
%close(n)
end

function [] = JxB(a,p1,p2, p3,q1,q2,q3,n,name,prename, list_value)
% MYMEAN Example of a local function.
close all

ndiv=300;
Np=max(size(a(:)));
p1_avg=mean(p1,3);
[Nx Ny Nz]=size(p1);
dp1=p1;
for k=1:Nz
    dp1(:,:,k)=p1(:,:,k)-p1_avg;
end
p2_avg=mean(p2,3);
dp2=p2;
for k=1:Nz
    dp2(:,:,k)=p2(:,:,k)-p2_avg;
end
p3_avg=mean(p3,3);
dp3=p3;
for k=1:Nz
    dp3(:,:,k)=p3(:,:,k)-p3_avg;
end

q1_avg=mean(q1,3);
dq1=q1;
for k=1:Nz
    dq1(:,:,k)=q1(:,:,k)-q1_avg;
end
q2_avg=mean(q2,3);
dq2=q2;
for k=1:Nz
    dq2(:,:,k)=q2(:,:,k)-q2_avg;
end
q3_avg=mean(q3,3);
dq3=q3;
for k=1:Nz
    dq3(:,:,k)=q3(:,:,k)-q3_avg;
end


[cp1, cp2, cp3] = cross_prod(dp1, dp2, dp3, dq1, dq2, dq3);
figure(n)

%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),cp1(:),ones(Np,1),0,min(a(:)),max(a(:)),min(cp1(:)),max(cp1(:)),ndiv);
%int1=(urange*nbinpq)./sum(nbinpq,1);
[int1, xrange] = spaziofasi_int(a(:),cp1(:),ones(Np,1),ndiv);


%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),cp2(:),ones(Np,1),0,min(a(:)),max(a(:)),min(cp2(:)),max(cp2(:)),ndiv);
%int2=urange*nbinpq./sum(nbinpq,1);
int2 = spaziofasi_int(a(:),cp2(:),ones(Np,1),ndiv);

%[totnum,nbinpq,xrange,urange]=spaziofasi2(a(:),cp3(:),ones(Np,1),0,min(a(:)),max(a(:)),min(cp3(:)),max(cp3(:)),ndiv);
%int3=urange*nbinpq./sum(nbinpq,1);
int3 = spaziofasi_int(a(:),cp3(:),ones(Np,1),ndiv);

figure(n)  
plot(xrange,-int1,xrange,-int2,xrange,-int3);

legend('x','y','z')
set(gca,'fontsize',[14])
print('-dpng','-r300',[prename '_sum_' name])
%close(n)
end


