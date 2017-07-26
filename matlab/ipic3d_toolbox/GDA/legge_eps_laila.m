clear all
close all

dir='/shared02/gianni/lailaO1/data6/'
cycle=48000
dir='/shared02/gianni/lailaO1/data5/'
cycle=40000
%dir='/shared02/gianni/lailaO2/data3/'
%cycle=37000
dir='/shared02/gianni/StefanO12/data2/'
cycle=50000

ncycle=num2str(cycle);

global blowup contours
blowup=0;
contours=1;

code_E = 2060.21;
code_B = 6.87213e-06;
code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;
e=1.6e-19;

%
% Prepare to read
%

addpath '/home/gianni/matlab2/matlab-parsek'
filename=[dir 'settings.hdf'];
Lx=hdf5read(filename,'/collective/Lx'); 
Ly=hdf5read(filename,'/collective/Ly');
Nx=hdf5read(filename,'/collective/Nxc'); 
Ny=hdf5read(filename,'/collective/Nyc');
x=[0 Lx]; 
y=[0 Ly];
xn=linspace(0,Lx,Nx);
yn=linspace(0,Ly,Ny);


%
% Read data
%


% read magnetic field

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny);

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny);

b2=Bx.^2+By.^2+Bz.^2;
b2inp=Bx.^2+By.^2;

%read ion slippage information

scippa=0
if (scippa)
file=[dir 'VNFi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
VNFix=fread(fid,'real*8');
fclose(fid);
VNFix=reshape(VNFix,Nx,Ny);

file=[dir 'VNFi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
VNFiy=fread(fid,'real*8');
fclose(fid);
VNFiy=reshape(VNFiy,Nx,Ny);

file=[dir 'VNFi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
VNFiz=fread(fid,'real*8');
fclose(fid);
VNFiz=reshape(VNFiz,Nx,Ny);

VNFi_unsm=sqrt(VNFix.^2+VNFiy.^2+VNFiz.^2);
Nsm=0
VNFi=smooth(VNFi_unsm,Nsm);
end

% read pressure tensor informaiton

file=[dir 'Pi_eps_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
EPS=fread(fid,'real*8');
fclose(fid);
EPS=reshape(EPS,Nx,Ny);


file=[dir 'Pi_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ppar=fread(fid,'real*8');
fclose(fid);
Ppar=reshape(Ppar,Nx,Ny);

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny);

file=[dir 'Pi_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper2=fread(fid,'real*8');
fclose(fid);
Pper2=reshape(Pper2,Nx,Ny);

%compute epsilon firehose locally

eps2inp=1-4*pi*(Ppar-sqrt(Pper1.*Pper2))./b2inp;


% plot only if ion slippage is requested

if(scippa)
immagine_dir(x,y,VNFi,['' ncycle],[0 .1],0)

xlim([150 400]) 
ylim([50 150])
name='VNFi'
print('-dpng','-r300',[name '.png'])


xlim([220 300]) 
ylim([102-15 102+15])
name='VNFi'
print('-dpng','-r300',[name '_bu.png'])

figure
% First index is x 

Nplot=2400:200:2800;

plot(yn,VNFi(Nplot,:))
%title(['x=' num2str(xn(Nplot))],'fontsize',[14])
xlim([50 150]) 
legend(num2str(round(xn(Nplot)')))
set(gca,'fontsize',[14])
xlabel('y/d_i','fontsize',[14])
ylabel('(u_i-u_{ExB})/c','fontsize',[14])
print('-depsc','-painters',[name '_cuts.eps'])
%hold on
%contour(xc,yc,fliplr(Ay'),50,'w')
end



% Plot of the epsilon of firehose

figure(1)

EPS_sm=immagine_dir(x,y,Bx,['\epsilon' ncycle],[0 0],0);

xlim([60 140]) 
ylim([0 Ly/2])

hold on 
%contour(xn,yn,EPS_sm',[ 0.25 0.25],'m')
Nplot=Nx/2-Nx/10:Nx/40:Nx/2-Nx/40;
Nplot=[Nplot Nx-Nplot]
for i=Nplot
xn(i)
plot( [xn(i) xn(i)], [0 Ly],'w')
end

name='Bx'
print('-dpng','-r300',[name '.png'])

test=0
if(test)
immagine_dir(x,y,eps2inp,['\epsilon' ncycle],[-.5 1.5],0)

%immagine_dir(x,y,eps2,['\epsilon' ncycle],[-10 10],0)

%xlim([150 400]) 
%ylim([50 150])
name='eps2'
print('-dpng','-r300',[name '.png'])
end

figure(2)
   
plot(yn,EPS(Nplot,:))
%title(['x=' num2str(xn(Nplot))],'fontsize',[14])
%xlim([80 120]) 
%hold on 
%plot([0 Ly], [.25 .25],'m--') 
%ylim([-2.5 1.5])
legend(num2str(round(xn(Nplot)')),'location','SouthEast')
set(gca,'fontsize',[14])
xlabel('y/d_i','fontsize',[14])
ylabel('\epsilon','fontsize',[14])
print('-depsc','-painters',[name '_cuts.eps'])


figure(3)
   
name='bx';
plot(yn,Bx(Nplot,:))
%title(['x=' num2str(xn(Nplot))],'fontsize',[14])
%xlim([80 120]) 
xlim([4 12])
hold on 
plot([0 Ly/2], [0 0],'m--') 
legend(num2str(round(xn(Nplot)')),'location','SouthEast')
set(gca,'fontsize',[14])
xlabel('y/d_i','fontsize',[14])
ylabel('B_x','fontsize',[14])
print('-depsc','-painters',[name '_cuts.eps'])


