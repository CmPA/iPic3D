global blowup contours
blowup=0;
contours=0;

dir='/shared02/gianni/maha2/data1/'
cycle=37000

code_E = 2060.21;
code_B = 6.87213e-06;
code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;


addpath '/home/gianni/matlab2/matlab-parsek'
filename=[ dir 'settings.hdf'];
Lx=hdf5read(filename,'/collective/Lx'); 
Ly=hdf5read(filename,'/collective/Ly');
Nx=hdf5read(filename,'/collective/Nxc'); 
Ny=hdf5read(filename,'/collective/Nyc');
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
x=[-15 -45]; 
y=[-9 3];


legge_B=0
legge_J=0
legge_J0=0
legge_E=0
legge_V=1

if(legge_B)
fid= fopen('B_x_cycle55000.gda','rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny);

fid= fopen('B_y_cycle55000.gda','rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny);

fid= fopen('B_z_cycle55000.gda','rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny);

Ay=zeros(size(Bx));
if (contours) 
Ay=vecpot_uniform(xc,yc,Bx,By);
end

immagine(x,y,Bz*code_B,'Bz50000',[-10 10],0)

clear Bx By Bz
end

if(legge_J)
close all

fid= fopen('Je_z_cycle19000.gda','rb');
Je=fread(fid,'real*8');
fclose(fid);
Je=reshape(Je,Nx,Ny);


fid= fopen('Ji_z_cycle19000.gda','rb');
Ji=fread(fid,'real*8');
fclose(fid);
Ji=reshape(Ji,Nx,Ny);


%h = fspecial('gaussian',[5 5], .5);
%J = filter2(h, Je+Ji);

immagine(x,y,(Je+Ji)*code_J,'Jz19000',[-2 2],3)


%coplot_uniform(x,y,(J')*code_J,Ay','x/d_i','y/d_i','')
%caxis([-5 1])
%pcolor(x,y,Bx')


%set(gcf, 'Renderer', 'zbuffer');
%print('-dpng','-r300','Jy_1900.png')


end

if(legge_J0)
close all

fid= fopen('Je_z_cycle200.gda','rb');
Je=fread(fid,'real*8');
fclose(fid);
Je=reshape(Je,Nx,Ny);

fid= fopen('Ji_z_cycle200.gda','rb');
Ji=fread(fid,'real*8');
fclose(fid);
Ji=reshape(Ji,Nx,Ny);

immagine(x,y,(Je+Ji)*code_J,'Jz200',[-2 2],3)
end

if(legge_E)
close all

fid= fopen('E_y_cycle19000.gda','rb');
E=fread(fid,'real*8');
fclose(fid);
E=reshape(E,Nx,Ny);


immagine(x,y,E*code_E,'Ey19000',[-1e-2 1e-2],6)
end

if(legge_V)
close all

filename=[dir 'Je_x_cycle19000.gda']
fid= fopen(filename,'rb')
Je=fread(fid,'real*8');
fclose(fid);
Je=reshape(Je,Nx,Ny);
fid= fopen([dir 'rho_0_cycle19000.gda'],'rb');
Ne=fread(fid,'real*8');
fclose(fid);
Ne=reshape(Ne,Nx,Ny);

fid= fopen([dir 'Ji_x_cycle19000.gda'],'rb');
Ji=fread(fid,'real*8');
fclose(fid);
Ji=reshape(Ji,Nx,Ny);
fid= fopen([dir 'rho_1_cycle19000.gda'],'rb');
Ni=fread(fid,'real*8');
fclose(fid);
Ni=reshape(Ni,Nx,Ny);

%h = fspecial('gaussian',[5 5], .5);
%J = filter2(h, Je+Ji);

immagine(x,y,-Je./Ne*code_V,'Vex19000',[-2e3 2e3],3)
immagine(x,y,-Ji./Ni*code_V,'Vix19000',[-800 800],3)

%coplot_uniform(x,y,(J')*code_J,Ay','x/d_i','y/d_i','')
%caxis([-5 1])
%pcolor(x,y,Bx')


%set(gcf, 'Renderer', 'zbuffer');
%print('-dpng','-r300','Jy_1900.png')


end
