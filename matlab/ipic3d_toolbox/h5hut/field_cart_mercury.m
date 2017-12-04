clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Reconnection_results/sim5_lowT_th1_256ppc_dt0p1'
namefile = 'data-Fields';

% folder_name = '/shared/gianni/tred74_2D/data_rhofromJ_periodic/'
% namefile = 'GEM-Fields';

folder_name = '/Users/gianni/Dropbox/Science/ucla/mercury/data/'
namefile = 'Mercury3D-P1-Fields';



Ly=30.8
Lx=64
Lz=30.8



for i=13000:1000:13000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    old=1
    if(old)
    bx = hdf5read(fn,'/Step#0/Block/Btx/0/');
    by = hdf5read(fn,'/Step#0/Block/Bty/0/');
    bz = hdf5read(fn,'/Step#0/Block/Btz/0/');
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    %AAz=vecpot(xc,yc,bx',by');AAz=AAz';

    iz=floor(Nz/2)
    AAz=vecpot(xc,yc,bx(:,:,iz),by(:,:,iz));AAz=AAz';
    
    else
    bx = hdf5read(fn,'/Step#0/Block/Bx/0/');
    by = hdf5read(fn,'/Step#0/Block/By/0/');
    bz = hdf5read(fn,'/Step#0/Block/Bz/0/');
    bx_ext = hdf5read(fn,'/Step#0/Block/Bx_ext/0/');
    by_ext = hdf5read(fn,'/Step#0/Block/By_ext/0/');
    bz_ext = hdf5read(fn,'/Step#0/Block/Bz_ext/0/');
    bx=bx+bx_ext;
    by=by+by_ext;
    bz=bz+bz_ext;
    end
    
  
    ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    ez = hdf5read(fn,'/Step#0/Block/Ez/0/');

    ex=permute(squeeze(ex(:,:,round(Nz/2))),[2 1]);
    ey=permute(squeeze(ey(:,:,round(Nz/2))),[2 1]);
    ez=permute(squeeze(ez(:,:,round(Nz/2))),[2 1]);
    
    bx=permute(squeeze(bx(:,:,round(Nz/2))),[2 1]);
    by=permute(squeeze(by(:,:,round(Nz/2))),[2 1]);
    bz=permute(squeeze(bz(:,:,round(Nz/2))),[2 1]);
    
    b = sqrt (bx.^2 +by.^2 + bz.^2);
    
   
     epar=dot(ex,ey,ez,bx,by,bz)./b;
     rho_i = hdf5read(fn,'/Step#0/Block/rho_1/0/'); %+ hdf5read(fn,'/Step#0/Block/rho_3/0/') ; 
     rho_e = hdf5read(fn,'/Step#0/Block/rho_0/0/'); %+ hdf5read(fn,'/Step#0/Block/rho_2/0/') ; 
     rho=rho_e+rho_i;
     rho=permute(squeeze(rho(:,:,round(Nz/2))),[2 1]);
     rho_i=permute(squeeze(rho_i(:,:,round(Nz/2))),[2 1]);
     rho_e=permute(squeeze(rho_e(:,:,round(Nz/2))),[2 1]);
    
    %p = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');
   % p = p +hdf5read(fn,'/Step#0/Block/Pyy_1/0/');
    %p = p +hdf5read(fn,'/Step#0/Block/Pzz_1/0/');
    global ex ey ez bx by bz xg yg  Lx Ly qom Rout
    
    [xg,yg]=meshgrid(0:Nx-1,0:Ny-1);
    xg=xg/Nx*Lx;
    yg=yg/Ny*Ly;
    
    dx=Lx/Nx;
    dy=Ly/Ny;
    
    
    
    h=figure(1)
    %set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    titolo=[ 'cycle=' num2str(i) '  B (color) B(contours)']
    
    
    range1=[-3 3]*0e-5; 
    range2=[0 0];
 

    h=figure(1)
    set(h,'Position',[123 336 1264 469])
        subplot(1,2,1)
        titolo=[ 'cycle=' num2str(i) '  Ex (color) B(contours)']
    coplot(i,xg,yg,bx,AAz,xlab,ylab,titolo,range1, range2)
    subplot(1,2,2)
            titolo=[ 'cycle=' num2str(i) '  Bx (color) B(contours)']
    coplot(i,xg,yg,by,AAz,xlab,ylab,titolo,range1, range2)
        print('-dpng', 'figure1.png')
  
end 