clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Reconnection_results/sim5_lowT_th1_256ppc_dt0p1'
namefile = 'data-Fields';

% folder_name = '/shared/gianni/tred74_2D/data_rhofromJ_periodic/'
% namefile = 'GEM-Fields';

folder_name = '/Users/gianni/Dropbox/Science/san_diego/6-coil-Nov-2017/'
namefile = 'SixCoils-Fields';



Lx=70;
Ly=70;
Lz=70;



for i=4000:1000:4000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    old=0
    if(old)
    bx = hdf5read(fn,'/Step#0/Block/Btx/0/');
    by = hdf5read(fn,'/Step#0/Block/Bty/0/');
    bz = hdf5read(fn,'/Step#0/Block/Btz/0/');
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
    
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    AAz=vecpot(xc,yc,bx',by');AAz=AAz';

    h=figure(1)
    %set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    titolo=[ 'cycle=' num2str(i) '  B (color) B(contours)']
    
    
    range1=[-3 3]*0e-5; 
    range2=[0 0];
 
    %coplot(i,xg,yg,epar,AAz,xlab,ylab,titolo,range1, range2)

    divE=divergence(xc,yc,ex,ey);
    range1=[-1 1]*0e-2; 
%     figure(1)
%     coplot(i,xg,yg,divE,AAz,xlab,ylab,titolo,range1, range2)
%     figure(2)
%     coplot(i,xg,yg,rho,AAz,xlab,ylab,titolo,range1, range2)
%     figure(3)
%     coplot(i,xg,yg,4*pi*rho-divE,AAz,xlab,ylab,titolo,range1, range2)
    h=figure(1)
    set(h,'Position',[123 336 1264 469])
        subplot(1,2,1)
        titolo=[ 'cycle=' num2str(i) '  Ex (color) B(contours)']
    coplot(i,xg,yg,ex,AAz,xlab,ylab,titolo,range1, range2)
    subplot(1,2,2)
            titolo=[ 'cycle=' num2str(i) '  Bx (color) B(contours)']
    coplot(i,xg,yg,bx,AAz,xlab,ylab,titolo,range1, range2)
        print('-dpng', 'figure1.png')
        
        h=figure(2)
    set(h,'Position',[123 336 1264 469])
        subplot(1,2,1)
        titolo=[ 'cycle=' num2str(i) '  Ey (color) B(contours)']
    coplot(i,xg,yg,ey,AAz,xlab,ylab,titolo,range1, range2)
    subplot(1,2,2)
            titolo=[ 'cycle=' num2str(i) '  By (color) B(contours)']
    coplot(i,xg,yg,by,AAz,xlab,ylab,titolo,range1, range2)
            print('-dpng', 'figure2.png')
            
        h=figure(3)
    set(h,'Position',[123 336 1264 469])
        subplot(1,2,1)
        titolo=[ 'cycle=' num2str(i) '  Ez (color) B(contours)']
    coplot(i,xg,yg,ez,AAz,xlab,ylab,titolo,range1, range2)
    subplot(1,2,2)
            titolo=[ 'cycle=' num2str(i) '  Bz (color) B(contours)']
    coplot(i,xg,yg,bz,AAz,xlab,ylab,titolo,range1, range2)
            print('-dpng', 'figure3.png')
            
            h=figure(4)
    set(h,'Position',[123 336 1264 469])
        subplot(1,2,1)
        titolo=[ 'cycle=' num2str(i) '  rho_e (color) B(contours)']
    coplot(i,xg,yg,rho_e*4*pi,AAz,xlab,ylab,titolo,range1, range2)
    subplot(1,2,2)
            titolo=[ 'cycle=' num2str(i) '  rho_i (color) B(contours)']
    coplot(i,xg,yg,rho_i*4*pi,AAz,xlab,ylab,titolo,range1, range2)
            print('-dpng', 'figure4.png')
end 