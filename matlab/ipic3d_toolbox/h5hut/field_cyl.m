clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
%folder_name = '/Users/gianni/Dropbox/Science/codes/build_cyl/wb2d_b/'
namefile = 'TwoCoils-Fields';


Lx=37.5;
Ly=75;


i=35000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    
    p = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');
    p = p +hdf5read(fn,'/Step#0/Block/Pyy_1/0/');
    p = p +hdf5read(fn,'/Step#0/Block/Pzz_1/0/');
    p=permute(squeeze(p(:,:,round(Nz/2))),[2 1]);
    
    
    bx = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(1).Datasets(1));
    by = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1));
    bz = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(5).Datasets(1));
    bx_ext = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(2).Datasets(1));
    by_ext = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(4).Datasets(1));
    bz_ext = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(6).Datasets(1));
    bx=bx+bx_ext;
    by=by+by_ext;
    bz=bz+bz_ext;
    
  
    ex = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(7).Datasets(1));
    ey = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(8).Datasets(1));
    ez = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(9).Datasets(1));

    ex=permute(squeeze(ex(:,:,round(Nz/2))),[2 1]);
    ey=permute(squeeze(ey(:,:,round(Nz/2))),[2 1]);
    ez=permute(squeeze(ez(:,:,round(Nz/2))),[2 1]);
    
    bx=permute(squeeze(bx(:,:,round(Nz/2))),[2 1]);
    by=permute(squeeze(by(:,:,round(Nz/2))),[2 1]);
    bz=permute(squeeze(bz(:,:,round(Nz/2))),[2 1]);
    
    b = sqrt (bx.^2 +by.^2 + bz.^2);
    
   
     epar=dot(ex,ey,ez,bx,by,bz)./b;

    
    global ex ey ez bx by bz xg yg  Lx Ly qom Rout
    
    [xg,yg]=meshgrid(0:Nx-1,0:Ny-1);
    xg=xg/Nx*Lx;
    yg=yg/Ny*Ly;
    
    dx=Lx/Nx;
    dy=Ly/Ny;
    
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    ath=vecpot_cyl(xc,yc,bx,by);
    
    h=figure(1)
    %set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    titolo=[ 'cycle=' num2str(i) '  p (color) Az(contours)']
    
    
    range1=[-3 3]*0e-4; 
    range2=[0 0];
 
     coplot(i,xg,yg,p,ath,xlab,ylab,titolo,range1, range2)

