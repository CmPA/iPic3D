clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Science/san_diego/high-res-steady-state'
namefile = 'TwoCoils2D-Fields';


Lx=25;
Ly=50;
qom_0 = -40;
qom_1 = 2.5;


i = 72000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    n = hdf5read(fn,'/Step#0/Block/rho_0/0/');
    jx = hdf5read(fn,'/Step#0/Block/Jx_0/0/');
    pxx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/');
    pxx = (pxx - jx.*jx./(n+1e-10))/qom_0;
    jy = hdf5read(fn,'/Step#0/Block/Jy_0/0/');
    pyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/');
    pyy = (pyy - jy.*jy./(n+1e-10))/qom_0;
    jz = hdf5read(fn,'/Step#0/Block/Jz_0/0/');
    pzz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/');
    pzz = (pzz - jz.*jz./(n+1e-10))/qom_0;
    
    pe = pxx + pyy + pzz;
    
    vetheta=jz./(n+1e-10);
    ve=sqrt(jx.^2+jz.^2+jz.^2)./abs(n+1e-10);
    
    n = hdf5read(fn,'/Step#0/Block/rho_1/0/');
    jx = hdf5read(fn,'/Step#0/Block/Jx_1/0/');
    pxx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');
    pxx = (pxx - jx.*jx./(n+1e-10))/qom_1;
    jy = hdf5read(fn,'/Step#0/Block/Jy_1/0/');
    pyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/');
    pyy = (pyy - jy.*jy./(n+1e-10))/qom_1;
    jz = hdf5read(fn,'/Step#0/Block/Jz_1/0/');
    pzz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/');
    pzz = (pzz - jz.*jz./(n+1e-10))/qom_1;
    
    pion = pxx + pyy + pzz;
    
    p = pe + pion;
    
    vtheta=jz./(n+1e-10);
    v=sqrt(jx.^2+jz.^2+jz.^2)./(n+1e-10);
    
    p=permute(squeeze(p(:,:,round(Nz/2))),[2 1]);
    
    n=permute(squeeze(n(:,:,round(Nz/2))),[2 1]);
    
    vtheta=permute(squeeze(vtheta(:,:,round(Nz/2))),[2 1]);
    v=permute(squeeze(v(:,:,round(Nz/2))),[2 1]);
    vetheta=permute(squeeze(vetheta(:,:,round(Nz/2))),[2 1]);
    ve=permute(squeeze(ve(:,:,round(Nz/2))),[2 1]);
    
    bx = hdf5read(fn,'/Step#0/Block/Bx/0/');
    by = hdf5read(fn,'/Step#0/Block/By/0/');
    bz = hdf5read(fn,'/Step#0/Block/Bz/0/');
    bx_ext = hdf5read(fn,'/Step#0/Block/Bx_ext/0/');
    by_ext = hdf5read(fn,'/Step#0/Block/By_ext/0/');
    bz_ext = hdf5read(fn,'/Step#0/Block/Bz_ext/0/');
    bx=bx+bx_ext;
    by=by+by_ext;
    bz=bz+bz_ext;
    
  
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
    ben = b.^2/2/4/pi*qom_1;
   
     epar=dot(ex,ey,ez,bx,by,bz)./b;

    
    global contours dx dy color_choice
    
    [xg,yg]=meshgrid(0:Nx-1,0:Ny-1);
    xg=xg/Nx*Lx;
    yg=yg/Ny*Ly;
    
    dx=Lx/Nx;
    dy=Ly/Ny;
    
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
   % ath=vecpot_uniform(xc',yc',bx',by');ath=ath';
    ath=vecpot_cyl(xc,yc,bx,by);
    
    h=figure(1)
    %set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    
    
    
    range1=[-1 1]*2e-3*5; 
    range2=[0 0];
    color_choice = 3
 
titolo=[ 'cycle=' num2str(i) '  ve_t(int) (color) A(contours)']

  coplot(i,xg,yg,vetheta,ath,xlab,ylab,titolo,range1, range2,1)
      range1=[0 1]*1e-3*5; 
    range2=[0 0];
    titolo=[ 'cycle=' num2str(i) '  ve (color) A(contours)']
    color_choice = 0
  coplot(i,xg,yg,ve,ath,xlab,ylab,titolo,range1, range2,2)
  
  
    titolo=[ 'cycle=' num2str(i) '  p*r (color) A(contours)']
        range1=[0 1]*2e-7; 
    coplot(i,xg,yg,p.*xg,ath,xlab,ylab,titolo,range1, range2,3)
  