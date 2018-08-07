clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Science/san_diego/high-res-steady-state'
folder_name = '/Users/gianni/Downloads/pressure-anisotropy'
folder_name = '/Users/giovannilapenta/Dropbox/Science/ucla/controlloE/picket-fence with divE-removed'
folder_name = '/Users/giovannilapenta/Dropbox/Science/ucla/controlloE/picket-fencewithdivE'
namefile = 'PF4-Fields';


Lx=40;
Ly=30;
qom_0 = -64;
qom_1 = 1;


i = 64000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    ne = hdf5read(fn,'/Step#0/Block/rho_0/0/');
    jx = hdf5read(fn,'/Step#0/Block/Jx_0/0/');
    pxx = hdf5read(fn,'/Step#0/Block/Pxx_0/0/');
    pexx = (pxx - jx.*jx./(ne-1e-10))/qom_0;
    jy = hdf5read(fn,'/Step#0/Block/Jy_0/0/');
    pyy = hdf5read(fn,'/Step#0/Block/Pyy_0/0/');
    peyy = (pyy - jy.*jy./(ne-1e-10))/qom_0;
    jz = hdf5read(fn,'/Step#0/Block/Jz_0/0/');
    pzz = hdf5read(fn,'/Step#0/Block/Pzz_0/0/');
    pezz = (pzz - jz.*jz./(ne-1e-10))/qom_0;
    
   
    pe = pexx + peyy + pezz;
    
    ni = hdf5read(fn,'/Step#0/Block/rho_1/0/');
    jx = hdf5read(fn,'/Step#0/Block/Jx_1/0/');
    pxx = hdf5read(fn,'/Step#0/Block/Pxx_1/0/');
    pixx = (pxx - jx.*jx./(ni+1e-10))/qom_1;
    jy = hdf5read(fn,'/Step#0/Block/Jy_1/0/');
    pyy = hdf5read(fn,'/Step#0/Block/Pyy_1/0/');
    piyy = (pyy - jy.*jy./(ni+1e-10))/qom_1;
    jz = hdf5read(fn,'/Step#0/Block/Jz_1/0/');
    pzz = hdf5read(fn,'/Step#0/Block/Pzz_1/0/');
    pizz = (pzz - jz.*jz./(ni+1e-10))/qom_1;
    
    pion = pixx + piyy + pizz;
    
    p = pe + pion;
    
    vtheta=jz./(ni+1e-10);
    
    p=permute(squeeze(p(:,:,round(Nz/2))),[2 1]);
    pe=permute(squeeze(pe(:,:,round(Nz/2))),[2 1]);
    pion=permute(squeeze(pion(:,:,round(Nz/2))),[2 1]);
    pexx=permute(squeeze(pexx(:,:,round(Nz/2))),[2 1]);
    peyy=permute(squeeze(peyy(:,:,round(Nz/2))),[2 1]);
    pezz=permute(squeeze(pezz(:,:,round(Nz/2))),[2 1]);

    pixx=permute(squeeze(pixx(:,:,round(Nz/2))),[2 1]);
    piyy=permute(squeeze(piyy(:,:,round(Nz/2))),[2 1]);
    pizz=permute(squeeze(pizz(:,:,round(Nz/2))),[2 1]);
    
    ne=permute(squeeze(ne(:,:,round(Nz/2))),[2 1]);
    ni=permute(squeeze(ni(:,:,round(Nz/2))),[2 1]);
    
    vtheta=permute(squeeze(vtheta(:,:,round(Nz/2))),[2 1]);
    
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

     
    lde = sqrt(pe./4./3.14./(ne.^2+1e-10));
    rhoe = sqrt(pe./ne./qom_0./b.^2);
    
    global contours dx dy
    
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
    titolo=[ 'cycle=' num2str(i) '  \lambda_{De}(int) (color) A(contours)']
    

    range1=[0 .1]; 
    range2=[0 0];
 
    ptot=p+ben;
  coplot(i,xg,yg,lde,ath,xlab,ylab,titolo,range1, range2)
  titolo=[ 'cycle=' num2str(i) '  \rho_e(int) (color) A(contours)']
    figure(3)
    subplot(2,2,1)

  coplot(i,xg,yg,ne+ni,ath,xlab,ylab,'net \rho',[-.5 .5], range2)
  load cm_new
colormap(cm_kbwrk)
  subplot(2,2,2)
  coplot(i,xg,yg,(ne+ni).*xg,ath,xlab,ylab,'r \rho',[-20 20], range2)
colormap(cm_kbwrk)
  %subplot(2,2,4)
  %coplot(i,xg,yg,cumsum((ne+ni).*xg,2)*dx,ath,xlab,ylab,titolo,[0 0], range2)
subplot(2,2,3)
  coplot(i,xg,yg,cumsum((ne+ni).*xg,2)*dx./xg/4/pi,ath,xlab,ylab,'(\int r \rho dr)/4\pir',[-.1 .1], range2)
subplot(2,2,4)
  coplot(i,xg,yg,ex,ath,xlab,ylab,'Ex',[0 0], range2)
print -dpng 'withDiv.png'
 
  
return
  figure(2)
  subplot(3,1,1)
  iy=round(Ny/2);
  plot(xg(iy,:),ben(iy,:))
  title('radial-axis (x)')
  hold on 
   plot(xg(iy,:),p(iy,:))
      plot(xg(iy,:),ptot(iy,:))
      xlabel('r')
   
  subplot(3,1,2)
  ix=1;
  plot(yg(:,ix),ben(:,ix))
  title('central-axis (y)')
  hold on 
   plot(yg(:,ix),p(:,ix))
    plot(yg(:,ix),ptot(:,ix))
         xlabel('z')

  subplot(3,1,3)
plot(diag(xg(1+(Ny-1)/2:end,:)),diag(ben(1+(Ny-1)/2:end,:)))
  title('upper diagonal')
  hold on 
  plot(diag(xg(1+(Ny-1)/2:end,:)),diag(p(1+(Ny-1)/2:end,:)))
  plot(diag(xg(1+(Ny-1)/2:end,:)),diag(ptot(1+(Ny-1)/2:end,:)))
  legend('P_B','P','P_{tot}')
  ylim([0 4e-7])
        xlabel('r')
        
        
      
