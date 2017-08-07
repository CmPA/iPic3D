clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;

folder_name = '/shared/gianni/WB/base'
namefile = 'TwoCoils-Fields';



Lx=37.5;
Ly=75;

qom =10;

% for initial vacuum field
i=0
% for WB later field
i=35000
%i=4000
%i=00

    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    
    
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
    
    electric_damping=1;
    
    ex=permute(squeeze(ex(:,:,round(Nz/2))),[2 1])*electric_damping;
    ey=permute(squeeze(ey(:,:,round(Nz/2))),[2 1])*electric_damping;
    ez=permute(squeeze(ez(:,:,round(Nz/2))),[2 1])*electric_damping;
    
    bx=permute(squeeze(bx(:,:,round(Nz/2))),[2 1]);
    by=permute(squeeze(by(:,:,round(Nz/2))),[2 1]);
    bz=permute(squeeze(bz(:,:,round(Nz/2))),[2 1]);
    
    b = sqrt (bx.^2 +by.^2 + bz.^2);
    
   
    magnetic_damping=1;
    
     ii=b<mean(b(:))/10;
     bx(ii)=bx(ii)*magnetic_damping;
     by(ii)=by(ii)*magnetic_damping;
     bz(ii)=bz(ii)*magnetic_damping;
     b = sqrt(bx.^2 +by.^2 + bz.^2);
   
     rho = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(30).Datasets(1));

     rho=permute(squeeze(rho(:,:,round(Nz/2))),[2 1])*electric_damping;
     
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    ath=vecpot_cyl(xc,yc,bx,by);
    
    global ex ey ez bx by bz xg yg  Lx Ly qom Rout
    
    [xg,yg]=meshgrid(0:Nx-1,0:Ny-1);
    xg=xg/Nx*Lx;
    yg=yg/Ny*Ly;
    
    h=figure(1)
    set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    titolo=[ 'cycle=' num2str(i) '  B (color) Ath(contours)']
    
    
    range1=[-10 -3]; 
    range2=[0 0];
    %coplot(i,xg,yg,log(b+1e-10),b,xlab,ylab,titolo,range1, range2)
    range1=[0 0];
     coplot(i,xg,yg,log(b),ath,xlab,ylab,titolo,range1, range2)

    hold on
    
    %print(['frame_' it '.png'],'-dpng')
    %pause(.01)
    
    Rout=20
    
    th = -pi/2:pi/50:pi/2;
    xunit = Rout * cos(th) ;
    yunit = Rout * sin(th) + Ly/2;
    plot(xunit, yunit);
    
    Npart=100;
    
    mean_t=0;
    traffic=0;
    
    for ip=1:Npart
        
        xp=[0.5 Ly/2 0]+rand(1,3)*3;

        vmono = .01;
        costh = rand(1,1)*2-1; sinth =1-costh^2;
        fi = 2*pi *rand(1,1);
        vp = vmono * [sinth*cos(fi) costh sinth*sin(fi)] ;
    
    opts=odeset('Events',@lostparticle); %,'OutputFcn',@odeplot);
     
    Tend=30000;
    dt=1; % this is used onyl for graphics, the actual time stpe is adaptive 
    [t,y]=ode45(@newton,0:dt:Tend ,[xp vp],opts);
    
    xout=y(:,1);
    yout=y(:,2);
    zout=y(:,3);
    r = sqrt(xout.^2+zout.^2);
    dx=diff(xout); dy=diff(yout);dz=diff(zout);
    traffic=traffic+sum(sqrt(dx.^2+dy.^2+dz.^2));
    
    mean_t =mean_t + t(end);
    if(Tend==t(end))
        % particle remains confined
        
        % this colors line by energy
%         col = sum(y(:,4:6).^2,2);  % This is the color, vary with x in this case.
%         surface([r,r],[yout,yout],[zeros(size(r)),zeros(size(r))],[col,col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%         plot(r(1),yout(1),'go',r(end),yout(end),'kp')
        plot(r,yout,r(1),yout(1),'go',r(end),yout(end),'kp','linew',2)

    else
        % particle lost
        
%                 col = sum(y(:,4:6).^2,2);  % This is the color, vary with x in this case.
%         surface([r,r],[yout,yout],[zeros(size(r)),zeros(size(r))],[col,col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%         plot(r(1),yout(1),'go',r(end),yout(end),'kx')
        plot(r,yout,r(1),yout(1),'go',r(end),yout(end),'kx','linew',2)
        
    end    
    pause(.1)
%     if (sqrt(r(end).^2+(yout(end)-Ly/2).^2)<Rout-Rout/100) 
%         return 
%     end
    disp(['particles processed', num2str(ip)])
    disp(['residence time=', num2str(mean_t/ip)])
    disp(['mean traffic =', num2str(traffic/ip)]) 
    end
       
    

