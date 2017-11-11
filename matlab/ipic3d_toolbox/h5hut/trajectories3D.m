clear all
close all
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;

folder_name = '/Users/gianni/Dropbox/Science/san_diego/6-coil-Nov-2017/'
namefile = 'SixCoils-Fields';

Lx=70;
Ly=70;
Lz=70;

CoilD=20;

methodflag=2;

qom =1; is= 1; % 0 electrons, 1 ions

% cycle
i=4000

    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3);
    
    iz=round(Nz/2)
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
     smorz_field = 1;
    bx = hdf5read(fn, 'Step#0/Block/Bx/0')/smorz_field;
    by = hdf5read(fn, 'Step#0/Block/By/0')/smorz_field;
    bz = hdf5read(fn, 'Step#0/Block/Bz/0')/smorz_field;
    bx_ext = hdf5read(fn, 'Step#0/Block/Bx_ext/0');
    by_ext = hdf5read(fn, 'Step#0/Block/By_ext/0');
    bz_ext = hdf5read(fn, 'Step#0/Block/Bz_ext/0');
    bx=(bx+bx_ext);
    by=(by+by_ext);
    bz=(bz+bz_ext);
    
    
    smorz_field = 1e10;
    ex = hdf5read(fn, 'Step#0/Block/Ex/0')/smorz_field;
    ey = hdf5read(fn, 'Step#0/Block/Ey/0')/smorz_field;
    ez = hdf5read(fn, 'Step#0/Block/Ez/0')/smorz_field;

    
   
    
%     ex=permute(ex),[2 1]);
%     ey=permute(squeeze(ey(:,:,round(Nz/2))),[2 1]);
%     ez=permute(squeeze(ez(:,:,round(Nz/2))),[2 1]);
%     
%     bx=permute(squeeze(bx(:,:,round(Nz/2))),[2 1]);
%     by=permute(squeeze(by(:,:,round(Nz/2))),[2 1]);
%     bz=permute(squeeze(bz(:,:,round(Nz/2))),[2 1]);
    
    b = sqrt (bx.^2 +by.^2 + bz.^2);bmax=max(b(:));
    e = sqrt (ex.^2 +ey.^2 + ez.^2);
   
%     magnetic_damping=1;
%     
%      ii=b<mean(b(:))/10;
%      bx(ii)=bx(ii)*magnetic_damping;
%      by(ii)=by(ii)*magnetic_damping;
%      bz(ii)=bz(ii)*magnetic_damping;
%      b = sqrt(bx.^2 +by.^2 + bz.^2);
   
     rho = hdf5read(hinfo.GroupHierarchy.Groups.Groups.Groups(30).Datasets(1));

%      rho=permute(squeeze(rho(:,:,round(Nz/2))),[2 1])*electric_damping;
     
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    az=vecpot(xc,yc,squeeze(bx(:,:,iz)),squeeze(by(:,:,iz)));
    
    global ex ey ez bx by bz Lx Ly Lz qom Rout
    global contours dx dy dz lprint
    global xn yn zn 
    global randi1 randi2
    
    randi1=rand(30,1);
    randi2=rand(30,1);
    
    contours = 1 ;
    lprint=0;
    
    [xg,yg]=meshgrid(0:Nx-1,0:Ny-1);
    xg=xg/(Nx-1)*Lx;
    yg=yg/(Ny-1)*Ly;
    
    [xn, yn, zn] = ndgrid(0:Nx-1,0:Ny-1,0:Nz-1);
    xn=xn/(Nx-1)*Lx;
    yn=yn/(Ny-1)*Ly;
    zn=zn/(Nz-1)*Lz;
    
%     V=.001*((xn-Lx/2).^2/Lx^2+(yn-Ly/2).^2/Ly^2+(zn-Lz/2).^2/Lz^2);
%     [exes,eyes,ezes]=gradient(permute(V,[1 2 3]));
    
    exes=-(xn-Lx/2)/Lx*bmax/10000;
    eyes=-(yn-Ly/2)/Lx*bmax/10000;
    ezes=-(zn-Lz/2)/Lz*bmax/10000;

       ex=ex+exes;ey=ey+eyes;ez=ez+ezes;
    h=figure(1)
    set(h,'Position',[677 70 627 910])
    
    xlab='x';
    ylab='y'
    titolo=[ 'cycle=' num2str(i) '  B (color) Ath(contours)']
    
        bmax= max(b(:))
    range2=[0 bmax];
    range1=log([bmax/100 bmax]);
    figure(1)
     coplot(i,xg,yg,log(squeeze(b(:,:,iz))),squeeze(b(:,:,iz)),xlab,ylab,titolo,range1, range2)


    hold on

    
    figure(2)

    isosurface(permute(xn, [2 1 3]),permute(yn, [2 1 3]),permute(zn, [2 1 3]),permute(b, [2 1 3]),bmax/4)
hold on
axis equal
grid on

    %print(['frame_' it '.png'],'-dpng')
    %pause(.01)
    
    Rout=Lx*.4
    
    th = -pi:pi/10:pi;
    xunit = Lx/2+ Rout * cos(th) ;
    yunit = Rout * sin(th) + Ly/2;
    plot(xunit, yunit,'w');
    
    Npart=10;
    
    mean_t=0;
    traffic=0;
    tic
    for ip=1:Npart
        
        random=0
        if(random)
        % random
        xp=[Lx/2+ 0.5 Ly/2 Lz/2+0]+rand(1,3)*3;
        vmono = .01;
        costh = rand(1,1)*2-1; sinth =1-costh^2;
        fi = 2*pi *rand(1,1);
        vp = vmono * [sinth*cos(fi) costh sinth*sin(fi)] ;
        else
       % deterministic
        xp=[Lx/2+0.5 Ly/2 Lz/2+0]+ ip * ones(1,3)/Npart*3;
        vmono = .01;
        costh = ip/Npart*2-1; sinth =1-costh^2;
        fi = 2*pi*ip/Npart;
        vp = vmono * [sinth*cos(fi) costh sinth*sin(fi)] ;
        end
    opts=odeset('Events',@lostparticle3D); %,'OutputFcn',@odeplot);
     
    Tend=30000;
    dt=100; % this is used onyl for graphics, the actual time stpe is adaptive 
    % Slow matlab intrinsic
    %[t,y]=ode45(@newton,0:dt:Tend ,[xp vp],opts);
    % Fast implementation 
    [t,y]=ode45(@newton_interp3D,0:dt:Tend ,[xp vp],opts);
    
    xout=y(:,1);
    yout=y(:,2);
    zout=y(:,3);
    vxout=y(:,4);
    vyout=y(:,5);
    vzout=y(:,6);
    
    
    
    dx=diff(xout); dy=diff(yout);dz=diff(zout);
    traffic=traffic+sum(sqrt(dx.^2+dy.^2+dz.^2));
    
    kex=vxout.^2;
    key=vyout.^2;
    kez=vzout.^2;
    ke=kex+key+kez;
    mean_t =mean_t + t(end);
    figure(1)
    if(Tend==t(end))
        % particle remains confined
        
        % this colors line by energy
%         col = sum(y(:,4:6).^2,2);  % This is the color, vary with x in this case.
%         surface([r,r],[yout,yout],[zeros(size(r)),zeros(size(r))],[col,col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%         plot(r(1),yout(1),'go',r(end),yout(end),'kp')
figure(1)
        plot(xout,yout,xout(1),yout(1),'go',xout(end),yout(end),'kp','linew',2)
        figure(2)
        plot3(xout,yout,zout)
        hold on

    else
        % particle lost
        
%                 col = sum(y(:,4:6).^2,2);  % This is the color, vary with x in this case.
%         surface([r,r],[yout,yout],[zeros(size(r)),zeros(size(r))],[col,col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);
%         plot(r(1),yout(1),'go',r(end),yout(end),'kx')
figure(1)
        plot(xout,yout,xout(1),yout(1),'go',xout(end),yout(end),'kx','linew',2)
         figure(2)
        plot3(xout,yout,zout)
        hold on
        
    end    
    
    pause(.1)
    disp(['particles processed', num2str(ip)])
    disp(['residence time=', num2str(mean_t/ip)])
    disp(['crossings =', num2str(mean_t/ip*(vmono/CoilD))]) 
    disp(['mean traffic =', num2str(traffic/ip)]) 
    end
      toc 
      figure(1)
      print('-dpng','-r300',['nnoE' it 'figure1.png'])
            figure(2)
      print('-dpng','-r300',['nnoE' it 'figure2.png'])
    

