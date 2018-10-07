clear all
close all
!rm *.png
!rm *.eps
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/shared/gianni/emc2paper/PF-g3-ss0-qom64-run/PF-g3-ss0-qom64-damp-re80k'
namefile = 'PF4-Partcl';
namefile_field = 'PF4-Fields';

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange 

Lx=45;
Ly=30;


qom=[-64,1];


i = 104000

    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];
    
    hinfo=hdf5info(fn);
    
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups(:).Name
    
    is=0;
    
    q = hdf5read(fn,['/Step#0/q_' num2str(is) '/']);
    x = hdf5read(fn,['/Step#0/x_' num2str(is) '/']);
    y = hdf5read(fn,['/Step#0/y_' num2str(is) '/']);
    u = hdf5read(fn,['/Step#0/u_' num2str(is) '/']);
    v = hdf5read(fn,['/Step#0/v_' num2str(is) '/']);
    w = hdf5read(fn,['/Step#0/w_' num2str(is) '/']);
    
    
    fn_field=[folder_name,'/',namefile_field,'_',it,'.h5'];

    hinfo=hdf5info(fn_field);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3);
    old=0
    if(old)
    bx = hdf5read(fn_field,'/Step#0/Block/Btx/0/');
    by = hdf5read(fn_field,'/Step#0/Block/Bty/0/');
    bz = hdf5read(fn_field,'/Step#0/Block/Btz/0/');
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    %AAz=vecpot(xc,yc,bx',by');AAz=AAz';

    iz=floor(Nz/2)
    AAz=vecpot(xc,yc,bx(:,:,iz),by(:,:,iz));AAz=AAz';
       else
    bx = hdf5read(fn_field,'/Step#0/Block/Bx/0/');
    by = hdf5read(fn_field,'/Step#0/Block/By/0/');
    bz = hdf5read(fn_field,'/Step#0/Block/Bz/0/');
    bx_ext = hdf5read(fn_field,'/Step#0/Block/Bx_ext/0/');
    by_ext = hdf5read(fn_field,'/Step#0/Block/By_ext/0/');
    bz_ext = hdf5read(fn_field,'/Step#0/Block/Bz_ext/0/');
    bx=bx+bx_ext;
    by=by+by_ext;
    bz=bz+bz_ext;
    
    bx=permute(squeeze(bx(:,:,round(Nz/2))),[2 1]);
    by=permute(squeeze(by(:,:,round(Nz/2))),[2 1]);
    bz=permute(squeeze(bz(:,:,round(Nz/2))),[2 1]);
    
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    %AAz=vecpot(xc,yc,bx',by');AAz=AAz';

    iz=floor(Nz/2)
    xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    ath=vecpot_cyl(xc,yc,bx,by);
    b=sqrt(bx.^2+by.^2+bz.^2);
    end
    %
    %   Extract subset
    %
radius = .2
        
        posy=Ly/4


    ii=abs(y-posy)<radius ;
    sum(ii)
    qsub=q(ii);
    xsub=x(ii);
    ysub=y(ii);
    usub=u(ii); 
    vsub=v(ii);
    wsub=w(ii);
    Nsub=max(size(xsub))

    if(is==0)
        vmax=.5;
    else
        vmax=.5/sqrt(64);
    end    
    vmin= -vmax;
    ndiv=100;
    [totnum,vdf_sp,xrange,urange]=spaziofasi2(xsub,usub,qsub,0,0,Lx,vmin,vmax,ndiv);
        %vdf_sp=vdf_sp./sum(vdf_sp(:));
        
    global color_choice symmetric_color titolo square labelT
      
   close all
    Nsm= 2
    square =1
    color_choice = 0
    symmetric_color = 0 
    labelT=''
    titolo=['Npart=' num2str(Nsub)];
    immagine_dir([0 Lx],[vmin vmax],log10(vdf_sp'), ...
             ['vdfXZ_' 'species_' num2str(is) '_' ], ...
             [0 5e-3]*0,Nsm,titolo,0,1,[1 1],'v_x/c','v_z/c','vdf');

    xlabel('x', 'fontsize',[14])
    ylabel('y', 'fontsize',[14])
    set(gca,'fontsize',[14])
    hold off
    print('-dpng',['spaziofasi' sprintf('%06.0f',is)])

  

