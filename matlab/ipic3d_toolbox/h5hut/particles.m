clear all
close all
!rm *.png
!rm *.eps
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Science/codes/build_ecsim/data/'
namefile = 'MHDUCLA-Partcl';
namefile_field = 'MHDUCLA-Fields';

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange 
Xgsmrange= [-38.2 -7];
Zgsmrange= [-11.25 4.35];
Ygsmrange= [-7.32 13.48];
Ly=50
Lx=200
Lz=1

qom=[-256,1];

volgorde=0
for i = 1000:1000:7000
volgorde=volgorde+1;

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
    
    
    fn_field=[folder_name,namefile_field,'_',it,'.h5'];

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
        xc=linspace(0, Lx, Nx);
    yc=linspace(0, Ly, Ny);
    %AAz=vecpot(xc,yc,bx',by');AAz=AAz';

    iz=floor(Nz/2)
    AAz=vecpot(xc,yc,bx(:,:,iz),by(:,:,iz));AAz=AAz';
    end
    %
    %   Extract subset
    %
        mp=938e6
        me=mp/256;
        en_level = 3*1e3/me


    e = (u.^2+v.^2+w.^2)/2;
    %for posy=Ly/2:radius:Ly
    %    posx=radius;
    ii=e>en_level & abs(x-Lx/2)<Lx/2-Lx/10;
    sum(ii)
    qsub=q(ii);
    xsub=x(ii);
    ysub=y(ii);
    usub=u(ii); 
    vsub=v(ii);
    wsub=w(ii);
    Nsub=max(size(xsub))
    if Nsub>100
    %plot(usub,vsub,'.','MarkerSize',[1])
    %title(num2str(size(xsub)))
    if(is==0)
        vmax=.1;
    else
        vmax=.11/4;
    end    
    vmin= -vmax;
    ndiv=50;
    vdf_sp=spaziofasi3D(-usub,wsub,vsub,qsub,vmin,vmax,ndiv);
        %vdf_sp=vdf_sp./sum(vdf_sp(:));
        
    global color_choice symmetric_color titolo square labelT
      
    figure(1)
    Nsm= 2
    square =1
    color_choice = 0
    symmetric_color = 0 
    labelT=''
    titolo=['e(keV)>' num2str(me*en_level/1e3) '   Npart=' num2str(Nsub)];
    immagine_dir([vmin vmax],[vmin vmax],log10(1e-10+squeeze(sum(vdf_sp,2))), ...
             ['vdfXZ_' 'species_' num2str(is) '_' num2str(volgorde)], ...
             [0 5e-3]*0,Nsm,titolo,0,1,'v_x/c','v_z/c','vdf');
    immagine_dir([vmin vmax],[vmin vmax],log10(1e-10+squeeze(sum(vdf_sp,3))), ...
             ['vdfXY_' 'species_' num2str(is) '_' num2str(volgorde)], ...
             [0 5e-3]*0,Nsm,titolo,0,1,'v_x/c','v_y/c','vdf');
    figure(2)
    contourf(gsmx(xc),gsmy2z(yc),AAz)
    hold on
    plot(gsmx(xsub),gsmy2z(ysub),'w.')
    title(['e(keV)>' num2str(me*en_level/1e3) ])
    axis equal
    set(gca,'xdir','reverse','TickDir','out')
    %xlim(gsmx([0 Lx]))
    %ylim(gsmy2z([0 Ly]))
    xlabel('x/R_E', 'fontsize',[14])
    ylabel('z/R_E', 'fontsize',[14])
    set(gca,'fontsize',[14])
    print('-dpng',['scatter' num2str(volgorde)])

    figure(3)
    [v,ev]=hist(e*me/1e3,100);
    Np=max(size(e));
    loglog(ev,v/Np)
        xlim([.1 10])
    xlabel('E(keV)', 'fontsize',[14])
    ylabel('Fraction of Particles', 'fontsize',[14])
    set(gca,'fontsize',[14])
    print('-dpng',['part' num2str(volgorde)])

    end 
    
end  
    
%figure(3)
%legend(num2str(4000:7000))