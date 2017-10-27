clear all
close all
!rm *.png
!rm *.eps
clc
addpath(genpath('../../ipic3d_toolbox'))
folder_name = pwd;
folder_name = '/Users/gianni/Dropbox/Science/san_diego/high-res-steady-state'
namefile = 'TwoCoils2D-Partcl';


Lx=25;
Ly=50;


i = 72000


    it=sprintf('%06.0f',i);
        
    fn=[folder_name,'/',namefile,'_',it,'.h5'];
    
    hinfo=hdf5info(fn);
    
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups(:).Name
    
    is=1;
    
    q = hdf5read(fn,['/Step#0/q_' num2str(is) '/']);
    x = hdf5read(fn,['/Step#0/x_' num2str(is) '/']);
    y = hdf5read(fn,['/Step#0/y_' num2str(is) '/']);
    u = hdf5read(fn,['/Step#0/u_' num2str(is) '/']);
    v = hdf5read(fn,['/Step#0/v_' num2str(is) '/']);
    w = hdf5read(fn,['/Step#0/w_' num2str(is) '/']);
    
    %
    %   Extract subset
    %
        radius = .2
        volgorde= 0
    for posx=0:radius:Lx
        posy=posx+Ly/2
    %for posy=Ly/2:radius:Ly
    %    posx=radius;
    ii=abs(y-posy)<radius & abs(x-posx)<radius;
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
        vmax=.02;
    else
        vmax=.02/4;
    end    
    vmin= -vmax;
    ndiv=50;
    vdf_sp=spaziofasi3D(usub,vsub,wsub,qsub,vmin,vmax,ndiv);
        vdf_sp=vdf_sp./sum(vdf_sp(:));
        
    global color_choice symmetric_color titolo square labelT
        
    Nsm= 2
    square =1
    color_choice = 3
    symmetric_color = 0 
    labelT=''
    titolo=['x=' num2str(posx) '   y=' num2str(posy) '   Npart=' num2str(Nsub)];
    immagine_dir([vmin vmax],[vmin vmax],log10(1e-10+squeeze(sum(vdf_sp,2))), ...
             ['vdfXZ_' 'species_' num2str(is) '_' num2str(volgorde)], ...
             [0 5e-3]*0,Nsm,titolo,0,1,'v_x','v_z','vdf');
             immagine_dir([vmin vmax],[vmin vmax],log10(1e-10+squeeze(sum(vdf_sp,3))), ...
             ['vdfXY_' 'species_' num2str(is) '_' num2str(volgorde)], ...
             [0 5e-3]*0,Nsm,titolo,0,2,'v_x','v_y','vdf');
    pause(.1)
    volgorde = volgorde +1; 
    end
    
    end     
    