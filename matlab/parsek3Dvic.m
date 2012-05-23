%  parsek3D       -  Assemble hdf5 data for Parsek3D 
%			  
%   parsek2D reads data saved in hdf5 from all processors and assemble them in the 
%   correct topology. The user needs to specify the following variables:
%   results_dir = path where the *.hdf file are saved;
%   variable list = string variable with quantities that will be read. Possible tags are:
%   B, E -> magnetic and electric field
%   A, phi -> potentials
%   J, rho -> current and number densities
%   pressure
%   x, v, ID -> position , velocity, ID
%   k_energy -> kinetic energy
%   E_energy, B_energy -> electric and magnetic energy
%

clear all
close all

global Lx Ly Lz Nx Ny Nz Dx Dy Dz Dt Th Ncycles Ns c Smooth
global PfaceXright PfaceXleft PfaceYright PfaceYleft PfaceZright PfaceZleft
global PHIfaceXright PHIfaceXleft PHIfaceYright PHIfaceYleft PHIfaceZright PHIfaceZleft 
global EMfaceXright EMfaceXleft EMfaceYright EMfaceYleft EMfaceZright EMfaceZleft

global Np Npcelx Npcely Npcelz NpMax qom
global uth vth wth u_drift v_drift w_drift
global XLEN YLEN ZLEN Nprocs periodicX periodicY periodicZ


results_dir='/l/users/s0198231/parsek3D/' ; % directory for results

variable_list='B E rho J';  % list of variable you want to load :
                                                  %  B, E -> magnetic and electric field
                                                  %  A, phi -> potentials
                                                  %  J, rho -> current and number densities
                                                  %  pressure
                                                  %  x, v, ID, q -> position , velocity, ID, q
                                                  %  k_energy -> kinetic energy
                                                  %  E_energy, B_energy -> electric and magnetic energy
                                                 

variable_list=' B E J k_energy e_energy b_energy ';
variable_list='k_energy e_energy b_energy ';
variable_list=' B E rho J';

setting_name=[results_dir 'settings.hdf']; % file with collective settings
processor_name=[results_dir 'proc']; % file with data for each processor (i.e. proc0, proc1, ...)


read_parsek_settings3D(setting_name);

nxc = Nx/XLEN + 2;
nyc = Ny/YLEN + 2;
nzc = Nz/ZLEN + 2;
nxn = nxc + 1;
nyn = nyc + 1;
nzz = nzc + 1;


for i=1:Nprocs
 cart=hdf5read([processor_name,num2str(i-1),'.hdf'],'/topology/cartesian_coord');
 MAP(cart(1)+1,cart(2)+1,cart(3)+1)=hdf5read([processor_name,num2str(i-1),'.hdf'],'/topology/cartesian_rank')+1;
end

for PROC=1:Nprocs
info=hdf5info([processor_name,num2str(PROC-1),'.hdf']);
nGroups=size(info.GroupHierarchy.Groups,2);
 
% read all the groups	   
for i=1:nGroups

if strcmp(info.GroupHierarchy.Groups(i).Name,'/energy')  
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);

    for ii=1:nnGroups

      if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/energy/electric') & strfind (variable_list, 'e_energy')
        E_nrg_time=uint16([]);
        ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
        for time=1:ntime
            E_nrg_slab(time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            E_nrg_time=[E_nrg_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/energy/electric/cycle_',''))];
        end
        [E_nrg_time,E_nrg_index]=intersect(E_nrg_time,sort(E_nrg_time),'rows');
      end
      
     if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/energy/magnetic') & strfind (variable_list, 'b_energy')
        B_nrg_time=uint16([]);
        ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
        for time=1:ntime
            B_nrg_slab(time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            B_nrg_time=[B_nrg_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/energy/magnetic/cycle_',''))];
        end
        [B_nrg_time,B_nrg_index]=intersect(B_nrg_time,sort(B_nrg_time),'rows');
     end

              
      if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/energy/kinetic') & strfind (variable_list, 'k_energy')
        nnnGroups=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups,2);
     
          for iii=1:nnnGroups
            for nspec=1:Ns
                if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/energy/kinetic/species_',num2str(nspec-1)])
                 k_nrg_time=uint16([]);
                 ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                    for time=1:ntime
                     eval(['k_nrg_slab' num2str(nspec-1) '(time,PROC)= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);         
                     k_nrg_time=[k_nrg_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/energy/kinetic/species_',num2str(nspec-1),'/cycle_'], ''))];
                    end
                    [k_nrg_time,k_nrg_index]=intersect(k_nrg_time,sort(k_nrg_time),'rows');
                end
            end
          end

      end
    end
   
    
end
% fields    
if strcmp(info.GroupHierarchy.Groups(i).Name,'/fields') 
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);
    
    for ii=1:nnGroups

       if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/Bx') & strfind (variable_list, 'B')
            Bx_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Bx_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Bx_time=[Bx_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/Bx/cycle_',''))];
           end
           [Bx_time,Bx_index]=intersect(Bx_time,sort(Bx_time),'rows');
       end

       if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/By') & strfind (variable_list, 'B')
            By_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            By_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            By_time=[By_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/By/cycle_',''))];
           end
           [By_time,By_index]=intersect(By_time,sort(By_time),'rows');
       end

       if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/Bz') & strfind (variable_list, 'B')
            Bz_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Bz_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Bz_time=[Bz_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/Bz/cycle_',''))];
           end
           [Bz_time,Bz_index]=intersect(Bz_time,sort(Bz_time),'rows');
       end

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/Ex') & strfind (variable_list, 'E')
            Ex_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Ex_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Ex_time=[Ex_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/Ex/cycle_',''))];
           end
           [Ex_time,Ex_index]=intersect(Ex_time,sort(Ex_time),'rows');
       end

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/Ey') & strfind (variable_list, 'E')
            Ey_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Ey_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Ey_time=[Ey_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/Ey/cycle_',''))];
           end
           [Ey_time,Ey_index]=intersect(Ey_time,sort(Ey_time),'rows');
       end
   
        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/fields/Ez') & strfind (variable_list, 'E')
            Ez_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Ez_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Ez_time=[Ez_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/fields/Ez/cycle_',''))];
           end
           [Ez_time,Ez_index]=intersect(Ez_time,sort(Ez_time),'rows');
       end

    end
end

if strcmp(info.GroupHierarchy.Groups(i).Name,'/potentials')
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);
    
    for ii=1:nnGroups
        % REMOVED THE VECTOR POTENTIAL by stef


        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/potentials/phi') & strfind (variable_list, 'phi')
            phi_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            phi_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            phi_time=[phi_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,'/potentials/phi/cycle_',''))];
           end
           [phi_time,phi_index]=intersect(phi_time,sort(phi_time),'rows');
        end

    end
end

if strcmp(info.GroupHierarchy.Groups(i).Name,'/moments')
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);
    
    for ii=1:nnGroups

        for nspec=1:Ns
              
            if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,['/moments/species_',num2str(nspec-1)])
                
                 nnnGroups=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups,2);

                 for iii=1:nnnGroups
                 
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/Jx']) & strfind (variable_list, 'J')
                    Jxs_time=uint16([]);
                    ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                    for time=1:ntime
                        eval(['Jxs_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                        Jxs_time=[Jxs_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/moments/species_',num2str(nspec-1),'/Jx','/cycle_'], ''))];
                    end
                    [Jxs_time,Jxs_index]=intersect(Jxs_time,sort(Jxs_time),'rows');
                  end

                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/Jy']) & strfind (variable_list, 'J')
                      Jys_time=uint16([]);
                      ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                        eval(['Jys_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                      Jys_time=[Jys_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/moments/species_',num2str(nspec-1),'/Jy','/cycle_'], ''))];
                     end
                     [Jys_time,Jys_index]=intersect(Jys_time,sort(Jys_time),'rows');
                  end
                                    
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/Jz']) & strfind (variable_list, 'J')
                    Jzs_time=uint16([]);
                      ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                        eval(['Jzs_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                     Jzs_time=[Jzs_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/moments/species_',num2str(nspec-1),'/Jz','/cycle_'], ''))];
                     end
                     [Jzs_time,Jzs_index]=intersect(Jzs_time,sort(Jzs_time),'rows');
                  end
                  
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/rho']) & strfind (variable_list, 'rho')
                      rhos_time=uint16([]);
                      ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                        eval(['rhos_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                        rhos_time=[rhos_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/moments/species_',num2str(nspec-1),'/rho','/cycle_'], ''))];
                     end
                     [rhos_time,rhos_index]=intersect(rhos_time,sort(rhos_time),'rows');
                  end

                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pXX']) & strfind (variable_list, 'pressure')
                     p_time=uint16([]);
                      ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                         eval(['pxx_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                         p_time=[p_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/moments/species_',num2str(nspec-1),'/pXX','/cycle_'], ''))];
                     end
                     [p_time,p_index]=intersect(p_time,sort(p_time),'rows');
                  end
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pXY']) & strfind (variable_list, 'pressure')
                    ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                         eval(['pxy_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                     end
                  end
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pXZ']) & strfind (variable_list, 'pressure')
                    ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                         eval(['pxz_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                     end
                  end
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pYY']) & strfind (variable_list, 'pressure')
                    ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                         eval(['pyy_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                     end
                  end
                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pYZ']) & strfind (variable_list, 'pressure')
                      ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                      for time=1:ntime
                          eval(['pyz_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                      end
                  end

                  if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/moments/species_',num2str(nspec-1),'/pZZ']) & strfind (variable_list, 'pressure')
                    ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                     for time=1:ntime
                         eval(['pzz_slab' num2str(nspec-1) '(:,:,:,time,PROC) = hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                     end
                  end
                  
                 end
            end
        end
        
                

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/moments/Jx') & strfind (variable_list, 'J')
            Jx_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Jx_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Jx_time=[Jx_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,['/moments/Jx/cycle_'], ''))];
           end
            [Jx_time,Jx_index]=intersect(Jx_time,sort(Jx_time),'rows');
        end

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/moments/Jy') & strfind (variable_list, 'J')
            Jy_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Jy_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Jy_time=[Jy_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,['/moments/Jy/cycle_'], ''))];
           end
           [Jy_time,Jy_index]=intersect(Jy_time,sort(Jy_time),'rows');
        end

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/moments/Jz') & strfind (variable_list, 'J')
            Jz_time=uint16([]);
            ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            Jz_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            Jz_time=[Jz_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,['/moments/Jz/cycle_'], ''))];
           end
           [Jz_time,Jz_index]=intersect(Jz_time,sort(Jz_time),'rows');
        end

        if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,'/moments/rho') & strfind (variable_list, 'rho')
           rho_time=uint16([]);
           ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Datasets,2);
           for time=1:ntime
            rho_slab(:,:,:,time,PROC)=hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time));
            rho_time=[rho_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Datasets(time).Name,['/moments/rho/cycle_'], ''))];
           end
           [rho_time,rho_index]=intersect(rho_time,sort(rho_time),'rows');
        end

    end
  end



if strcmp(info.GroupHierarchy.Groups(i).Name,'/particles')
    nnGroups=size(info.GroupHierarchy.Groups(i).Groups,2);

    for ii=1:nnGroups
  %nspec=ii;
        for nspec=1:Ns
          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Name,['/particles/species_',num2str(nspec-1)])
              nnnGroups=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups,2);
              
             for iii=1:nnnGroups
                 
                         
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/x']) & strfind (variable_list, 'x') 
                                x_time=uint16([]);
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                eval(['x_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);
                                x_time=[x_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/x','/cycle_'], ''))];
                                end
                                [x_time,x_index]=intersect(x_time,sort(x_time),'rows');
                          end
                          
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/y']) & strfind (variable_list, 'x')
                                y_time=uint16([]);   
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                    eval(['y_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);                  
                                y_time=[y_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/y','/cycle_'], ''))];
                                end
                                [y_time,y_index]=intersect(y_time,sort(y_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/z']) & strfind (variable_list, 'x')
                                z_time=uint16([]);
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                    eval(['z_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);                 
                                    z_time=[z_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/z','/cycle_'], ''))];
                                end
                                [z_time,z_index]=intersect(z_time,sort(z_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/u']) & strfind (variable_list, 'v')
                                u_time=uint16([]); 
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                     eval(['u_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);                  
                                u_time=[u_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/u','/cycle_'], ''))];
                                end
                                [u_time,u_index]=intersect(u_time,sort(u_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/v']) & strfind (variable_list, 'v')
                                v_time=uint16([]);
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                     eval(['v_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);         
                                 v_time=[v_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/v','/cycle_'], ''))];
                                end
                                [v_time,v_index]=intersect(v_time,sort(v_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/w']) & strfind (variable_list, 'v')
                                w_time=uint16([]);
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                     eval(['w_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);     
                                 w_time=[w_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/w','/cycle_'], ''))];
                                end
                                [w_time,w_index]=intersect(w_time,sort(w_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/q']) & strfind (variable_list, 'q')
                                q_time=uint16([]);
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                     eval(['q_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);     
                                 q_time=[q_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/q','/cycle_'], ''))];
                                end
                                [q_time,q_index]=intersect(q_time,sort(q_time),'rows');
                          end
                          if strcmp(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Name,['/particles/species_',num2str(nspec-1),'/ID']) & strfind (variable_list, 'ID')
                                ID_time=uint16([]);  
                                ntime=size(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets,2);
                                for time=1:ntime
                                      eval(['ID_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) '= hdf5read(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time));']);              
                                ID_time=[ID_time;str2double(regexprep(info.GroupHierarchy.Groups(i).Groups(ii).Groups(iii).Datasets(time).Name,['/particles/species_',num2str(nspec-1),'/ID','/cycle_'], ''))];
                                end
                               [ID_time,ID_index]=intersect(ID_time,sort(ID_time),'rows');
                          end
                          
                     
                      
                  
                  
                                    
                  
              end
          end
        end
        
    end
end
    
end

end

clear info nGroups nnGroups nnnGroups cart
%%%%%%%%%%%%%%%%%%%%% end of reading from .hdf file

%% assemble pieces..
%return

if exist('E_nrg_slab','var')
    E_energy=sum(E_nrg_slab,2);
    clear E_nrg_slab 

 for h=1:size(E_nrg_index)
    EE_energy(h)=E_energy(E_nrg_index(h));
 end
E_energy=EE_energy;
clear EE_energy
end

if exist('B_nrg_slab','var')
    B_energy=sum(B_nrg_slab,2);
    clear B_nrg_slab 

 for h=1:size(B_nrg_index)
    BB_energy(h)=B_energy(B_nrg_index(h));
 end
B_energy=BB_energy;
clear BB_energy
end
%Bx
if exist('Bx_slab','var')
    Bx_slab=(Bx_slab(1:end-1,:,:,:,:)+Bx_slab(2:end,:,:,:,:))*0.5;
    Bx_slab=(Bx_slab(:,1:end-1,:,:,:)+Bx_slab(:,2:end,:,:,:))*0.5;
	Bx_slab=(Bx_slab(:,:,1:end-1,:,:)+Bx_slab(:,:,2:end,:,:))*0.5;
    Bx=[];Bx_row=[]; Bx_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Bx_row=cat(2,Bx_row,Bx_slab(:,:,:,:,MAP(i,j,k)));
		  end
          Bx_col=cat(1,Bx_col,Bx_row);
          Bx_row=[];
          end
		  Bx=cat(3,Bx,Bx_col);
		  Bx_col=[];
		end
		 
    clear Bx_row Bx_col Bx_slab

 for h=1:size(Bx_index)
    BBx(:,:,:,h)=Bx(:,:,:,Bx_index(h));
 end
Bx=BBx; clear BBx Bx_index
Bx=squeeze(Bx);
end %endif

%By
if exist('By_slab','var')
    By_slab=(By_slab(1:end-1,:,:,:,:)+By_slab(2:end,:,:,:,:))*0.5;
    By_slab=(By_slab(:,1:end-1,:,:,:)+By_slab(:,2:end,:,:,:))*0.5;
	By_slab=(By_slab(:,:,1:end-1,:,:)+By_slab(:,:,2:end,:,:))*0.5;
    By=[];By_row=[];By_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                By_row=cat(2,By_row,By_slab(:,:,:,:,MAP(i,j,k)));
        end
        By_col=cat(1,By_col,By_row);
        By_row=[];
      end
	  By=cat(3,By,By_col);
	  By_col=[];  
	 end
clear By_row By_col By_slab

for h=1:size(By_index)
    BBy(:,:,:,h)=By(:,:,:,By_index(h));
end
By=BBy; clear BBy By_index
By=squeeze(By);
end
%Bz
if exist('Bz_slab','var')
    Bz_slab=(Bz_slab(1:end-1,:,:,:,:)+Bz_slab(2:end,:,:,:,:))*0.5;
    Bz_slab=(Bz_slab(:,1:end-1,:,:,:)+Bz_slab(:,2:end,:,:,:))*0.5;
	Bz_slab=(Bz_slab(:,:,1:end-1,:,:)+Bz_slab(:,:,2:end,:,:))*0.5;
    Bz=[];Bz_row=[]; Bz_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Bz_row=cat(2,Bz_row,Bz_slab(:,:,:,:,MAP(i,j,k)));
        end
         Bz_col=cat(1,Bz_col,Bz_row);
        Bz_row=[];
      end
	  Bz=cat(3,Bz,Bz_col);
	  Bz_col=[];  
	 end
clear Bz_row Bz_col Bz_slab

for h=1:size(Bz_index)
    BBz(:,:,:,h)=Bz(:,:,:,Bz_index(h));
end
Bz=BBz; clear BBz Bz_index 
Bz=squeeze(Bz);
end
%Ex
if exist('Ex_slab','var')
    Ex_slab=(Ex_slab(1:end-1,:,:,:,:)+Ex_slab(2:end,:,:,:,:))*0.5;
    Ex_slab=(Ex_slab(:,1:end-1,:,:,:)+Ex_slab(:,2:end,:,:,:))*0.5;
	Ex_slab=(Ex_slab(:,:,1:end-1,:,:)+Ex_slab(:,:,2:end,:,:))*0.5;
    Ex=[];Ex_row=[];Ex_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Ex_row=cat(2,Ex_row,Ex_slab(:,:,:,:,MAP(i,j,k)));
        end
        Ex_col=cat(1,Ex_col,Ex_row);
        Ex_row=[];
      end
	  Ex=cat(3,Ex,Ex_col);
	  Ex_col=[];  
	 end
clear Ex_row Ex_col  Ex_slab
for h=1:size(Ex_index)
    BEx(:,:,:,h)=Ex(:,:,:,Ex_index(h));
end
Ex=BEx; clear BEx Ex_index 
Ex=squeeze(Ex);
end
%Ey
if exist('Ey_slab','var')
    Ey_slab=(Ey_slab(1:end-1,:,:,:,:)+Ey_slab(2:end,:,:,:,:))*0.5;
    Ey_slab=(Ey_slab(:,1:end-1,:,:,:)+Ey_slab(:,2:end,:,:,:))*0.5;
	Ey_slab=(Ey_slab(:,:,1:end-1,:,:)+Ey_slab(:,:,2:end,:,:))*0.5;
    Ey=[];Ey_row=[]; Ey_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Ey_row=cat(2,Ey_row,Ey_slab(:,:,:,:,MAP(i,j,k)));
        end
        Ey_col=cat(1,Ey_col,Ey_row);
        Ey_row=[];
      end
	  Ey=cat(3,Ey,Ey_col);
	  Ey_col=[];  
	 end
clear Ey_row Ey_col Ey_slab
for h=1:size(Ey_index)
    BEy(:,:,:,h)=Ey(:,:,:,Ey_index(h));
end
Ey=BEy; clear BEy Ey_index 
Ey=squeeze(Ey);
end

if exist('Ez_slab','var')
    Ez_slab=(Ez_slab(1:end-1,:,:,:,:)+Ez_slab(2:end,:,:,:,:))*0.5;
    Ez_slab=(Ez_slab(:,1:end-1,:,:,:)+Ez_slab(:,2:end,:,:,:))*0.5;
	Ez_slab=(Ez_slab(:,:,1:end-1,:,:)+Ez_slab(:,:,2:end,:,:))*0.5;
    Ez=[];Ez_row=[]; Ez_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Ez_row=cat(2,Ez_row,Ez_slab(:,:,:,:,MAP(i,j,k)));
        end
        Ez_col=cat(1,Ez_col,Ez_row);
        Ez_row=[];
      end
	  Ez=cat(3,Ez,Ez_col);
	  Ez_col=[];  
	 end
    clear Ez_row Ez_col Ez_slab
    for h=1:size(Ez_index)
        BEz(:,:,:,h)=Ez(:,:,:,Ez_index(h));
    end
    Ez=BEz; clear BEz Ez_index
    Ez=squeeze(Ez);
end

% REMOVED the vector POTENTIAL by stef
% CALCULATED FORM POST PROCESSINg

% PHI CALCULATED on centres: averaging not need
if exist('phi_slab','var')
    phi=[];phi_row=[]; phi_col = [];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                phi_row=cat(2,phi_row,phi_slab(:,:,:,:,MAP(i,j,k)));
        end
       phi_col=cat(1,phi_col,phi_row);
	   phi_row=[];
      end
	  phi=cat(3,phi,phi_col);
	  phi_col=[];  
	 end
clear phi_row phi_col phi_slab

for h=1:size(phi_index)
    Bphi(:,:,:,h)=phi(:,:,:,phi_index(h));
end
phi=Bphi; clear Bphi phi_index 
phi=squeeze(phi);
end

%Jx
if exist('Jx_slab','var')
    Jx_slab=(Jx_slab(1:end-1,:,:,:,:)+Jx_slab(2:end,:,:,:,:))*0.5;
    Jx_slab=(Jx_slab(:,1:end-1,:,:,:)+Jx_slab(:,2:end,:,:,:))*0.5;
    Jx=[];Jx_row=[];Jx_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Jx_row=cat(2,Jx_row,Jx_slab(:,:,:,:,MAP(i,j,k)));
        end
       Jx_col=cat(1,Jx_col,Jx_row);
	   Jx_row=[];
      end
	  Jx=cat(3,Jx,Jx_col);
	  Jx_col=[];  
	 end
clear Jx_row Jx_slab
for h=1:size(Jx_index)
    BJx(:,:,h)=Jx(:,:,:,Jx_index(h));
end
Jx=BJx; clear BJx Jx_index 
Jx=squeeze(Jx);
end
%Jy
if exist('Jy_slab','var')
    Jy_slab=(Jy_slab(1:end-1,:,:,:,:)+Jy_slab(2:end,:,:,:,:))*0.5;
    Jy_slab=(Jy_slab(:,1:end-1,:,:,:)+Jy_slab(:,2:end,:,:,:))*0.5;
	Jy_slab=(Jy_slab(:,:,1:end-1,:,:)+Jy_slab(:,:,2:end,:,:))*0.5;
    Jy=[];Jy_row=[]; Jy_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Jy_row=cat(2,Jy_row,Jy_slab(:,:,:,:,MAP(i,j,k)));
        end
        Jy_col=cat(1,Jy_col,Jy_row);
	   Jy_row=[];
      end
	  Jy=cat(3,Jy,Jy_col);
	  Jy_col=[];  
	 end
clear Jy_row Jy_col Jy_slab
for h=1:size(Jy_index)
    BJy(:,:,:,h)=Jy(:,:,:,Jy_index(h));
end
Jy=BJy; clear BJy Jy_index
Jy=squeeze(Jy);
end
% Jz
if exist('Jz_slab','var')
    Jz_slab=(Jz_slab(1:end-1,:,:,:,:)+Jz_slab(2:end,:,:,:,:))*0.5;
    Jz_slab=(Jz_slab(:,1:end-1,:,:,:)+Jz_slab(:,2:end,:,:,:))*0.5;
	Jz_slab=(Jz_slab(:,:,1:end-1,:,:)+Jz_slab(:,:,2:end,:,:))*0.5;
    Jz=[];Jz_row=[]; Jz_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                Jz_row=cat(2,Jz_row,Jz_slab(:,:,:,:,MAP(i,j,k)));
        end
       Jz_col=cat(1,Jz_col,Jz_row);
	   Jz_row=[];
      end
	  Jz=cat(3,Jz,Jz_col);
	  Jz_col=[];  
	 end
clear Jz_row Jz_col Jz_slab
for h=1:size(Jz_index)
    BJz(:,:,:,h)=Jz(:,:,:,Jz_index(h));
end
Jz=BJz; clear BJz Jz_index
Jz=squeeze(Jz);
end

if exist('rho_slab','var')
    rho_slab=(rho_slab(1:end-1,:,:,:,:)+rho_slab(2:end,:,:,:,:))*0.5;
    rho_slab=(rho_slab(:,1:end-1,:,:,:)+rho_slab(:,2:end,:,:,:))*0.5;
	rho_slab=(rho_slab(:,:,1:end-1,:,:)+rho_slab(:,:,2:end,:,:))*0.5;
    rho=[];rho_row=[]; rho_col=[];
    for i=1:XLEN
        for j=1:YLEN
		  for k=1:ZLEN
                rho_row=cat(2,rho_row,rho_slab(:,:,:,:,MAP(i,j,k)));
        end
        rho_col=cat(1,rho_col,rho_row);
	    rho_row=[];
      end
	  rho=cat(3,rho,rho_col);
	  rho_col=[];  
	 end
clear rho_row rho_col rho_slab
for h=1:size(rho_index)
    Brho(:,:,:,h)=rho(:,:,:,rho_index(h));
end
rho=Brho; clear Brho rho_index
rho=squeeze(rho);
end

for nspec=1:Ns
    
   if exist(['k_nrg_slab' num2str(nspec-1)],'var')
      eval(['k_energy' num2str(nspec-1) '=sum(k_nrg_slab' num2str(nspec-1) ',2);']);
       
      for h=1:size(k_nrg_index)
        eval(['kk_energy(h)=k_energy' num2str(nspec-1) '(k_nrg_index(h));']);
      end
        eval(['k_energy' num2str(nspec-1) '=kk_energy;']);
        clear (['k_nrg_slab' num2str(nspec-1) ]); clear kk_energy
     if (nspec==Ns)
          k_energy_total=zeros(size(k_energy0));
          for nspec=1:Ns
             k_energy_total=eval(['k_energy_total+k_energy' num2str(nspec-1)]);
          end
     end
  
   
   end
            
	%Jxs
    if exist(['Jxs_slab' num2str(nspec-1)],'var')
        eval(['Jxs_slab=Jxs_slab' num2str(nspec-1) ';']);
        Jxs_slab=(Jxs_slab(1:end-1,:,:,:,:)+Jxs_slab(2:end,:,:,:,:))*0.5;
        Jxs_slab=(Jxs_slab(:,1:end-1,:,:,:)+Jxs_slab(:,2:end,:,:,:))*0.5;
        Jxs_slab=(Jxs_slab(:,:,1:end-1,:,:)+Jxs_slab(:,:,2:end,:,:))*0.5;
		Jxs=[];Jxs_row=[]; Jxs_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    Jxs_row=cat(2,Jxs_row,Jxs_slab(:,:,:,:,MAP(i,j,k)));
            end
            Jxs_col=cat(1,Jxs_col,Jxs_row);
	        Jxs_row=[];
      end
	  Jxs=cat(3,Jxs,Jxs_col);
	  Jxs_col=[];  
	 end
        clear Jxs_row Jxs_col Jxs_slab
        for h=1:size(Jxs_index)
          BJxs(:,:,:,h)=Jxs(:,:,:,Jxs_index(h));
        end
        Jxs=BJxs; clear BJxs 
        Jxs=squeeze(Jxs);
        eval(['Jxs' num2str(nspec-1) '=Jxs;']);
        clear (['Jxs_slab' num2str(nspec-1) ]); clear  Jxs  
    end

   
    if exist(['Jys_slab' num2str(nspec-1)],'var')
        eval(['Jys_slab=Jys_slab' num2str(nspec-1) ';']);
        Jys_slab=(Jys_slab(1:end-1,:,:,:,:)+Jys_slab(2:end,:,:,:,:))*0.5;
        Jys_slab=(Jys_slab(:,1:end-1,:,:,:)+Jys_slab(:,2:end,:,:,:))*0.5;
		Jys_slab=(Jys_slab(:,:,1:end-1,:,:)+Jys_slab(:,:,2:end,:,:))*0.5;
        Jys=[];Jys_row=[]; Jys_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    Jys_row=cat(2,Jys_row,Jys_slab(:,:,:,:,MAP(i,j,k)));
            end
             Jys_col=cat(1,Jys_col,Jys_row);
	        Jys_row=[];
      end
	  Jys=cat(3,Jys,Jys_col);
	  Jys_col=[];  
	 end
        clear Jys_row Jys_col Jys_slab
        for h=1:size(Jys_index)
            BJys(:,:,:,h)=Jys(:,:,:,Jys_index(h));
        end
        Jys=BJys; clear BJys
        Jys=squeeze(Jys);
        eval(['Jys' num2str(nspec-1) '=Jys;']);
        clear (['Jys_slab' num2str(nspec-1) ]); clear  Jys
    end
     %Jzs
     if exist(['Jzs_slab' num2str(nspec-1)],'var')
         eval(['Jzs_slab=Jzs_slab' num2str(nspec-1) ';']);
         Jzs_slab=(Jzs_slab(1:end-1,:,:,:,:)+Jzs_slab(2:end,:,:,:,:))*0.5;
         Jzs_slab=(Jzs_slab(:,1:end-1,:,:,:)+Jzs_slab(:,2:end,:,:,:))*0.5;
		 Jzs_slab=(Jzs_slab(:,:,1:end-1,:,:)+Jzs_slab(:,:,1:end-1,:,:))*0.5;
         Jzs=[];Jzs_row=[]; Jzs_col=[];
         for i=1:XLEN
             for j=1:YLEN
			   for k=1:ZLEN
                    Jzs_row=cat(2,Jzs_row,Jzs_slab(:,:,:,:,MAP(i,j,k)));
            end
			Jzs_col=cat(1,Jzs_col,Jzs_row);
	        Jzs_row=[]; 
      end
	  Jzs=cat(3,Jzs,Jzs_col);
	  Jzs_col=[];  
	 end
        clear Jzs_row Jzs_col Jzs_slab
        
        for h=1:size(Jzs_index)
          BJzs(:,:,:,h)=Jzs(:,:,:,Jzs_index(h));
        end
        Jzs=BJzs; clear BJzs 
        Jzs=squeeze(Jzs);
        eval(['Jzs' num2str(nspec-1) '=Jzs;']);
        clear (['Jzs_slab' num2str(nspec-1) ]); clear  Jzs  
    end

%return   
    if exist(['rhos_slab' num2str(nspec-1)],'var')
        eval(['rhos_slab=rhos_slab' num2str(nspec-1) ';']);
        rhos_slab=(rhos_slab(1:end-1,:,:,:,:)+rhos_slab(2:end,:,:,:,:))*0.5;
        rhos_slab=(rhos_slab(:,1:end-1,:,:,:)+rhos_slab(:,2:end,:,:,:))*0.5;
	rhos_slab=(rhos_slab(:,:,1:end-1,:,:)+rhos_slab(:,:,2:end,:,:))*0.5;
        rhos=[];rhos_row=[];rhos_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    rhos_row=cat(2,rhos_row,rhos_slab(:,:,:,:,MAP(i,j,k)));
            end
            rhos_col=cat(1,rhos_col,rhos_row);
	        rhos_row=[]; 
      end
	  rhos=cat(3,rhos,rhos_col);
	  rhos_col=[];  
	 end
        clear rhos_row rhos_col rhos_slab
        for h=1:size(rhos_index)
          Brhos(:,:,:,h)=rhos(:,:,:,rhos_index(h));
        end
        rhos=Brhos; clear Brhos 
        rhos=squeeze(rhos);
        eval(['rhos' num2str(nspec-1) '=rhos;']);
        clear (['rhos_slab' num2str(nspec-1) ]); clear  rhos  
    end
%return
 %pxx
   if exist(['pxx_slab' num2str(nspec-1)],'var')
       eval(['pxx_slab=pxx_slab' num2str(nspec-1) ';']);
       pxx_slab=(pxx_slab(1:end-1,:,:,:,:)+pxx_slab(2:end,:,:,:,:))*0.5;
       pxx_slab=(pxx_slab(:,1:end-1,:,:,:)+pxx_slab(:,2:end,:,:,:))*0.5;
	   pxx_slab=(pxx_slab(:,:,1:end-1,:,:)+pxx_slab(:,:,2:end,:,:))*0.5;
       pxx=[];pxx_row=[];pxx_row=[]; pxx_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pxx_row=cat(2,pxx_row,pxx_slab(:,:,:,:,MAP(i,j,k)));
            end
            pxx_col=cat(1,pxx_col,pxx_row);
	        pxx_row=[]; 
      end
	  pxx=cat(3,pxx,pxx_col);
	  pxx_col=[];  
	 end
        clear pxx_row pxx_slab pxx_col
        for h=1:size(p_index)
          Bpxx(:,:,:,h)=pxx(:,:,:,p_index(h));
        end
        pxx=Bpxx; clear Bpxx 
        pxx=squeeze(pxx); 
        eval(['pxx' num2str(nspec-1) '=pxx;']);
        clear (['pxx_slab' num2str(nspec-1) ]); clear  pxx  
   end
   % pxy
   if exist(['pxy_slab' num2str(nspec-1)],'var')
       eval(['pxy_slab=pxy_slab' num2str(nspec-1) ';']);
       pxy_slab=(pxy_slab(1:end-1,:,:,:,:)+pxy_slab(2:end,:,:,:,:))*0.5;
       pxy_slab=(pxy_slab(:,1:end-1,:,:,:)+pxy_slab(:,2:end,:,:,:))*0.5;
	   pxy_slab=(pxy_slab(:,:,1:end-1,:,:)+pxy_slab(:,:,2:end,:,:))*0.5;
       pxy=[];pxy_row=[]; pxy_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pxy_row=cat(2,pxy_row,pxy_slab(:,:,:,:,MAP(i,j,k)));
            end
            pxy_col=cat(1,pxy_col,pxy_row);
	        pxy_row=[]; 
      end
	  pxy=cat(3,pxy,pxy_col);
	  pxy_col=[];  
	 end
        clear pxy_row pxy_col pxy_slab
       
        for h=1:size(p_index)
          Bpxy(:,:,:,h)=pxy(:,:,:,p_index(h));
        end
        pxy=Bpxy; clear Bpxy 
        pxy=squeeze(pxy);
        eval(['pxy' num2str(nspec-1) '=pxy;']);
        clear (['pxy_slab' num2str(nspec-1) ]); clear  pxy  
   end
   % pxz
   if exist(['pxz_slab' num2str(nspec-1)],'var')
        eval(['pxz_slab=pxz_slab' num2str(nspec-1) ';']);
        pxz_slab=(pxz_slab(1:end-1,:,:,:,:)+pxz_slab(2:end,:,:,:,:))*0.5;
        pxz_slab=(pxz_slab(:,1:end-1,:,:,:)+pxz_slab(:,2:end,:,:,:))*0.5;
		pxz_slab=(pxz_slab(:,:,1:end-1,:,:)+pxz_slab(:,:,2:end,:,:))*0.5;
        pxz=[];pxz_row=[]; pxz_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pxz_row=cat(2,pxz_row,pxz_slab(:,:,:,:,MAP(i,j,k)));
            end
			pxz_col=cat(1,pxz_col,pxz_row);
	        pxz_row=[]; 
      end
	  pxz=cat(3,pxz,pxz_col);
	  pxz_col=[];  
	 end
        clear pxz_row pxz_row pxz_slab
        for h=1:size(p_index)
          Bpxz(:,:,:,h)=pxz(:,:,:,p_index(h));
        end
        pxz=Bpxz; clear Bpxz 
        pxz=squeeze(pxz);
        eval(['pxz' num2str(nspec-1) '=pxz;']);
        clear (['pxz_slab' num2str(nspec-1) ]); clear  pxz  
   end
   
  if exist(['pyz_slab' num2str(nspec-1)],'var')
        eval(['pyz_slab=pyz_slab' num2str(nspec-1) ';']);
        pyz_slab=(pyz_slab(1:end-1,:,:,:,:)+pyz_slab(2:end,:,:,:,:))*0.5;
        pyz_slab=(pyz_slab(:,1:end-1,:,:,:)+pyz_slab(:,2:end,:,:,:))*0.5;
        pyz=[];pyz_row=[];pyz_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pyz_row=cat(2,pyz_row,pyz_slab(:,:,:,:,MAP(i,j,k)));
            end
           pyz_col=cat(1,pyz_col,pyz_row);
	        pyz_row=[]; 
      end
	  pyz=cat(3,pyz,pyz_col);
	  pyz_col=[];  
	 end
        clear pyz_row pyz_col pyz_slab
        for h=1:size(p_index)
          Bpyz(:,:,:,h)=pyz(:,:,:,p_index(h));
        end
        pyz=Bpyz; clear Bpyz         
        pyz=squeeze(pyz);
        eval(['pyz' num2str(nspec-1) '=pyz;']);
        clear (['pyz_slab' num2str(nspec-1) ]); clear  pyz  
  end

  if exist(['pyy_slab' num2str(nspec-1)],'var')
        eval(['pyy_slab=pyy_slab' num2str(nspec-1) ';']);
        pyy_slab=(pyy_slab(1:end-1,:,:,:,:)+pyy_slab(2:end,:,:,:,:))*0.5;
        pyy_slab=(pyy_slab(:,1:end-1,:,:,:)+pyy_slab(:,2:end,:,:,:))*0.5;
		pyy_slab=(pyy_slab(:,:,1:end-1,:,:)+pyy_slab(:,:,2:end,:,:))*0.5;
        pyy=[];pyy_row=[];pyy_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pyy_row=cat(2,pyy_row,pyy_slab(:,:,:,:,MAP(i,j,k)));
            end
            pyy_col=cat(1,pyy_col,pyy_row);
	        pyy_row=[]; 
      end
	  pyy=cat(3,pyy,pyy_col);
	  pyy_col=[];  
	 end
        clear pyy_row pyy_col pyy_slab
        for h=1:size(p_index)
          Bpyy(:,:,:,h)=pyy(:,:,:,p_index(h));
        end
        pyy=Bpyy; clear Bpyy
        pyy=squeeze(pyy);
        eval(['pyy' num2str(nspec-1) '=pyy;']);
        clear (['pyy_slab' num2str(nspec-1) ]); clear  pyy  
  end
    
  if exist(['pzz_slab' num2str(nspec-1)],'var')
        eval(['pzz_slab=pzz_slab' num2str(nspec-1) ';']);
        pzz_slab=(pzz_slab(1:end-1,:,:,:,:)+pzz_slab(2:end,:,:,:,:))*0.5;
        pzz_slab=(pzz_slab(:,1:end-1,:,:,:)+pzz_slab(:,2:end,:,:,:))*0.5;
		pzz_slab=(pzz_slab(:,:,1:end-1,:,:)+pzz_slab(:,:,2:end,:,:))*0.5;
        pzz=[];pzz_row=[]; pzz_col=[];
        for i=1:XLEN
            for j=1:YLEN
			  for k=1:ZLEN
                    pzz_row=cat(2,pzz_row,pzz_slab(:,:,:,:,MAP(i,j,k)));
            end
             pzz_col=cat(1,pzz_col,pzz_row);
	        pzz_row=[]; 
      end
	  pzz=cat(3,pzz,pzz_col);
	  pz_col=[];  
	 end
        clear pzz_row ppz_col pzz_slab
        for h=1:size(p_index)
          Bpzz(:,:,:,h)=pzz(:,:,:,p_index(h));
        end
        pzz=Bpzz; clear Bpzz 
        pzz=squeeze(pzz);
        eval(['pzz' num2str(nspec-1) '=pzz;']);
        clear (['pzz_slab' num2str(nspec-1) ]); clear  pzz  
  end
  
  if exist(['x_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['xx.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['xx.time' num2str(time) '=cat(1,xx.time' num2str(time) ',x_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['xxsize(time)=size(xx.time' num2str(time) ',1);'])  
        end
	m=min(xxsize);
        eval(['x' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['x' num2str(nspec-1) '(:,time) =xx.time' num2str(time) '(1:m);'])
        end
   clear (['x_slab' num2str(nspec-1)]); clear xx xxsize
    for h=1:size(x_index)
          eval(['xx(:,h)=x' num2str(nspec-1) '(:,x_index(h));'])
    end
    eval(['x' num2str(nspec-1) '= xx;'])
    clear xx m
  end
  
  if exist(['y_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['yy.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['yy.time' num2str(time) '=cat(1,yy.time' num2str(time) ',y_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['yysize(time)=size(yy.time' num2str(time) ',1);'])  
        end
	m=min(yysize);
        eval(['y' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['y' num2str(nspec-1) '(:,time) =yy.time' num2str(time) '(1:m);'])
        end
    clear (['y_slab' num2str(nspec-1)]); clear yy yysize  
    for h=1:size(y_index)
          eval(['yy(:,h)=y' num2str(nspec-1) '(:,y_index(h));'])
    end
    eval(['y' num2str(nspec-1) '= yy;'])
    clear yy m
  end
  if exist(['z_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['zz.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['zz.time' num2str(time) '=cat(1,zz.time' num2str(time) ',z_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['zzsize(time)=size(zz.time' num2str(time) ',1);'])  
        end
	m=min(zzsize);
        eval(['z' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['z' num2str(nspec-1) '(:,time) =zz.time' num2str(time) '(1:m);'])
        end
    clear (['z_slab' num2str(nspec-1)]); clear zz zzsize
    for h=1:size(z_index)
          eval(['zz(:,h)=z' num2str(nspec-1) '(:,z_index(h));'])
    end
    eval(['z' num2str(nspec-1) '= zz;'])
    clear zz m
  end
  if exist(['u_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['uu.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['uu.time' num2str(time) '=cat(1,uu.time' num2str(time) ',u_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['uusize(time)=size(uu.time' num2str(time) ',1);'])  
        end
	m=min(uusize);
        eval(['u' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['u' num2str(nspec-1) '(:,time) =uu.time' num2str(time) '(1:m);'])
        end
    clear (['u_slab' num2str(nspec-1)]); clear uu uusize
    for h=1:size(u_index)
          eval(['uu(:,h)=u' num2str(nspec-1) '(:,u_index(h));'])
    end
    eval(['u' num2str(nspec-1) '= uu;'])
    clear uu m
  end
  if exist(['v_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['vv.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['vv.time' num2str(time) '=cat(1,vv.time' num2str(time) ',v_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['vvsize(time)=size(vv.time' num2str(time) ',1);'])  
        end
	m=min(vvsize);
        eval(['v' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['v' num2str(nspec-1) '(:,time) =vv.time' num2str(time) '(1:m);'])
        end
    clear (['v_slab' num2str(nspec-1)]); clear vv vvsize
    for h=1:size(v_index)
          eval(['vv(:,h)=v' num2str(nspec-1) '(:,v_index(h));'])
    end
    eval(['v' num2str(nspec-1) '= vv;'])
    clear vv m
  end
  if exist(['w_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['ww.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['ww.time' num2str(time) '=cat(1,ww.time' num2str(time) ',w_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['wwsize(time)=size(ww.time' num2str(time) ',1);'])  
        end
	m=min(wwsize);
        eval(['w' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['w' num2str(nspec-1) '(:,time) =ww.time' num2str(time) '(1:m);'])
        end
    clear (['w_slab' num2str(nspec-1)]); clear ww wwsize
    for h=1:size(w_index)
          eval(['ww(:,h)=w' num2str(nspec-1) '(:,w_index(h));'])
    end
    eval(['w' num2str(nspec-1) '= ww;'])
    clear ww m
  end
  if exist(['q_slab' num2str(nspec-1)],'var')
     
        for time=1:ntime
            eval(['qq.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['qq.time' num2str(time) '=cat(1,qq.time' num2str(time) ',q_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['qqsize(time)=size(qq.time' num2str(time) ',1);'])  
        end
	m=min(qqsize);
        eval(['q' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['q' num2str(nspec-1) '(:,time) =qq.time' num2str(time) '(1:m);'])
        end
    clear (['q_slab' num2str(nspec-1)]); clear qq qqsize
    for h=1:size(q_index)
          eval(['qq(:,h)=q' num2str(nspec-1) '(:,q_index(h));'])
    end
    eval(['q' num2str(nspec-1) '= qq;'])
    clear qq m
  end
  if exist(['ID_slab' num2str(nspec-1)],'var')
 
        for time=1:ntime
            eval(['IDID.time' num2str(time) '=[];']);
         for PROC=1:Nprocs
             eval(['IDID.time' num2str(time) '=cat(1,IDID.time' num2str(time) ',ID_slab' num2str(nspec-1) '.time' num2str(time) '.PROC' num2str(PROC) ');']);
         end
          eval(['IDIDsize(time)=size(IDID.time' num2str(time) ',1);'])  
        end
	m=min(IDIDsize);
        eval(['ID' num2str(nspec-1) '=NaN(m,ntime);'])
        for time=1:ntime
            eval(['ID' num2str(nspec-1) '(:,time) =IDID.time' num2str(time) '(1:m);'])
        end
    clear (['ID_slab' num2str(nspec-1)]); clear IDID IDIDsize
    for h=1:size(ID_index)
          eval(['IDID(:,h)=ID' num2str(nspec-1) '(:,ID_index(h));'])
    end
    eval(['ID' num2str(nspec-1) '= IDID;'])
    clear IDID m
  end
  
end

clear i ii iii j k h *_index
