%results_dir='/nobackup/gianni/run3D2/' ; % directory for results
%results_dir='/nobackup/gianni/3Dbackbg/' ; % directory for results
results_dir='/shared/gianni/tred54/' ; % directory for results
%results_dir='/nobackup/gianni/run3D1/' ; % directory for results

processor_name=[results_dir 'proc'];
info=hdf5info([processor_name,num2str(0),'.hdf']);
nslicesintime=size(info.GroupHierarchy.Groups(1).Groups(1).Datasets,2)

!rm film*/*
ig=0;
for reader_counter=1:nslicesintime
	time_counter_list=reader_counter
variable_list='B';
	parsek3D_select
vthe=uth(1)
vthi=uth(2)
va=wci;
wci=Bx0
time=double(Bx_time)*wci*Dt
        info=hdf5info([processor_name,num2str(0),'.hdf']);
        ig=str2double(regexprep(info.GroupHierarchy.Groups(1).Groups(1).Datasets(reader_counter).Name,'/fields/Bx/cycle_',''))
%	sec3Dxy
	sec3DxzB
%	sec3Dyz	
end 
