function read_parsek_settings(filename)

global Lx Ly Lz Nx Ny Nz Dx Dy Dz Dt Th Ncycles Ns c Smooth
global PfaceXright PfaceXleft PfaceYright PfaceYleft PfaceZright PfaceZleft
global PHIfaceXright PHIfaceXleft PHIfaceYright PHIfaceYleft PHIfaceZright PHIfaceZleft 
global EMfaceXright EMfaceXleft EMfaceYright EMfaceYleft EMfaceZright EMfaceZleft

global Np Npcelx Npcely Npcelz NpMax qom 
global uth vth wth u_drift v_drift w_drift
global XLEN YLEN ZLEN Nprocs periodicX periodicY periodicZ

Lx=hdf5read(filename,'/collective/Lx');
Ly=hdf5read(filename,'/collective/Ly');
Lz=hdf5read(filename,'/collective/Lz');
Nx=hdf5read(filename,'/collective/Nxc');
Ny=hdf5read(filename,'/collective/Nyc');
Nz=hdf5read(filename,'/collective/Nzc');
Dx=hdf5read(filename,'/collective/Dx');
Dy=hdf5read(filename,'/collective/Dy');
Dz=hdf5read(filename,'/collective/Dz');
Dt=hdf5read(filename,'/collective/Dt');
Th=hdf5read(filename,'/collective/Th');
Ncycles=hdf5read(filename,'/collective/Ncycles');
Ns=hdf5read(filename,'/collective/Ns');
c=hdf5read(filename,'/collective/c');
Smooth=hdf5read(filename,'/collective/Smooth');

PfaceXright=int8(hdf5read(filename,'/collective/bc/PfaceXright'));
PfaceXleft=int8(hdf5read(filename,'/collective/bc/PfaceXleft'));
PfaceYright=int8(hdf5read(filename,'/collective/bc/PfaceYright'));
PfaceYleft=int8(hdf5read(filename,'/collective/bc/PfaceYleft'));
PfaceZright=int8(hdf5read(filename,'/collective/bc/PfaceZright'));
PfaceZleft=int8(hdf5read(filename,'/collective/bc/PfaceZleft'));

PHIfaceXright=int8(hdf5read(filename,'/collective/bc/PHIfaceXright'));
PHIfaceXleft=int8(hdf5read(filename,'/collective/bc/PHIfaceXleft'));
PHIfaceYright=int8(hdf5read(filename,'/collective/bc/PHIfaceYright'));
PHIfaceYleft=int8(hdf5read(filename,'/collective/bc/PHIfaceYleft'));
PHIfaceZright=int8(hdf5read(filename,'/collective/bc/PHIfaceZright'));
PHIfaceZleft=int8(hdf5read(filename,'/collective/bc/PHIfaceZleft'));

EMfaceXright=int8(hdf5read(filename,'/collective/bc/EMfaceXright'));
EMfaceXleft=int8(hdf5read(filename,'/collective/bc/EMfaceXleft'));
EMfaceYright=int8(hdf5read(filename,'/collective/bc/EMfaceYright'));
EMfaceYleft=int8(hdf5read(filename,'/collective/bc/EMfaceYleft'));
EMfaceZright=int8(hdf5read(filename,'/collective/bc/EMfaceZright'));
EMfaceZleft=int8(hdf5read(filename,'/collective/bc/EMfaceZleft'));

for i=0:Ns-1
   
    Np(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/Np']);
    Npcelx(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/Npcelx']);
    Npcely(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/Npcely']);
	Npcelz(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/Npcelz']); 
    NpMax(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/NpMax']);
    qom(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/qom']);
    uth(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/uth']);
    vth(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/vth']);
    wth(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/wth']);
    u_drift(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/u0']);
    v_drift(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/v0']);
    w_drift(i+1)=hdf5read(filename,['/collective/species_',num2str(i),'/w0']);


end

XLEN=hdf5read(filename,'/topology/XLEN');
YLEN=hdf5read(filename,'/topology/YLEN');
ZLEN=hdf5read(filename,'/topology/ZLEN');
Nprocs=hdf5read(filename,'/topology/Nprocs');
periodicX=hdf5read(filename,'/topology/periodicX');
periodicY=hdf5read(filename,'/topology/periodicY');
periodicZ=hdf5read(filename,'/topology/periodicZ');


return
