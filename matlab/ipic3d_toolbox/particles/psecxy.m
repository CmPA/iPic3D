results_dir='/shared/gianni/run120.2/';
filename=[results_dir 'settings.hdf'];
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
B0x=double(hdf5read(filename,'/collective/B0x'))
Dt=double(hdf5read(filename,'/collective/Dt'))
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
Nprocs=double(hdf5read(filename,'/topology/Nprocs'))
ipx=XLEN/2+5;
ipy=YLEN/2;
ip=(ipx-1)*YLEN+ipy
xmin=double(0+(ipx-1)*Lx/XLEN)
xmax=double(xmin+Lx/XLEN)
ymin=double(0+(ipy-1)*Ly/YLEN)
ymax=double(ymin+2*double(Ly/YLEN))

info=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

[nt dummy]=size(info.GroupHierarchy.Groups(2).Groups(4).Datasets(:));
%info2=hdf5info([results_dir 'proc1040.hdf']);
info2=hdf5info([results_dir 'proc' num2str(ip+1) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(2).Groups(4).Datasets(:));
for it=1:nt
ex1=squeeze(hdf5read(info.GroupHierarchy.Groups(2).Groups(4).Datasets(it)));
ex2=squeeze(hdf5read(info2.GroupHierarchy.Groups(2).Groups(4).Datasets(it)));
ig=str2num(regexprep(info.GroupHierarchy.Groups(2).Groups(4).Datasets(it).Name,'/fields/Ex/cycle_',''))
[nx ny]=size(ex1);
dx=double(xmax-xmin)/double(nx-1);
dy=double(ymax-ymin)/double(2*ny-2);
[x,y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
%ex=[(ex1(1:end-1,1:end-1)+ex1(2:end,2:end)+ex1(1:end-1,2:end)+ex1(1:end-1,2:end))/4; (ex2(1:end-1,1:end-1)+ex2(2:end,2:end)+ex2(1:end-1,2:end)+ex2(1:end-1,2:end))/4];
%pcolor(ex)
pcolor(x,y,[ex1(1:end-1,:);ex2])
colorbar
shading interp
xlabel('x')
ylabel('y')
title(['E_x(\omega_{ci}t=' num2str(ig*Dt*B0x) ')'])
set(gcf,'Renderer','zbuffer')
print('-dpng',['pfilmEx/' num2str(ig,'%5.5i')])
saveas(gcf,['pfilmEx/' num2str(ig,'%5.5i') '.fig'])
pause(.1)
end
