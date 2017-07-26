info=hdf5info('/shared/gianni/tred26/proc111.hdf');
nt=max(size(info.GroupHierarchy.Groups(1).Groups(1).Datasets(:)))
for it=1:nt
% fieldz
info.GroupHierarchy.Groups(1).Groups(1).Datasets(it).Name
bx=hdf5read(info.GroupHierarchy.Groups(1).Groups(1).Datasets(it));
by=hdf5read(info.GroupHierarchy.Groups(1).Groups(2).Datasets(it));
bz=hdf5read(info.GroupHierarchy.Groups(1).Groups(3).Datasets(it));
% moments
Jxs0=hdf5read(info.GroupHierarchy.Groups(2).Groups(2).Groups(1).Datasets(it));
Jys0=hdf5read(info.GroupHierarchy.Groups(2).Groups(2).Groups(2).Datasets(it));
Jzs0=hdf5read(info.GroupHierarchy.Groups(2).Groups(2).Groups(3).Datasets(it));
rhos0=hdf5read(info.GroupHierarchy.Groups(2).Groups(2).Groups(4).Datasets(it));
Jxs2=hdf5read(info.GroupHierarchy.Groups(2).Groups(4).Groups(1).Datasets(it));
Jys2=hdf5read(info.GroupHierarchy.Groups(2).Groups(4).Groups(2).Datasets(it));
Jzs2=hdf5read(info.GroupHierarchy.Groups(2).Groups(4).Groups(3).Datasets(it));
rhos2=hdf5read(info.GroupHierarchy.Groups(2).Groups(4).Groups(4).Datasets(it));


ux=(Jxs2+Jxs2)./(rhos0+rhos2);
uy=(Jys2+Jys2)./(rhos0+rhos2);
uz=(Jys2+Jys2)./(rhos0+rhos2);
bb=sqrt(bx.^2+by.^2+bz.^2);
upar=(ux.*bx+uy.*by+uz.*bz)./bb;
uperx=ux-upar.*bx./bb;
upery=uy-upar.*by./bb;
uperz=uz-upar.*bz./bb;
uperp=sqrt(uperx.^2+upery.^2+uperz.^2);

%plot(squeeze(uperp(1,:,:)))
pcolor(squeeze(uperp(:,1,:)))
colorbar
pause
end
