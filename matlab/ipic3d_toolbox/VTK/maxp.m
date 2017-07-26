read=1
if(read)
[V,Vx,Vy,Vz,dx,dy,dz]=read_vtk_3d('B_xyz_cycle40000.vtk');
numvar=10;
filename='Pi_xyz_cycle40000.vtk';
[Pxx,Pxy,Pxz,Pyy,Pyz,Pzz,Ppar,Pper1,Pper2,Peps]=read_vtk_multiscalar_3d(filename,numvar);
end
[nx ny nz]= size(Pxx)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
ymax(i,k)=round(sum(w.*y)./sum(w));
%ymax(i,k)=round(sum(y./V(i,:,k).^2)./sum(1./V(i,:,k).^2));
%[dum j] = min(Vx(i,:,k).^2);
%ymax(i,k)=j;
Vxmax(i,k)=Vx(i,round(ymax(i,k)),k);
%Vxmax(i,k)=Vx(i,end/2,k);
Vymax(i,k)=Vy(i,round(ymax(i,k)),k);
Vzmax(i,k)=Vz(i,round(ymax(i,k)),k);
%Pparmax(i,k)=Pxx(i,end/2,k);
end
end
subplot(3,1,1)
%surf(ymax)
pcolor(ymax')
shading interp
colorbar
subplot(3,1,2)
pcolor(Vymax')
caxis([-1 1]*1e-3)
shading interp
colorbar
subplot(3,1,3)
pcolor(Vzmax')
colorbar
shading interp
set(gcf,'Renderer','zbuffer');
print -dpng caz.png

