[r, th, fi]=ndgrid(0:vmax/30:vmax*.9, 0:pi/10:pi, 0:pi/10:2*pi);
v2=r.*sin(th).*sin(fi);
v3=r.*sin(th).*cos(fi);
v1=r.*cos(th);

[vx,vy,vz]=ndgrid(1:ndiv,1:ndiv,1:ndiv);

vx=(vx-(ndiv-1)/2)/ndiv*vmax*2;
vy=(vy-(ndiv-1)/2)/ndiv*vmax*2;
vz=(vz-(ndiv-1)/2)/ndiv*vmax*2;
vv = vx.^2 +vy.^2 +vz.^2;

vdf_sp_r=interpn(vx,vy,vz,vv.*vdf_sp, v1,v2,v3);

for nsl=1:0
subplot(2,1,1)
imagesc(squeeze(vdf_sp_r(nsl,:,:)));title(num2str(r(nsl,1,1)))
subplot(2,1,2)
plot(squeeze(mean(vdf_sp_r(nsl,:,:),2)))
pause
end
close all
imagesc([0 max(fi(:))],[0 max(v1(:))],(squeeze(sum(vdf_sp_r,2))))
title([pos_label '  Agyrotropy = ' num2str(agyrotropy,3)])
xlabel('gyrophase','fontsize',[14])
ylabel('v/c','fontsize',[14])
set(gca,'fontsize',[14])
%caxis([0 5]*1e-3)
cb=colorbar('south')
cb.Color='g'

  %load rainbow_cm
  colormap(cmap)
  
filename=['EF_thdependence_' 'species_' num2str(is) '_' num2str(ipy*nsuby+suby) '.png'];
print(filename,'-dpng');
