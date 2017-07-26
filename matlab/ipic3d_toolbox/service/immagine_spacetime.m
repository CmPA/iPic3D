function [J]=immagine_spacetime(x,y,J1,name,asseX,asseY,lmt, Nsm,Nfig)
 
J=smooth_1d(J1,Nsm);

h=figure(Nfig)
%set(h,'Position',[167 451 976 351])
Ncut=max(Nsm*3,1);
if(lmt(1) == lmt(2))
lmt(1)=min(J1(:));
lmt(2)=max(J1(:));
end
imagesc(x,y,J(2:end-1,Ncut:end-Ncut)',lmt)
%colormap hot
colorbar
set(gca,'fontsize',[14])
xlabel(asseX,'fontsize',[14])
ylabel(asseY,'fontsize',[14])
title(name,'fontsize',[14])
%axis image
axis xy
%set(gca,'xdir','reverse','TickDir','out')
set(gca,'TickDir','out')
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
saveas(gcf,[name '.fig'])


h=figure(100+Nfig)
pcolor(y(Ncut:end-Ncut),x(2:end-1),J(2:end-1,Ncut:end-Ncut))
caxis(lmt)
shading interp
colorbar
set(gca,'fontsize',[14])
ylabel(asseX,'fontsize',[14])
xlabel(asseY,'fontsize',[14])
title(name,'fontsize',[14])
print('-dpng','-r300',[name '_pcolor.png'])
saveas(gcf,[name '_pcolor.fig'])

