function [J]=immagine_spacetime(x,y,J1,name,asse,lmt, Nsm,Nfig)
 
J=smooth_1d(J1,Nsm);

h=figure(Nfig)
%set(h,'Position',[167 451 976 351])
Ncut=max(Nsm*3,1);
if(lmt(1) == lmt(2))
lmt(1)=min(J1(:));
lmt(2)=max(J1(:));
end
imagesc(x,y,J(1:end,Ncut:end-Ncut)',lmt)
%colormap hot
colorbar
set(gca,'fontsize',[14])
xlabel(asse,'fontsize',[14])
ylabel('t(s)','fontsize',[14])
title(name,'fontsize',[14])
%axis image
axis xy
%set(gca,'xdir','reverse','TickDir','out')
set(gca,'TickDir','out')
%print('-depsc','-r300',[name '.eps'])
%print('-dpng','-r300',[name '.png'])
%saveas(gcf,[name '.fig'])
