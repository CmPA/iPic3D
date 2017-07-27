function [J]=common_image_due(x,y,J1,Az,Psi,nlabel,name,clmt, Nsm, nfig)

global color_choice symmetric_color labelx labely Ncycle


J=smooth(J1,Nsm);

if(nargin>8)
h=figure(nfig)
else
h=figure(1)
end

set(h, 'Position', [167 26 515 776])

if(clmt(1) == clmt(2))
clmt(1)=min(J(:));
clmt(2)=max(J(:));
if(symmetric_color==1)
clmt(1)=-max(abs(clmt))
clmt(2)=-clmt(1)
elseif (symmetric_color==-1)
clmt(1)=-min(abs(clmt))
clmt(2)=-clmt(1)
end
end


colorcontour='k--'
if(color_choice==0)
        colormap jet
elseif (color_choice==1)
        load cm_new
        colormap(cm_kbwrk)
elseif (color_choice==2)
        colormap hot
elseif (color_choice==-1)
        colormap parula
        colorcontour='w'
end

Ncut=max(Nsm*3,1)

%imagesc([min(x(:)) max(x(:))],[min(y(:)) max(y(:))],J(Ncut:end-Ncut,Ncut:end-Ncut)',clmt)
imagesc(x(1,:),y(:,1),J(Ncut:end-Ncut,Ncut:end-Ncut)',clmt)
hold on
contour(x,y,Az',30,'k','linewidth',2)
contour(x,y,Psi',21,'k:','linewidth',2)
caxis(clmt)


colorbar
set(gca,'fontsize',[14])

xlabel(labelx,'fontsize',[14])
ylabel(labely,'fontsize',[14])
title([ nlabel '   ' Ncycle],'fontsize',[14])


axis image
axis xy
%set(gca,'xdir','reverse','TickDir','out')
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name Ncycle '.png'])
%saveas(gcf,[name '.fig'])
