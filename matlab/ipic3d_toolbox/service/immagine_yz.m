function []=immagine(x,y,J1,name, lmt, Nsm, Ncycle, xpos, nfig,labelx,labely)
if(nargin>6)
cycle=str2num(Ncycle)
time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D
ntime=num2str(time,'%5.2f')
else
ntime=[]
end 
J=smooth(J1,Nsm);
if(nargin>7)
h=figure(nfig)
else
h=figure(1)
end
ar=(max(x(:))-min(x(:)))/(max(x(:))-min(x(:)));
set(h,'Position',[167 451 167+351*ar 351])
if(lmt(1) == lmt(2))
lmt(1)=min(J(:));
lmt(2)=max(J(:));
end
Ncut=max(Nsm*3,1)
imagesc(x,y,J(Ncut:end-Ncut,Ncut:end-Ncut)',lmt)
%colormap hot
%load cm_new
%colormap(cm_kbwrk)
colormap jet
colorbar
set(gca,'fontsize',[14])
if(nargin>6) title([name  '  time (s) = ' ntime '   x(GSM)=' num2str(gsmx(xpos))], 'fontsize',[14]), end
if(nargin>8)
xlabel(labelx,'fontsize',[14])
ylabel(labely,'fontsize',[14])
else
xlabel('x/R_E','fontsize',[14])
ylabel('z/R_E','fontsize',[14])
end
axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%saveas(gcf,[name '.fig'])
