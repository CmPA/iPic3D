function [J]=immagine(x,y,J1,name,lmt, Nsm, Ncycle, Ygsm, nfig, labelx,labely)
global color_choice symmetric_color Dt
if(nargin>6)
    Ncycle
cycle=str2num(Ncycle)
% for HRmaha3D1:
 time=60*(cycle/75000.0*Dt/.125) %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D

%ADD initial time of the RUN
time=time+(03*60+48)*60
% Prepare string
ntime = datestr(time/86400,'HH:MM:SS UT')
%ntime=num2str(time,'%5.2f')

else
ntime=[]
end 
J=smooth(J1,Nsm);
if(nargin>8)
h=figure(nfig)
else
h=figure(1)
end
ar=(max(x(:))-min(x(:)))/(max(x(:))-min(x(:)));
set(h,'Position',[167 451 167+351*ar 351])
if(lmt(1) == lmt(2))
lmt(1)=min(J(:));
lmt(2)=max(J(:));
if(symmetric_color==1)
lmt(1)=-max(abs(lmt))
lmt(2)=-lmt(1)
elseif (symmetric_color==-1)
lmt(1)=-min(abs(lmt))
lmt(2)=-lmt(1)
end
end
Ncut=max(Nsm*3,1)
imagesc(x,y,J(Ncut:end-Ncut,Ncut:end-Ncut)',lmt)
if(color_choice==0)
        colormap jet
elseif (color_choice==1)
        load cm_new
        colormap(cm_kbwrk)
elseif (color_choice==2)
        colormap hot
end
colorbar
set(gca,'fontsize',[14])
if(nargin>9)
xlabel(labelx,'fontsize',[14])
ylabel(labely,'fontsize',[14])
else
xlabel('x/R_E','fontsize',[14])
ylabel('z/R_E','fontsize',[14])
end

if(nargin>7 & max(size(Ygsm))>0)
        title(['time (s) = ' ntime '  Ygsm=', num2str(Ygsm)], 'fontsize',[14])
elseif(nargin>6)
        title(['time (s) = ' ntime], 'fontsize',[14]),
end

axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%saveas(gcf,[name '.fig'])
