function [J]=immagine_dir(x,y,J1,name,lmt, Nsm, Ncycle, Ygsm, nfig, nplot, labelx,labely, labelc)
global color_choice symmetric_color titolo square

load cm_new
load cm_multi4

J=smooth(J1,Nsm);

h=figure(nfig)
if (nplot>1) 
    subplot(1,3,nplot) 
end
if(nplot==1) 
    set(h,'Position',[167 51 3*500*max(x)/max(y) 500]);
end

Ncut=max(Nsm*3,1);

if(lmt(1) == lmt(2))
lmt(1)=min(J(:));
lmt(2)=max(J(:));
if(symmetric_color==1)
lmt(1)=-max(abs(lmt));
lmt(2)=-lmt(1);
elseif (symmetric_color==-1)
lmt(1)=-min(abs(lmt));
lmt(2)=-lmt(1);
elseif (symmetric_color==2)
lmt(2)=(sqrt(mean(J(:).^2))+max(abs(J(:))))/2;
lmt(1)=-lmt(2);
elseif (symmetric_color==-2)
lmt(2)=(mean(abs(J(:))) + max(abs(J(:))))/2;
lmt(1)=0.0;
end
end

lmt
imagesc(x,y,J(Ncut:end-Ncut,Ncut:end-Ncut)',lmt);

        if(color_choice==0)
        colormap jet
        elseif (color_choice==1)
        colormap(cm_kbwrk)
        elseif (color_choice==2)
        colormap(cool)
        elseif (color_choice==3)
        colormap(cm_cool_hot_2)
        elseif (color_choice==4)
        colormap(flipud(hot))
        end
        
        
        c = colorbar
        ylabel(c,labelc)
        set(gca,'fontsize',[14])
set(gca,'fontsize',[14])
xlabel(labelx,'fontsize',[14])
ylabel(labely,'fontsize',[14])

%if isempty(titolo)
%        title(['\omega_{ci}t= ' Ncycle], 'fontsize',[14])       
%else
        title(titolo, 'fontsize',[14])
%end

        
if (square) 
    axis xy
    axis square 
else
axis xy
end
%set(gca,'xdir','reverse','TickDir','out')
%set(gca,'TickDir','out')
%print('-depsc','-r300',[name '.eps'])
%print('-dpng','-r300',[name '.png'])
%saveas(gcf,[name '.fig'])
