function []=immagine_vv(x,y,J1, Nsm )

J=smooth(J1,Nsm);
Ncut=1
imagesc(x,y,log10(J(Ncut:end-Ncut,Ncut:end-Ncut)) )
%colormap hot
%load cm_new
%colormap(cm_kbwrk)
colormap jet

end
