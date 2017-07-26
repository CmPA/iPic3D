% need to use parsek3D to read it

[nx ny nz nt]=size(Bx);

indexf = 1;
for it=1:nt
% electron density
figure(1)
title('electron density');
rhos0temp = rhos0(:,:,:,it);
h1 = vol3d('cdata',rhos0temp,'texture','3D');
view(-75.0,30); 
%alphamap(.06 .*alphamap); % control the non-transparent region
vol3d(h1);
F10(indexf)=getframe(gcf);

%different view
figure(10)
h10 = vol3d('cdata',permute(rhos0temp,[3 2 1]),'texture','3D');
view(3); 
%alphamap(.06 .*alphamap); % control the non-transparent region
vol3d(h10);
F1(indexf)=getframe(gcf);


% Jz for electron should mimic the Az
Jzs0temp = Jzs0(:,:,:,it);
figure(20)
h20 = vol3d('cdata',permute(Jzs0temp,[3 2 1]),'texture','3D');
view(3); 
alphamap(.06 .*alphamap);
vol3d(h20)
F2(indexf)= getframe(gcf);

% Magnetic Intensity

moduloB = Bx(:,:,:,it).^2+By(:,:,:,it).^2 + Bz(:,:,:,it).^2;
figure(30)
h30 = vol3d('cdata',permute(moduloB,[3 2 1]),'texture','3D');
view(3);
%alphamap(.06 .*alphamap);
vol3d(h30)
F3(indexf)= getframe(gcf);

% electric field out of plane
Eztemp = Ez(:,:,:,it);
figure(40)
h40 = vol3d('cdata',permute(Eztemp,[3 2 1]),'texture','3D');
view(3); 
alphamap(.06 .*alphamap);
vol3d(h40)
F4(indexf)= getframe(gcf);


% contour of the midplane
figure(5)
subplot(3,1,1)
pcolor(rhos0temp(:,:,nz/2));
shading interp
subplot(3,1,2)
pcolor(Jzs0temp(:,:,nz/2));
shading interp
subplot(3,1,3)
pcolor(Eztemp(:,:,nz/2));
shading interp
F5(indexf) =getframe(gcf);

close all
indexf = indexf + 1

end %end of the iteration

% make movies
movie2avi(F1,'rho_el_3D.avi','fps',[2],'quality',[100])
movie2avi(F10,'rho_el_3D_diffview.avi','fps',[2],'quality',[100])
movie2avi(F2,'Jz_el_3D.avi','fps',[2],'quality',[100])
movie2avi(F3,'B_intensity_3D.avi','fps',[2],'quality',[100])
movie2avi(F4,'Ez_3D.avi',fps',[2],'quality',[100])
movie2avi(F5,'midplane_contour.avi',fps',[2],'quality',[100])


