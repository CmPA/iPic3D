
close all


for cycle=Ncyc_ini:1000:Ncyc_max
ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

global blowup contours
blowup=0;
contours=1;



x=fliplr(Xgsmrange);
y=Zgsmrange;


close all
file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny);

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny);

Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);
else
xc=linspace(-45,-15,Nx);
yc=linspace(-9,3,Ny);
[xc,yc]=meshgrid(xc,yc);
end

immagine(x,y,By*code_B,['Bn' ncycle1],[-10 10],0,ncycle1,mean(Ygsmrange),1,'x/R_E','z/R_E')
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Bn_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,Bx*code_B,['Bl' ncycle1],[-25 25],0,ncycle1,mean(Ygsmrange),1,'x/R_E','z/R_E')
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
        %ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Bl_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,Bz*code_B,['Bm' ncycle1],[-10 10],0,ncycle1,mean(Ygsmrange),1,'x/R_E','z/R_E')
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(gsmx(xc),gsmy2z(yc),Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
%starty=linspace(0,Ly,30);
%startx=ones(size(starty))*Ly/4;
%streamline(Bx*dx/dy,By,startx,starty)

hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
end

set(gcf, 'Renderer', 'zbuffer');
name=['Bm_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])
end
