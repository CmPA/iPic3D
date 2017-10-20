addpath(genpath('../../ipic3d_toolbox'))
global contours
contours = 1;
for i=18000:1000:18000
disp(['CYCLE = ' num2str(i)])
ncyc=num2str(i);
read=1;
if (read) 
    [V,Bx,By,Bz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGdB2dt_AVG_xy_cycle' ncyc '.vtk']);
    [V,dVx,dVy,dVz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGdeldB2dt_AVG_xy_cycle' ncyc '.vtk']);
Bx=Bx+dVx;By=By+dVy;Bz=Bz+dVz;
[nx,ny]=size(Bx);
B2cut=Bx(5:end-5,5:end-5)/4/pi+By(5:end-5,5:end-5)/4/pi+Bz(5:end-5,5:end-5)/4/pi;
clear Bx By Bz

    [V,JEex,JEey,JEez,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGJdotEe_AVG_xy_cycle' ncyc '.vtk']);
    [V,dVx,dVy,dVz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGdelJdotEe_AVG_xy_cycle' ncyc '.vtk']);
JEex=JEex+dVx; JEey=JEey+dVy; JEez=JEez+dVz;
JEecut=JEex(5:end-5,5:end-5)+JEey(5:end-5,5:end-5)+JEez(5:end-5,5:end-5);
clear JEex JEey JEez

    [V,JEix,JEiy,JEiz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGJdotEi_AVG_xy_cycle' ncyc '.vtk']);
    [V,dVx,dVy,dVz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGdelJdotEi_AVG_xy_cycle' ncyc '.vtk']);
JEix=JEix+dVx; JEiy=JEiy+dVy; JEiz=JEiz+dVz;
JEicut=JEix(5:end-5,5:end-5)+JEiy(5:end-5,5:end-5)+JEiz(5:end-5,5:end-5);
clear JEix JEiy JEiz
    [divS,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data1/AVGdivS_AVG_xy_cycle' ncyc '.vtk']);
divScut=divS(5:end-5,5:end-5)/4/pi;
% ATTENTION: in the first version S had signied switched. Corrected on
% April 18. After that it is wiht the proper sign
clear dummyx dummy dummyz

end

[xx yy]=meshgrid(1:ny,1:nx);

xxcut=xx(5:end-5,5:end-5)*dx;
yycut=yy(5:end-5,5:end-5)*dy;

close all
figure(1)
subplot(4,1,1)
pcolor(xxcut,yycut,B2cut)
shading interp
colorbar
subplot(4,1,2)
pcolor(xxcut,yycut,JEecut)
shading interp
colorbar
subplot(4,1,3)
pcolor(xxcut,yycut,JEicut)
shading interp
colorbar
subplot(4,1,4)
pcolor(xxcut,yycut,divScut)
shading interp
colorbar
set(gcf, 'Renderer', 'zbuffer');
print('-depsc','-r300',['energy_terms_xy_' ncyc '.eps'])
%print('-dpng','-r300',['energy_terms_xy_' ncyc '.png'])

xxx=mean(xxcut);
figure(2)
B2int=sum(B2cut);
plot(xxx,B2int)
hold on 
JEeint=sum(JEecut);
plot(xxx,JEeint,'r')
JEiint=sum(JEicut);
plot(xxx,JEiint,'m')
divSint=sum(divScut);
plot(xxx,divSint,'k')

%plot(xxx,-B2int+JEeint+JEiint,'g:')
plot(xxx,divSint-B2int+JEeint+JEiint,'y')
grid on
legend('dB^2dt','Je.E','Ji.E','divS','dE^2dt','Location','BestOutside')
print('-depsc',['energy_terms_x_' ncyc '.eps'])


figure(3)
subplot(2,1,1)
plot(xxx(round(end/2):end),cumsum(B2int(round(end/2):end)))
hold on
plot(xxx(round(end/2):end),cumsum(JEeint(round(end/2):end)),'r')
plot(xxx(round(end/2):end),cumsum(JEiint(round(end/2):end)),'m')
plot(xxx(round(end/2):end),cumsum(divSint(round(end/2):end)),'k')
legend('dB^2dt','Je.E','Ji.E','divS','location','BestOutside')
subplot(2,1,2)
plot(xxx(round(end/2):-1:1),cumsum(B2int(round(end/2):-1:1)))
hold on
plot(xxx(round(end/2):-1:1),cumsum(JEeint(round(end/2):-1:1)),'r')
plot(xxx(round(end/2):-1:1),cumsum(JEiint(round(end/2):-1:1)),'m')
plot(xxx(round(end/2):-1:1),cumsum(divSint(round(end/2):-1:1)),'k')
legend('dB^2dt','Je.E','Ji.E','divS','location','BestOutside')
print('-depsc',['energy_terms_cumsum_' ncyc '.eps'])

end
