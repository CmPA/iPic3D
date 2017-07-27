addpath('/home/gianni/matlab2/matlab-parsek');
global contours

close all

contours = 1;
for i=20000:1000:20000
disp(['CYCLE = ' num2str(i)])
ncyc=num2str(i);
read=1;

dire ='/shared02/gianni/tred60/data/'

if (read) 
   
    [V,JEex,JEey,JEez,dx,dy,dz]=read_vtk([dire 'AVGJdotEe_AVG_xy_cycle' ncyc '.vtk']);
    [V,dVx,dVy,dVz,dx,dy,dz]=read_vtk([dire 'AVGdelJdotEe_AVG_xy_cycle' ncyc '.vtk']);
JEex=JEex+dVx; JEey=JEey+dVy; JEez=JEez+dVz;
[nx ny]=size(JEex);
JEecut=JEex(5:end-5,5:end-5)+JEey(5:end-5,5:end-5)+JEez(5:end-5,5:end-5);
clear JEex JEey JEez

    [V,JEix,JEiy,JEiz,dx,dy,dz]=read_vtk([dire 'AVGJdotEi_AVG_xy_cycle' ncyc '.vtk']);
    [V,dVx,dVy,dVz,dx,dy,dz]=read_vtk([dire 'AVGdelJdotEi_AVG_xy_cycle' ncyc '.vtk']);
JEix=JEix+dVx; JEiy=JEiy+dVy; JEiz=JEiz+dVz;
JEicut=JEix(5:end-5,5:end-5)+JEiy(5:end-5,5:end-5)+JEiz(5:end-5,5:end-5);
clear JEix JEiy JEiz

[divQe,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGdivQe_AVG_xy_cycle' ncyc '.vtk']);
divQecut=divQe(5:end-5,5:end-5); 

[divQbulke,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGdivQbulke_AVG_xy_cycle' ncyc '.vtk']);
divQbulkecut=divQbulke(5:end-5,5:end-5); 

[vdivPe,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvdivPe_AVG_xy_cycle' ncyc '.vtk']);
vdivPecut=vdivPe(5:end-5,5:end-5);

[vgradUbulke,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvgradUbulke_AVG_xy_cycle' ncyc '.vtk']);
vgradUbulkecut=vgradUbulke(5:end-5,5:end-5);

[vgradUthe,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvgradUthe_AVG_xy_cycle' ncyc '.vtk']);
vgradUthecut=vgradUthe(5:end-5,5:end-5);

[divQi,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGdivQi_AVG_xy_cycle' ncyc '.vtk']);
divQicut=divQi(5:end-5,5:end-5);

[divQbulki,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGdivQbulki_AVG_xy_cycle' ncyc '.vtk']);
divQbulkicut=divQbulki(5:end-5,5:end-5);

[vdivPi,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvdivPi_AVG_xy_cycle' ncyc '.vtk']);
vdivPicut=vdivPi(5:end-5,5:end-5);

[vgradUbulki,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvgradUbulki_AVG_xy_cycle' ncyc '.vtk']);
vgradUbulkicut=vgradUbulki(5:end-5,5:end-5);

[vgradUthi,dummyx,dummyy,dummyz,dx,dy,dz]=read_vtk([dire 'AVGvgradUthi_AVG_xy_cycle' ncyc '.vtk']);
vgradUthicut=vgradUthi(5:end-5,5:end-5);

clear dummyx dummy dummyz

end

[xx yy]=meshgrid(1:ny,1:nx);

xxcut=xx(5:end-5,5:end-5)*dx;
yycut=yy(5:end-5,5:end-5)*dy;

cuts2D=0
if(cuts2D)

close all
figure(1)
subplot(4,1,1)
pcolor(xxcut,yycut,JEecut)
shading interp
colorbar
subplot(4,1,2)
pcolor(xxcut,yycut,JEicut)
shading interp
colorbar
subplot(4,1,3)
pcolor(xxcut,yycut,divQecut)
shading interp
colorbar
subplot(4,1,4)
pcolor(xxcut,yycut,divQicut)
shading interp
colorbar
set(gcf, 'Renderer', 'zbuffer');
%print('-depsc','-r300',['parenergy_terms_xy_' ncyc '.eps'])
print('-dpng','-r300',['parenergy_terms_xy_' ncyc '.png'])

close all
figure(1)
subplot(4,1,1)
pcolor(xxcut,yycut,divQbulkecut)
shading interp
colorbar
subplot(4,1,2)
pcolor(xxcut,yycut,divQbulkicut)
shading interp
colorbar
subplot(4,1,3)
pcolor(xxcut,yycut,vdivPecut)
shading interp
colorbar
subplot(4,1,4)
pcolor(xxcut,yycut,vdivPicut)
shading interp
colorbar
set(gcf, 'Renderer', 'zbuffer');
print('-dpng','-r300',['parenergy2_terms_xy_' ncyc '.png'])

end

xxx=mean(xxcut);

close all
figure(2)
JEeint=sum(JEecut);

divQbulkeint=sum(divQbulkecut);

vdivPeint=sum(vdivPecut);

divQeint=sum(divQecut);

vgradUbulkeint=sum(vgradUbulkecut);

vgradUtheint=sum(vgradUthecut);


euler_bal_bulk_e = -vdivPeint-divQbulkeint+JEeint;
euler_bal_th_e = vdivPeint+divQbulkeint-divQeint;

lagrange_bal_bulk_e = vgradUbulkeint-vdivPeint-divQbulkeint+JEeint;
lagrange_bal_th_e = vgradUtheint+vdivPeint+divQbulkeint-divQeint;

% bulk balance electrons
plot(xxx,JEeint,'r')
hold on
plot(xxx,divQbulkeint,'m')
plot(xxx,vdivPeint,'g')
plot(xxx,vgradUbulkeint,'y')
plot(xxx,euler_bal_bulk_e,'b')
plot(xxx,lagrange_bal_bulk_e,'k')
ylim([-1e-7 2e-7])
grid on
legend('J.E','divQbulk','v\nabla P','vgradUbulk','balance-Euler','balance-Lagrange','Location','BestOutside')
title('Bulk Energy Electrons')
print('-depsc',['parenergy_terms_x_bulk_e_' ncyc '.eps'])


close all
% bulk balance electrons -smooth
plot(xxx,smooth(JEeint),'r')
hold on
plot(xxx,smooth(divQbulkeint),'m')
plot(xxx,smooth(vdivPeint),'g')
plot(xxx,smooth(vgradUbulkeint),'y')
plot(xxx,smooth(euler_bal_bulk_e),'b')
plot(xxx,smooth(lagrange_bal_bulk_e),'k')
ylim([-1e-7 2e-7])
grid on
legend('J.E','divQbulk','v\nabla P','vgradUbulk','balance-Euler','balance-Lagrange','Location','BestOutside')
title('Bulk Energy Electrons')
print('-depsc',['parenergy_terms_x_bulk_smooth_e_' ncyc '.eps'])

close all
% thermal balance electrons
plot(xxx,divQeint-divQbulkeint,'m')
hold on 
plot(xxx,vdivPeint,'g')
plot(xxx,vgradUtheint,'y')
plot(xxx,euler_bal_th_e,'b')
plot(xxx,lagrange_bal_th_e,'k')
ylim([-.4e-6 .4e-6])
grid on
title('Thermal Energy Electrons')
legend('divQenth','v\nabla P','vgradUth','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_th_e_' ncyc '.eps'])

close all
% thermal balance electrons -smooth
plot(xxx,smooth(divQeint-divQbulkeint),'m')
hold on 
plot(xxx,smooth(vdivPeint),'g')
plot(xxx,smooth(vgradUtheint),'y')
plot(xxx,smooth(euler_bal_th_e),'b')
plot(xxx,smooth(lagrange_bal_th_e),'k')
ylim([-.4e-6 .4e-6])
grid on
title('Thermal Energy Electrons')
legend('divQenth','v\nabla P','vgradUth','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_th__smooth_e_' ncyc '.eps'])


close all
figure(2)
JEiint=sum(JEicut);

divQbulkiint=sum(divQbulkicut);

vdivPiint=sum(vdivPicut);

divQiint=sum(divQicut);

vgradUbulkiint=sum(vgradUbulkicut);

vgradUthiint=sum(vgradUthicut);


euler_bal_bulk_i = -vdivPiint-divQbulkiint+JEiint;
euler_bal_th_i = vdivPiint+divQbulkiint-divQiint;

lagrange_bal_bulk_i = vgradUbulkiint-vdivPiint-divQbulkiint+JEiint;
lagrange_bal_th_i = vgradUthiint+vdivPiint+divQbulkiint-divQiint;

close all
% bulk balance ions
plot(xxx,JEiint,'r')
hold on
plot(xxx,divQbulkiint,'m')
plot(xxx,vdivPiint,'g')
plot(xxx,vgradUbulkiint,'y')
plot(xxx,euler_bal_bulk_i,'b')
plot(xxx,lagrange_bal_bulk_i,'k')
ylim([-1e-6 2e-6])
grid on
title('Bulk Energy Ions')
legend('Ji.E','divQbulk','v\nabla P','vgradUbulk','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_bulk_i_' ncyc '.eps'])


close all
% thermal balance ions
plot(xxx,divQiint-divQbulkiint,'m')
hold on 
plot(xxx,vdivPiint,'g')
plot(xxx,vgradUthiint,'y')
plot(xxx,euler_bal_th_i,'b')
plot(xxx,lagrange_bal_th_i,'k')
ylim([-3e-6 3e-6])
grid on
title('Thermal Energy Ions')
legend('divQenth','v\nabla P','vgradUth','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_th_i_' ncyc '.eps'])


% The cumulative sum below needs to be updated for the split bilk/thermal
close all
figure(3)
subplot(2,1,1)
plot(xxx(round(end/2):end),cumsum(JEeint(round(end/2):end)),'r')
hold on
plot(xxx(round(end/2):end),cumsum(divQbulkeint(round(end/2):end)),'m')
plot(xxx(round(end/2):end),cumsum(vdivPeint(round(end/2):end)),'g')
plot(xxx(round(end/2):end),cumsum(vgradUbulkeint(round(end/2):end)),'y')
plot(xxx(round(end/2):end),cumsum(euler_bal_bulk_e(round(end/2):end)),'b')
plot(xxx(round(end/2):end),cumsum(lagrange_bal_bulk_e(round(end/2):end)),'k')
subplot(2,1,2)
plot(xxx(round(end/2):-1:1),cumsum(JEeint(round(end/2):-1:1)),'r')
hold on
plot(xxx(round(end/2):-1:1),cumsum(divQbulkeint(round(end/2):-1:1)),'m')
plot(xxx(round(end/2):-1:1),cumsum(vdivPeint(round(end/2):-1:1)),'g')
plot(xxx(round(end/2):-1:1),cumsum(vgradUbulkeint(round(end/2):-1:1)),'y')
plot(xxx(round(end/2):-1:1),cumsum(euler_bal_bulk_e(round(end/2):-1:1)),'b')
plot(xxx(round(end/2):-1:1),cumsum(lagrange_bal_bulk_e(round(end/2):-1:1)),'k')
title('Bulk Energy Electrons')
legend('J.E','divQbulk','v\nabla P','vgradUbulk','balance-Euler','balance-Lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_bulk_cumsum_e_' ncyc '.eps'])

close all
figure(3)
subplot(2,1,1)
plot(xxx(round(end/2):end),cumsum(JEiint(round(end/2):end)),'r')
hold on
plot(xxx(round(end/2):end),cumsum(divQbulkiint(round(end/2):end)),'m')
plot(xxx(round(end/2):end),cumsum(vdivPiint(round(end/2):end)),'g')
plot(xxx(round(end/2):end),cumsum(vgradUbulkiint(round(end/2):end)),'y')
plot(xxx(round(end/2):end),cumsum(euler_bal_bulk_i(round(end/2):end)),'b')
plot(xxx(round(end/2):end),cumsum(lagrange_bal_bulk_i(round(end/2):end)),'k')
subplot(2,1,2)
plot(xxx(round(end/2):-1:1),cumsum(JEiint(round(end/2):-1:1)),'r')
hold on
plot(xxx(round(end/2):-1:1),cumsum(divQbulkiint(round(end/2):-1:1)),'m')
plot(xxx(round(end/2):-1:1),cumsum(vdivPiint(round(end/2):-1:1)),'g')
plot(xxx(round(end/2):-1:1),cumsum(vgradUbulkiint(round(end/2):-1:1)),'y')
plot(xxx(round(end/2):-1:1),cumsum(euler_bal_bulk_i(round(end/2):-1:1)),'b')
plot(xxx(round(end/2):-1:1),cumsum(lagrange_bal_bulk_i(round(end/2):-1:1)),'k')
title('Bulk Energy Ions')
legend('J.E','divQbulk','v\nabla P','vgradUbulk','balance-Euler','balance-Lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_bulk_cumsum_i_' ncyc '.eps'])

close all
figure(3)
subplot(2,1,1)
plot(xxx(round(end/2):end),cumsum(divQeint(round(end/2):end)-divQbulkeint(round(end/2):end)),'m')
hold on 
plot(xxx(round(end/2):end),cumsum(vdivPeint(round(end/2):end)),'g')
plot(xxx(round(end/2):end),cumsum(vgradUtheint(round(end/2):end)),'y')
plot(xxx(round(end/2):end),cumsum(euler_bal_th_e(round(end/2):end)),'b')
plot(xxx(round(end/2):end),cumsum(lagrange_bal_th_e(round(end/2):end)),'k')
subplot(2,1,2)
plot(xxx(round(end/2):-1:1),cumsum(divQeint(round(end/2):-1:1)-divQbulkeint(round(end/2):-1:1)),'m')
hold on 
plot(xxx(round(end/2):-1:1),cumsum(vdivPeint(round(end/2):-1:1)),'g')
plot(xxx(round(end/2):-1:1),cumsum(vgradUtheint(round(end/2):-1:1)),'y')
plot(xxx(round(end/2):-1:1),cumsum(euler_bal_th_e(round(end/2):-1:1)),'b')
plot(xxx(round(end/2):-1:1),cumsum(lagrange_bal_th_e(round(end/2):-1:1)),'k')
title('Thermal Energy Electrons')
legend('divQenth','v\nabla P','vgradUth','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_th_cumsum_e_' ncyc '.eps'])

close all
figure(3)
subplot(2,1,1)
plot(xxx(round(end/2):end),cumsum(divQiint(round(end/2):end)-divQbulkiint(round(end/2):end)),'m')
hold on 
plot(xxx(round(end/2):end),cumsum(vdivPiint(round(end/2):end)),'g')
plot(xxx(round(end/2):end),cumsum(vgradUtheint(round(end/2):end)),'y')
plot(xxx(round(end/2):end),cumsum(euler_bal_th_i(round(end/2):end)),'b')
plot(xxx(round(end/2):end),cumsum(lagrange_bal_th_i(round(end/2):end)),'k')
subplot(2,1,2)
plot(xxx(round(end/2):-1:1),cumsum(divQiint(round(end/2):-1:1)-divQbulkiint(round(end/2):-1:1)),'m')
hold on 
plot(xxx(round(end/2):-1:1),cumsum(vdivPiint(round(end/2):-1:1)),'g')
plot(xxx(round(end/2):-1:1),cumsum(vgradUthiint(round(end/2):-1:1)),'y')
plot(xxx(round(end/2):-1:1),cumsum(euler_bal_th_i(round(end/2):-1:1)),'b')
plot(xxx(round(end/2):-1:1),cumsum(lagrange_bal_th_i(round(end/2):-1:1)),'k')
title('Thermal Energy Electrons')
legend('divQenth','v\nabla P','vgradUth','balance-euler','balance-lagrange','Location','BestOutside')
print('-depsc',['parenergy_terms_x_th_cumsum_i_' ncyc '.eps'])
end
