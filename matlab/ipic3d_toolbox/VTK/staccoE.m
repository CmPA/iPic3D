clear all
dir='/shared02/gianni/maha2/xt8'
file=[dir '/data1/E_xt.vtk']
[B,Bx,By,Bz]=read_vtk(file);    
file=[dir '/data2/E_xt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd]; 
Bx=[Bx;Bxadd]; 
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data3/E_xt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd];
Bx=[Bx;Bxadd];
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data4/E_xt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd];
Bx=[Bx;Bxadd];
By=[By;Byadd];
Bz=[Bz;Bzadd];

t=[0 60];
x=[-15 -45];

code_E = 2060.21;
code_B = 6.87213e-06;
code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;
e=1.6e-19;

imagesc(x,t,Bx*code_E)                
title('E_x','fontsize',14)
xlabel('x/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','ex.eps')
print('-dpng','-r300','ex.png')

imagesc(x,t,By*code_E)                
title('E_z','fontsize',14)
xlabel('x/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','ey.eps')
print('-dpng','-r300','ey.png')

imagesc(x,t,Bz*code_E)                
title('E_y','fontsize',14)
xlabel('x/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','ez.eps')
print('-dpng','-r300','ez.png')
