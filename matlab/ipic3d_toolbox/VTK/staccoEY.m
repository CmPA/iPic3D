clear all
dir='/shared02/gianni/maha2/yt8'
file=[dir '/data1/E_yt.vtk']
[B,Bx,By,Bz]=read_vtk(file);    
file=[dir '/data2/E_yt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd]; 
Bx=[Bx;Bxadd]; 
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data3/E_yt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd];
Bx=[Bx;Bxadd];
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data4/E_yt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd];
Bx=[Bx;Bxadd];
By=[By;Byadd];
Bz=[Bz;Bzadd];

t=[0 60];
x=[3 -9];

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

grafic=1
if(grafic)
imagesc(x,t,Bx*code_E)                
title('E_x','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','exYeps')
print('-dpng','-r300','exYpng')

imagesc(x,t,By*code_E)                
title('E_z','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','eyYeps')
print('-dpng','-r300','eyYpng')

imagesc(x,t,Bz*code_E)                
title('E_y','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','ezYeps')
print('-dpng','-r300','ezYpng')
end

fftok=0
if(fftok)
byfft=fft2(By);
[n1 n2]=size(byfft);
byfftg=byfft(1:round(n1/2)+1,1:round(n2/2)+1);
imagesc(log(abs(byfftg)))
print -depsc fftEy
end
