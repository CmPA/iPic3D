clear all
dir='/shared02/gianni/maha2/yt8'
file=[dir '/data1/B_yt.vtk']
[B,Bx,By,Bz]=read_vtk(file);    
file=[dir '/data2/B_yt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd]; 
Bx=[Bx;Bxadd]; 
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data3/B_yt.vtk']
[Badd,Bxadd,Byadd,Bzadd]=read_vtk(file);
B=[B;Badd];
Bx=[Bx;Bxadd];
By=[By;Byadd];
Bz=[Bz;Bzadd];
file=[dir '/data4/B_yt.vtk']
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

grafic=1;
if(grafic) 
imagesc(x,t,Bx*code_B)                
title('B_x','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','bxYeps')
print('-dpng','-r300','bxYpng')

imagesc(x,t,By*code_B)                
title('B_z','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','byYeps')
print('-dpng','-r300','byYpng')

imagesc(x,t,Bz*code_B)                
title('B_y','fontsize',14)
xlabel('z/R_E','fontsize',14)
ylabel('t(s)','fontsize',14)
colorbar
pbaspect([.5  1 1])
%axis image
axis xy
set(gca,'xdir','reverse','TickDir','out')
print('-depsc','-r300','bzYeps')
print('-dpng','-r300','bzYpng')
end

fftok=0
if(fftok)
byfft=fft2(By);
[n1 n2]=size(byfft);
byfftg=byfft(1:round(n1/2)+1,1:round(n2/2)+1);
imagesc(log(abs(byfftg)))
print -depsc fftBy
end
