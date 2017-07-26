addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

% This works i the run has saved the particles
%results_dir='/shared/gianni/lhdi6redux/';
%results_dir='/shared02/gianni/tred60/data/';
results_dir='/shared/gianni/laila11b/';
results_dir='/shared02/gianni/lailaO1/data5/';
Nspecies=4
species_plot=3
results_dir='/shared02/gianni/maha2/data4/';
Nspecies=2
species_plot=2



filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
ZLEN=1
%ZLEN=double(hdf5read(filename,'/topology/ZLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
Lz=1
%Lz=double(hdf5read(filename,'/collective/Lz'))
%B0x=double(hdf5read(filename,'/collective/B0x'))
B0x=double(hdf5read(filename,'/collective/Bx0'))
Nxc=double(hdf5read(filename,'/collective/Nxc'))
Nyc=double(hdf5read(filename,'/collective/Nyc'))
Nzc=1
%Nzc=double(hdf5read(filename,'/collective/Nzc'))
if (species_plot==1) 
vth=double(hdf5read(filename,'/collective/species_0/uth'))
else
vth=double(hdf5read(filename,'/collective/species_1/uth'))
end
% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))

mp = 1.6726e-27;
me= mp/256;
e= 1.6022e-19;
mu0=4*pi*1e-7;
eps0=8.8542*1.e-12;
c=1/sqrt(mu0*eps0);
TeoTi=1/5;
energy_code=me*c^2/e/TeoTi;

bins=0:(vth/5).^2:(5*vth)^2;
bins=bins*energy_code;
histpar=[];
histperp=[];
xhist=[];
bhist=[];
bxhist=[];
byhist=[];
bzhist=[];

ipx=XLEN/2;
for  ipx=1:XLEN
%ipy=YLEN/2+2;
%ipy=YLEN/2-1*3;
ipy=YLEN/2+1;
if(ZLEN>1)
ipz=ZLEN/2;
else
ipz=1;
end
%ipy=10

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

xmin=double(0+(ipx-1)*Lx/XLEN)
xmax=double(xmin+Lx/XLEN)
ymin=double(0+(ipy-1)*Ly/YLEN)
ymax=double(ymin+double(Ly/YLEN))
zmin=double(0+(ipz-1)*Lz/ZLEN)
zmax=double(zmin+double(Lz/ZLEN))

nx=Nxc/XLEN+1;
ny=Nyc/YLEN+1;
nz=Nzc/ZLEN+1;

dx=double(xmax-xmin)/double(nx-1);
dy=double(ymax-ymin)/double(ny-1);
dz=double(zmax-zmin)/double(nz-1);

[yg,xg]=ndgrid(ymin:dy:ymax,xmin:dx:xmax);


close all

vmax_graph=[.25 .06 .25 .06]/3;
for is=2:Nspecies

vmax = vmax_graph(is);
partfile='part';
partfile='restart';
if exist([results_dir partfile num2str(ip) '.hdf']) ==2 

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


%infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

%[ntf dummy]=size(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(:));

for it=nt:nt
%ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(it)));   
%ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(it)));   


Bx=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(1).Datasets(it)));
By=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(2).Datasets(it)));
Bz=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(3).Datasets(it)));
% for runs with ipic3D
Bx=squeeze(Bx(1,:,:));
By=squeeze(By(1,:,:));
Bz=squeeze(Bz(1,:,:));


%
% 1 = q
% 2 = u
% 3 = v
% 4 = w
% 5 = x
% 6 = y
% 7 = z

q=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(1).Datasets(it)));
x=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(it)));
y=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(6).Datasets(it)));
%z=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(7).Datasets(it)));
bxp=interp2(xg,yg,Bx,x,y);
byp=interp2(xg,yg,By,x,y);
bzp=interp2(xg,yg,Bz,x,y);
bp=sqrt(bxp.^2+byp.^2+bzp.^2);
bp2D=sqrt(bxp.^2+byp.^2);
perp2x=bzp.*bxp./(bp.*bp2D);
perp2y=bzp.*byp./(bp.*bp2D);
perp2z=-bp2D./bp;

u=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(2).Datasets(it)));
v=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(3).Datasets(it)));
w=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(4).Datasets(it)));

upar=(u.*bxp+v.*byp+w.*bzp)./(bp+1e-10);
uperp1=(byp.*u-bxp.*v)./bp2D;

uperp2=perp2x.*u+perp2y.*v+perp2z.*w;

epar=.5*upar.^2;
eperp=.5*(uperp1.^2+uperp2.^2);

epar=epar*energy_code;
eperp=eperp*energy_code;

if (is == species_plot) 
histpar=[histpar; hist(epar,bins)];
histperp=[histperp; hist(eperp,bins)]; 
xhist=[xhist;mean(x)];
bhist=[bhist;mean(bp)];
bxhist=[bxhist;mean(bxp)];
byhist=[byhist;mean(byp)];
bzhist=[bzhist;mean(bzp)];
end

max(y)
min(y)

%ii=y<7;
%q=q(ii);
%x=x(ii);
%y=y(ii);
%u=u(ii);
% q2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(2).Datasets(it)));
% x2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(6).Datasets(it)));
% u2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(3).Datasets(it)));
% q=[q;q2];
% x=[x;x2];
% u=[u;u2];
%figure(1)
%plot(x,y,'.','MarkerSize',1)
%axis square
%set(gcf,'Renderer','zbuffer');
%print('-dpng',['filename '.png'])
%print -dpng figure1


figure(1)
dbins=371;
RE=6371000;
box=1.6671*RE;
npart=32*32*14*14;
n0=.25e6;
c=300000e3;
[totnum,nbinu,xrange,urange]=spaziofasi(upar'/TeoTi,uperp1'/TeoTi,q',0);
dv=mean(diff(urange))*c;
dx=mean(diff(xrange))*c;
nbinunorm=nbinu/box^3/dv/dx/npart*box^3*n0;
%pcolor(xrange(1:end-1)',urange',log10(nbinunorm(1:end-1,:))')
pcolor(xrange(1:end-1)',urange',(nbinunorm(1:end-1,:))')
%lighting phong
shading interp
%camlight(0,90) % luce dall'alto
%view(2) %visione dall'alto
%xlim([-vmax vmax]/TeoTi)
%ylim([-vmax vmax]/TeoTi)	 
	 axis square 
%title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
title(['x(GSM)=' num2str(-(xmin+xmax)/2/Lx*30-15) '  z(GSM)=' num2str((ymin+ymax)/2/Ly*12-9) ])
xlabel('U_{||}/c','fontsize',[14])
         ylabel('U_{\perp 1}/c','fontsize',[14])
	 colorbar
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','-r600',['SP' num2str(is) 'phasespaceUV' num2str(ipx) '.png'])
print('-depsc','-r600',['SP' num2str(is) 'phasespaceUV' num2str(ipx) '.eps'])
saveas(gcf,['SP' num2str(is) 'phasespaceUV' num2str(ipx) '.fig']);
filesave=[num2str(is) 'UV' num2str(ipx)];
save(filesave,'xrange', 'urange','nbinu')

others=0
if(others)

figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(upar',uperp2',q',0);
surf(xrange(1:end-1)',urange',log10(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
%	 xlim([-vmax vmax])
%	 ylim([-vmax vmax])	 
	 axis square 
title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
	 xlabel('U_{||}/c','fontsize',[14])
	 ylabel('U_{\perp 2}/c','fontsize',[14])
	 colorbar
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','-r300',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.png'])
print('-depsc','-painters',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.eps'])


	 figure(1)
	 [totnum,nbinu,xrange,urange]=spaziofasi(uperp1',uperp2',q',0);
				surf(xrange(1:end-1)',urange',log10(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
					 lighting phong
					 shading interp
					 camlight(0,90) % luce dall'alto
					 view(2) %visione dall'alto
					% xlim([-vmax vmax])
					 %ylim([-vmax vmax])	 
					 axis square 
					 title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
					 xlabel('U_{\perp 1}/c','fontsize',[14])
					 ylabel('U_{\perp 2}/c','fontsize',[14])
					 colorbar
set(gca,'fontsize',[14])
					 set(gcf,'Renderer','zbuffer');
					 print('-dpng','-r300',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.png'])
					 print('-painters','-r300',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.eps'])
					 
end

% This below works only in 2D
%[nx ny]=size(ex);
%dx=double(xmax-xmin)/double(nx-1);
%dy=double(ymax-ymin)/double(ny-1);
%[xg,yg]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

%figure(3)
%pcolor(xg,yg,ex)
%shading interp
%set(gcf,'Renderer','zbuffer');
%print -dpng figure3
%pause
end

end

end

end

figure
subplot(3,1,1:2)
[n1 n2]=size(histperp)
for i=1:n1
	histperp2(i,:) = histperp(i,:).*bins;
end					 
xhistRe=-15-xhist/100*30;
%1e4 to go to cm^-2
dbins=371;
RE=6371000;
box=1.6671*RE;
npart=32*32*14*14;
n0=.25e6;
histperp2norm=histperp2/box^3/npart*(n0*box^3)*1e4/dbins;
pcolor(xhistRe,bins,log10(histperp2norm'))
ylabel('energy','fontsize',[14])
title('Perp Energy Flux Ef(E)','fontsize',[14])
colorbar('NorthOutside')
shading interp
set(gca,'fontsize',[14])
set(gca,'xdir','reverse','TickDir','out')
subplot(3,1,3)
code_B = 6.87213e-06;
code_B = code_B*1e9;
bhist=bhist*code_B;
bxhist=bxhist*code_B;
byhist=byhist*code_B;
bzhist=bzhist*code_B;
plot(xhistRe,bhist,xhistRe,-bxhist,xhistRe,bzhist,xhistRe,byhist)
axis tight
xlabel('x/R_E','fontsize',[14])
				legend('B','B_x','B_y','B_z')
set(gca,'fontsize',[14])
set(gca,'xdir','reverse','TickDir','out')
set(gcf,'Renderer','zbuffer');
print('-dpng','-r600','EnergyFluxPerp.png')
print('-depsc','-r600','EnergyFluxPerp.eps')
saveas(gcf,'EnergyFluxPerp.fig');
save EnFluxPerp xhist bins histperp2norm


figure
subplot(2,1,1)
for i=1:n1
	histpar2(i,:) = histpar(i,:).*bins;
end					 
pcolor(xhist,bins,log(histpar2'))
xlabel('x/d_i','fontsize',[14])
ylabel('energy','fontsize',[14])
title('Par Energy Flux Ef(E)','fontsize',[14])
shading interp
set(gca,'fontsize',[14])
subplot(2,1,2)
plot(xhist,bhist,xhist,-bxhist,xhist,byhist,xhist,bzhist)
					  legend('B','B_x','B_y','B_z')
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','EnergyFluxPar.png')
print('-deps','-painters','EnergyFluxPar.eps')


ii=[n1/4+3 n1/2-6 n1/2-4  n1/2-2];

figure
semilogy(bins,histpar(ii,:))
legend(num2str(xhist(ii)))
title('Parallel','fontsize',[14])
xlabel('Energy','fontsize',[14])
ylabel('pdf(E)','fontsize',[14]) 
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp.png')
print('-depsc','-r300','tmp.eps')


figure
alfa= flipud(cumsum(flipud(histpar')))';
for i=1:n1
alfa(i,:)=alfa(i,:)/alfa(i,1);
end

loglog(bins,alfa(ii,:))
%xlim([1e-4 1e-1])
legend(num2str(xhist(ii))) 
title('Parallel','fontsize',[14])
xlabel('Energy','fontsize',[14])
ylabel('cdf(E)','fontsize',[14]) 
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp2.png')
print('-depsc','-painters','tmp2.eps')

figure
semilogy(bins,histperp(ii,:))
legend(num2str(-15-30*xhist(ii)/100)) 
title('Perp','fontsize',[14])
xlabel('Energy','fontsize',[14])
ylabel('pdf(E)','fontsize',[14]) 
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp3.png')
print('-depsc','-painters','tmp3.eps')



figure
alfa= flipud(cumsum(flipud(histperp')))';
for i=1:n1
alfa(i,:)=alfa(i,:)/alfa(i,1);
end
loglog(bins,alfa(ii,:),bins.*byhist(ii(1))./byhist(ii(2)),alfa(ii(2),:),'--',bins,alfa(ii(2),:),'x')
%xlim([1e-4 1e-1])
legend(num2str(xhist(ii))) 
title('Perp','fontsize',[14])
xlabel('Energy','fontsize',[14])
ylabel('cdf(E)','fontsize',[14]) 
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp4.png')
print('-depsc','-painters','tmp4.eps')


close all
figure
semilogy(histperp2(:,1:30:end/2)) 
%legend(num2str(bins(1:30:end/2)'*96/.012),'Location','NorthOutside')
legend(num2str(bins(1:30:end/2)'),'Location','NorthOutside')
title('Energy in eV','fontsize',[14])
xlabel('x/d_i','fontsize',[14])
ylabel('f(E)','fontsize',[14])
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','pdf1.png')

print -dpng caz
