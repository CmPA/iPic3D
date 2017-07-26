dir='/nobackup/glapenta/1mar08C/data2/'
preamble

Xgsm=-16
ipxc= round(XLEN * (Xgsmrange(2)-Xgsm) / (Xgsmrange(2)-Xgsmrange(1)))+1
Ygsm=12
ipzc = round(ZLEN * (Ygsm-Ygsmrange(1)) / (Ygsmrange(2)-Ygsmrange(1)) )+1
Zgsm=1
ipyc = round(YLEN * (Zgsm-Zgsmrange(1)) / (Zgsmrange(2)-Zgsmrange(1)) )+1

for ipx=max(1,ipxc-3):min(ipxc+3,XLEN)
for ipy=ipyc:ipyc
for ipz=ipzc-1:ipzc+1

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;


nome=[dir 'VirtualSatelliteTraces' num2str(ip) '.txt']

if(exist(nome)~=2)
disp(['transfer ' nome])
return
end




fid=fopen(nome);

ennes=fscanf(fid,'%f',4)
nsat=ennes(1)
nsatx=ennes(2);
nsaty=ennes(3);
nsatz=ennes(4);

for i=1:27
x=fscanf(fid,'%f',3); 
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
end

fscanf(fid,'%f');
a=fread(fid,'float');

fclose(fid)
n=max(size(a));

aaa=a;

a=reshape(a,16,n/16)';

skip=0;
time=a(:,1+skip);
bx=a(:,3+skip);
by=a(:,4+skip);
bz=a(:,5+skip);
ex=a(:,6+skip);
ey=a(:,7+skip);
ez=a(:,8+skip);
rhoe=a(:,9+skip)*4*pi;
jxe=a(:,10+skip);
jye=a(:,11+skip);
jze=a(:,12+skip);
rhoi=a(:,13+skip)*4*pi;
jxi=a(:,14+skip);
jyi=a(:,15+skip);
jzi=a(:,16+skip);


mratio=256;

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

[n m]=size(bx);
n2=floor(n/27)
cycle=reshape(time(1:n2*27),27,n2);
time=60*(cycle(1,:)/75000.0*Dt/.125);

t=(time+initial_time)/86400;

bx=reshape(bx(1:n2*27),27,n2);
by=reshape(by(1:n2*27),27,n2);
bz=reshape(bz(1:n2*27),27,n2);
ex=reshape(ex(1:n2*27),27,n2);
ey=reshape(ey(1:n2*27),27,n2);
ez=reshape(ez(1:n2*27),27,n2);
rhoe=reshape(rhoe(1:n2*27),27,n2);
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;

for isat=1:27
plot(t,-ex(isat,:)*code_E,'k',t,ey(isat,:)*code_E,'r',t,ez(isat,:)*code_E,'g')
pbaspect([5 1 1])
datetick('x','HH:MM:SS')
legend('Ex','Ez','Ey','location','northwest')
title(['X=' num2str(gsmx(xp(isat))) '    Y=' num2str(gsmz2y(zp(isat))) '    Z=' num2str(gsmy2z(yp(isat)))])
print('-dpng',['E_VP_' num2str(ip) '_' num2str(isat) '.png'])

plot(t,-bx(isat,:),'k',t,by(isat,:),'r',t,bz(isat,:),'g')
pbaspect([5 1 1])
datetick('x','HH:MM:SS')
legend('Bx','Bz','By','location','northwest')
title(['X=' num2str(gsmx(xp(isat))) '    Y=' num2str(gsmz2y(zp(isat))) '    Z=' num2str(gsmy2z(yp(isat)))])
print('-dpng',['B_VP_' num2str(ip) '_' num2str(isat) '.png'])
end

end
end
end