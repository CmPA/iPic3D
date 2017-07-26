function [ez]=get_probe_ez(ip)

results_dir='/data1/gianni/gda/HRmaha3D3/';
nome=[results_dir '1_VirtualSatelliteTraces' num2str(ip) '.txt']
%system(['gunzip ' nome])
if(exist(nome)~=2)
system(['scp giovanni@128.97.76.22:/rd/giovanni/HRmaha3d3/VS/data1/' 'VirtualSatelliteTraces' num2str(ip) '.txt' ' ' nome ])
end

fid=fopen(nome);
for i=1:27
x=fscanf(fid,'%f',3);
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
end


a=fscanf(fid,'%f',[16 inf])';

fclose(fid);

nome=[results_dir '2_VirtualSatelliteTraces' num2str(ip) '.txt']
%system(['gunzip ' nome])
if(exist(nome)~=2)
system(['scp giovanni@128.97.76.22:/rd/giovanni/HRmaha3d3/VS/data2/' 'VirtualSatelliteTraces' num2str(ip) '.txt' ' ' nome])
end
fid=fopen(nome);
for i=1:27
x=fscanf(fid,'%f',3);
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
end

b= fscanf(fid,'%f',[16 inf])';

size_a=size(a)
size_b=size(b)
aa=a;
a=a(1:14000*27,:);

a=[a ; b];

fclose(fid);


skip=0;
bx=a(:,3+skip);
by=a(:,4+skip);
bz=a(:,5+skip);
ex=a(:,6+skip);
ey=a(:,7+skip);
ez=a(:,8+skip);
jxe=a(:,9+skip);
jye=a(:,10+skip);
jze=a(:,11+skip);
jxi=a(:,12+skip);
jyi=a(:,13+skip);
jzi=a(:,14+skip);
rhoe=a(:,15+skip)*4*pi;
rhoi=a(:,16+skip)*4*pi;

[n m]=size(bx);
n1=floor(n/27)
%n1=14000;
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
t=linspace(0,n1,n1);
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;
