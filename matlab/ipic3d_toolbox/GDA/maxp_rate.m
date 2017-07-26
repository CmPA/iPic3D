%maxp_common


Ezhist=[]
Ezloc=[];
Ezlocmin=[];
Ezlocmax=[];
Ezfluct=[];
Ezfluct2=[];
Byhist=[]

for cycle=Ncyc_ini:1000:Ncyc_max


ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

read=1
if(read)

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);

end

Ygsm = 6;
Xgsm = -35;
Ygsm2 = 5;
izrec = Nz-round(Nz*(9-Ygsm)/12)
izrec2 = Nz-round(Nz*(9-Ygsm2)/12)
ixrec = round(-Nx*(15+Xgsm)/30)

[nx ny nz]= size(Pper1)
y=1:ny;
for i=1:nx
for k=1:nz
w=Pper1(i,:,k);
%w=1.0./sqrt(Bx(i,:,k).^2+1e-10);

ymax(i,k)=round(sum(w.*y)./sum(w));
%ymax(i,k)=round(sum(y./V(i,:,k).^2)./sum(1./V(i,:,k).^2));
%[dum j] = min(Vx(i,:,k).^2);
%ymax(i,k)=j;
Bxmax(i,k)=Bx(i,round(ymax(i,k)),k);
%Vxmax(i,k)=Vx(i,end/2,k);
Bymax(i,k)=By(i,round(ymax(i,k)),k);
Bzmax(i,k)=Bz(i,round(ymax(i,k)),k);
Ezmax(i,k)=Ez(i,round(ymax(i,k)),k);
%Pparmax(i,k)=Pxx(i,end/2,k);
end
end

Ezhist=[Ezhist mean(Ezmax(:,izrec2-10:izrec2+10),2)];

Byhist=[Byhist mean(Bymax,2)];

end


Ygsm = 6;
Xgsm = -28;
Zgsm = -2;


[xd yd zd ixrec iyrec izrec ipx ipy ipz ip] = gsm2code(Xgsm, Ygsm, Zgsm);
Zgsm=gsmy2z(ymax(ixrec,iyrec)*dy)
[xd yd zd ixrec iyrec izrec ipx ipy ipz ip] = gsm2code(Xgsm, Ygsm, Zgsm);

ez=get_probe_ez(ip);
rate_local1=tsmovavg(ez(1,1:48000),'s',300)./B0^2;
time=60*([1:48000]/75000.0)*4;
figure(1)

plot(time,rate_local1,'r')
hold on
rate_local=tsmovavg(ez(1,1:48000),'s',3000)./B0^2;
plot(time,rate_local,'b')


cycle=Ncyc_ini:1000:Ncyc_max;
time2=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D;
plot(time2,max(Ezhist(ixrec-50:ixrec+50,:))./B0^2,'g')

figure
Ygsm = 6;
Xgsm = -35;
Ygsm2 = 5;
izrec = Nz-round(Nz*(9-Ygsm)/12)
izrec2 = Nz-round(Nz*(9-Ygsm2)/12)
ixrec = round(-Nx*(15+Xgsm)/30)
rate_local=rate_local1-rate_local+ interp1(time2,max(Ezhist(ixrec-50:ixrec+50,:))./B0^2,time);
plot(time,rate_local,'b',time,ones(size(time))*.1,'r')
hold on
plot(time2,max(Ezhist(ixrec-50:ixrec+50,:))./B0^2,'g')

ylim([-0.05 0.4])
set(gca,'fontsize',14)
xlabel('t(s)','fontsize',14)
ylabel('E_{rec}/B_0/V_{A0}','fontsize',14)
set(gcf, 'Renderer', 'zbuffer');
print('-depsc','rec_rate_3d')
%saveas('rec_rate_3d')
