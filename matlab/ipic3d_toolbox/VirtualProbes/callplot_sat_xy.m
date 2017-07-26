
plot_sat((1:n1)*Dt*B0x,xxp,EEX,yp(isaty),'EX')

clear EEX

plot_sat((1:n1)*Dt*B0x,xxp,EEY,yp(isaty),'EY')
clear EEY

plot_sat((1:n1)*Dt*B0x,xxp,EEZ,yp(isaty),'EZ')
clear EEZ

plot_sat((1:n1)*Dt*B0x,xxp,BBX,yp(isaty),'BX')

plot_sat((1:n1)*Dt*B0x,xxp,BBY,yp(isaty),'BY')

plot_sat((1:n1)*Dt*B0x,xxp,BBZ,yp(isaty),'BZ')


Navg=100
dBBX=diff(BBX,1,2)/Dt;
dBBXavg=tsmovavg(dBBX(:,1:npmax-1),'s',Navg);
plot_sat((1.5:n1-.5)*Dt*B0x,xxp,dBBX,yp(isaty),'dBXdt')
clear dBBX


dBBY=diff(BBY,1,2)/Dt;
dBBYavg=tsmovavg(dBBY(:,1:npmax-1),'s',Navg);
plot_sat((1.5:n1-.5)*Dt*B0x,xxp,dBBY,yp(isaty),'dBYdt')
clear dBBY

dBBZ=diff(BBZ,1,2)/Dt;
dBBZavg=tsmovavg(dBBZ(:,1:npmax-1),'s',Navg);
plot_sat((1.5:n1-.5)*Dt*B0x,xxp,dBBZ,yp(isaty),'dBZdt')
clear dBBZ

plot_sat((1:n1)*Dt*B0x,xxp,NE,yp(isaty),'NE')

figure(8)
pcolor((1:n1)*Dt*B0x,xxp,JX)
shading interp
colorbar
title(['u_{xe}   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEXsat_yz.png' )
saveas(gcf,['UEXsat_yz.fig'])

[sx,sy] = meshgrid(6.5,1:2:40);   
ut=ones(size(JY));
h=streamline(stream2((1:n1)*Dt*B0x,xxp,ut*Dt*B0x,JX,sx,sy,[1]));
set(h,'Color','white');

set(gcf,'Renderer','zbuffer');
print('-dpng','Flowtraces_xy.png' )
clear JX

close(8)

plot_sat((1:n1)*Dt*B0x,xxp,JY,yp(isaty),'UEY')
clear JY
plot_sat((1:n1)*Dt*B0x,xxp,JZ,yp(isaty),'UEZ')
clear JZ

plot_sat((1:n1)*Dt*B0x,xxp,sqrt((BBX.^2+BBY.^2+BBZ.^2)./abs(NE))./B0x,yp(isaty),'VA')

plot_sat((1:n1)*Dt*B0x,xxp,BBZ./BBX,yp(isaty),'BZoBX')
clear BBZ
plot_sat((1:n1)*Dt*B0x,xxp,BBY./BBX,yp(isaty),'BYoBX')
clear BBX
clear BBY

plot_sat((1:n1)*Dt*B0x,xxp,SSX,yp(isaty),'SX')
clear SSX
plot_sat((1:n1)*Dt*B0x,xxp,SSY,yp(isaty),'SY')
clear SSY
plot_sat((1:n1)*Dt*B0x,xxp,SSZ,yp(isaty),'SZ')
clear SSZ

plot_sat((Navg+1:n1-1)*Dt*B0x,xxp,dBBXavg(:,Navg+1:n1-1),yp(isaty),'dBXavg')
clear dBBXavg
plot_sat((Navg+1:n1-1)*Dt*B0x,xxp,dBBYavg(:,Navg+1:n1-1),yp(isaty),'dBYavg')
clear dBBYavg
plot_sat((Navg+1:n1-1)*Dt*B0x,xxp,dBBZavg(:,Navg+1:n1-1),yp(isaty),'dBZavg')
clear dBBZavg



if(detrended)
plot_sat((Ndetrend+1:n1)*Dt*B0x,xxp,SSXdet(:,Ndetrend+1:n1),yp(isaty),'SXdet')
clear SSXdet
plot_sat((2*Ndetrend+1:n1)*Dt*B0x,xxp,SSXdetavg(:,2*Ndetrend+1:n1),yp(isaty),'SXdetavg')
clear SSXdetavg
end

plot_sat((1:n1)*Dt*B0x,xxp,EDOTJ,yp(isaty),'EDOTJ')


Navg=100;
EDOTJ=tsmovavg(EDOTJ(:,1:npmax),'s',Navg);
EDOTJ=(2*EDOTJ(2:end-1,:)+EDOTJ(1:end-2,:)+EDOTJ(3:end,:))/4;

plot_sat((2*Navg+1:npmax-1)*Dt*B0x,xxp,EDOTJ(:,2*Navg+1:npmax-1),yp(isaty),'EDOTJavg')

