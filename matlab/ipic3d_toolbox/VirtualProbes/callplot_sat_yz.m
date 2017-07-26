if (exist('OHMX'))
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMX,mean(yyp),'OHM_x')
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMY,mean(yyp),'OHM_y')
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMZ,mean(yyp),'OHM_z')
end
%plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMiX,mean(yyp),'OHMi_x')
%plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMiY,mean(yyp),'OHMi_y')
%plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,OHMiZ,mean(yyp),'OHMi_z')

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,(EEX.*BBX+EEY.*BBY+EEZ.*BBZ)./sqrt(BBX.*BBX+BBY.*BBY+BBZ.*BBZ),mean(yyp),'Epar')

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,EEX,mean(yyp),'EY')
clear EEX

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,EEY,mean(yyp),'EY')
clear EEY

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,EEZ,mean(yyp),'EZ')
clear EEZ

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,BBX,mean(yyp),'BX')

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,BBY,mean(yyp),'BY')

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,BBZ,mean(yyp),'BZ')


Navg=100
  npmax=max(size(BBX));
dBBX=diff(BBX,1,2)/Dt;
dBBXavg=tsmovavg(dBBX(:,1:npmax-1),'s',Navg);
plot_sat((cycle0+.5:cycle0+n1-1.5)*Dt*B0x,xxp,dBBX,mean(yyp),'dBXdt')
clear dBBX


dBBY=diff(BBY,1,2)/Dt;
dBBYavg=tsmovavg(dBBY(:,1:npmax-1),'s',Navg);
plot_sat((cycle0+.5:cycle0+n1-1.5)*Dt*B0x,xxp,dBBY,mean(yyp),'dBYdt')
clear dBBY

dBBZ=diff(BBZ,1,2)/Dt;
dBBZavg=tsmovavg(dBBZ(:,1:npmax-1),'s',Navg);
plot_sat((cycle0+.5:cycle0+n1-1.5)*Dt*B0x,xxp,dBBZ,mean(yyp),'dBZdt')
clear dBBZ

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,NE,mean(yyp),'NE')

figure(8)
pcolor((cycle0:cycle0+n1-1)*Dt*B0x,xxp,JX)
shading interp
colorbar
title(['u_{xe}   y=   ' num2str(mean(yyp))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEXsat_yz.png' )
saveas(gcf,['UEXsat_yz.fig'])

[sx,sy] = meshgrid(6.5,1:2:40);   
ut=ones(size(JY));
h=streamline(stream2((cycle0:cycle0+n1-1)*Dt*B0x,xxp,ut*Dt*B0x,JX,sx,sy,[1]));
set(h,'Color','white');

set(gcf,'Renderer','zbuffer');
print('-dpng','Flowtraces_xy.png' )


close(8)

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,JY.*NE+JiX,mean(yyp),'JEX')
clear JX
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,JY.*NE+JiY,mean(yyp),'JEY')
clear JY
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,JZ.*NE+JiZ,mean(yyp),'JEZ')
clear JZ

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,sqrt((BBX.^2+BBY.^2+BBZ.^2)./abs(NE))./B0x,mean(yyp),'VA')

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,BBZ./BBX,mean(yyp),'BZoBX')
clear BBZ
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,BBY./BBX,mean(yyp),'BYoBX')
clear BBX
clear BBY

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,SSX,mean(yyp),'SX')
clear SSX
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,SSY,mean(yyp),'SY')
clear SSY
plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,SSZ,mean(yyp),'SZ')
clear SSZ

plot_sat((Navg+cycle0:cycle0+n1-1-1)*Dt*B0x,xxp,dBBXavg(:,Navg+1:n1-1),mean(yyp),'dBXavg')
clear dBBXavg
plot_sat((Navg+cycle0:cycle0+n1-1-1)*Dt*B0x,xxp,dBBYavg(:,Navg+1:n1-1),mean(yyp),'dBYavg')
clear dBBYavg
plot_sat((Navg+cycle0:cycle0+n1-1-1)*Dt*B0x,xxp,dBBZavg(:,Navg+1:n1-1),mean(yyp),'dBZavg')
clear dBBZavg



if(detrended)
plot_sat((Ndetrend+cycle0:cycle0+n1-1)*Dt*B0x,xxp,SSXdet(:,Ndetrend+cycle0:cycle0+n1-1),mean(yyp),'SXdet')
clear SSXdet
plot_sat((2*Ndetrend+cycle0:cycle0+n1-1)*Dt*B0x,xxp,SSXdetavg(:,2*Ndetrend+cycle0:cycle0+n1-1),mean(yyp),'SXdetavg')
clear SSXdetavg
end

plot_sat((cycle0:cycle0+n1-1)*Dt*B0x,xxp,EDOTJ,mean(yyp),'EDOTJ')


Navg=100;
EDOTJ=tsmovavg(EDOTJ(:,1:npmax),'s',Navg);
EDOTJ=(2*EDOTJ(2:end-1,:)+EDOTJ(1:end-2,:)+EDOTJ(3:end,:))/4;

plot_sat((2*Navg+1:npmax-1)*Dt*B0x,xxp,EDOTJ(:,2*Navg+1:npmax-1),mean(yyp),'EDOTJavg')

