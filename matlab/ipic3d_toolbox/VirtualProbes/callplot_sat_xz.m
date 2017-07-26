
plot_satx((1:n1)*Dt*B0x,yyp,EEX,xp(isatx),'EX')
clear EEX

plot_satx((1:n1)*Dt*B0x,yyp,EEY,xp(isatx),'EY')
clear EEY

plot_satx((1:n1)*Dt*B0x,yyp,EEZ,xp(isatx),'EZ')
clear EEZ

plot_satx((1:n1)*Dt*B0x,yyp,BBX,xp(isatx),'BX')

plot_satx((1:n1)*Dt*B0x,yyp,BBY,xp(isatx),'BY')

plot_satx((1:n1)*Dt*B0x,yyp,BBZ,xp(isatx),'BZ')


Navg=100
dBBX=diff(BBX,1,2)/Dt;
dBBXavg=tsmovavg(dBBX(:,1:npmax-1),'s',Navg);
plot_satx((1.5:n1-.5)*Dt*B0x,yyp,dBBX,xp(isatx),'dBXdt')
clear dBBX


dBBY=diff(BBY,1,2)/Dt;
dBBYavg=tsmovavg(dBBY(:,1:npmax-1),'s',Navg);
plot_satx((1.5:n1-.5)*Dt*B0x,yyp,dBBY,xp(isatx),'dBYdt')
clear dBBY

dBBZ=diff(BBZ,1,2)/Dt;
dBBZavg=tsmovavg(dBBZ(:,1:npmax-1),'s',Navg);
plot_satx((1.5:n1-.5)*Dt*B0x,yyp,dBBZ,xp(isatx),'dBZdt')
clear dBBZ

plot_satx((1:n1)*Dt*B0x,yyp,NE,xp(isatx),'NE')

%plot_satx((1:n1)*Dt*B0x,yyp,JX,xp(isatx),'UEX')
%clear JX
%plot_satx((1:n1)*Dt*B0x,yyp,JY,xp(isatx),'UEY')
%clear JY
%plot_satx((1:n1)*Dt*B0x,yyp,JZ,xp(isatx),'UEZ')
%clear JZ

plot_satx((1:n1)*Dt*B0x,yyp,sqrt((BBX.^2+BBY.^2+BBZ.^2)./abs(NE))./B0x,xp(isatx),'VA')

aaa=angle(BBX+i*BBZ);
plot_satx((1:n1)*Dt*B0x,yyp,aaa,xp(isatx),'AngleBZoBX')
return
clear BBZ
aaa=angle(BBX+i*BBY);
plot_satx((1:n1)*Dt*B0x,yyp,aaa,xp(isatx),'AngleBYoBX')
clear BBX
clear BBY

plot_satx((1:n1)*Dt*B0x,yyp,SSX,xp(isatx),'SX')
clear SSX
plot_satx((1:n1)*Dt*B0x,yyp,SSY,xp(isatx),'SY')
clear SSY
plot_satx((1:n1)*Dt*B0x,yyp,SSZ,xp(isatx),'SZ')
clear SSZ

plot_satx((Navg+1:n1-1)*Dt*B0x,yyp,dBBXavg(:,Navg+1:n1-1),xp(isatx),'dBXavg')
clear dBBXavg
plot_satx((Navg+1:n1-1)*Dt*B0x,yyp,dBBYavg(:,Navg+1:n1-1),xp(isatx),'dBYavg')
clear dBBYavg
plot_satx((Navg+1:n1-1)*Dt*B0x,yyp,dBBZavg(:,Navg+1:n1-1),xp(isatx),'dBZavg')
clear dBBZavg



if(detrended)
plot_satx((Ndetrend+1:n1)*Dt*B0x,yyp,SSXdet(:,Ndetrend+1:n1),xp(isatx),'SXdet')
clear SSXdet
plot_satx((2*Ndetrend+1:n1)*Dt*B0x,yyp,SSXdetavg(:,2*Ndetrend+1:n1),xp(isatx),'SXdetavg')
clear SSXdetavg
end

%plot_satx((1:n1)*Dt*B0x,yyp,EDOTJ,xp(isatx),'EDOTJ')


%Navg=100;
%EDOTJ=tsmovavg(EDOTJ(:,1:npmax),'s',Navg);
%EDOTJ=(2*EDOTJ(2:end-1,:)+EDOTJ(1:end-2,:)+EDOTJ(3:end,:))/4;
%plot_satx((2*Navg+1:npmax)*Dt*B0x,yyp,EDOTJ(:,2*Navg+1:npmax),xp(isatx),'EDOTJavg')

