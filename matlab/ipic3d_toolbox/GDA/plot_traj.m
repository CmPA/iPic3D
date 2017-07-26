name=['autotraject' num2str(counter)]
%load([name '.mat'])
Lx=100;
Ly=30;
xgsm=-y(:,1)/Lx*30-15;
zgsm=-y(:,2)/Ly*12-9;
subplot(2,2,1)
%plot(y(:,1),y(:,2),y(1,1),y(1,2),'ro',y(end,1),y(end,2),'kx')
plot(xgsm,zgsm,xgsm(1),zgsm(1),'ro',xgsm(end),zgsm(end),'kx')
set(gca,'xdir','reverse')
axis tight
xlabel('x/RE')
ylabel('z/RE')
subplot(2,2,2)
%plot(y(:,1),y(:,4),y(1,1),y(1,4),'ro',y(end,1),y(end,4),'kx')
plot(xgsm,y(:,4),xgsm(1),y(1,4),'ro',xgsm(end),y(end,4),'kx')
set(gca,'xdir','reverse')
axis tight
xlabel('x/RE')
ylabel('vxgsm/c')
subplot(2,2,3)
plot(y(:,4),y(:,5),y(1,4),y(1,5),'ro',y(end,4),y(end,5),'kx')
xlabel('vx/c')
ylabel('vzgsm/c')
subplot(2,2,4)
plot(y(:,4),y(:,6),y(1,4),y(1,6),'ro',y(end,4),y(end,6),'kx')
xlabel('vxgsm/c')
ylabel('vygsm/c')
print('-dpdf', [name '.pdf']) 
