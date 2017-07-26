
time=ig*wci*Dt

figs

coplot(axis1,axis2,bx,-ay,label1,label2,['B_x(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'bx' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'bx' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,by,-ay,label1,label2,['B_y(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'by' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'by' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,bz,-ay,label1,label2,['B_z(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'bz' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'bz' num2str(ig,'%8.8i') '.fig']); end
close all


