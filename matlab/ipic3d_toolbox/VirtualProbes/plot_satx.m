function plot_satx(t,x,var,xsat,filename)
figure(1)
pcolor(t,x,var) 
shading interp
colorbar
title([filename '   x=   ' num2str(xsat)])
xlabel('\omega_{ci}t','fontsize',[14])
ylabel('y/d_i','fontsize',[14])
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng',[filename '.png'])
saveas(gcf,[filename '.fig'])
save([filename '.mat'],'t', 'x', 'var')
close(1)
