function plot_sat(t,x,var,ysat,filename)
figure(1)
%pcolor(t,x,var) 
pcolor(x,t,var') 
shading interp
colorbar
title([filename '   y=   ' num2str(ysat)])
%xlabel('\omega_{ci}t','fontsize',[14])
%ylabel('x/d_i','fontsize',[14])
ylabel('\omega_{ci}t','fontsize',[14])
xlabel('x/d_i','fontsize',[14])
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng',[filename '.png'])
saveas(gcf,[filename '.fig'])
save([filename '.mat'],'t', 'x', 'var')
close(1)
