bb=sqrt(bx.^2+by.^2+bz.^2);
epar=(ex.*bx+ey.*by+ez.*bz)./(bb);
upar=((jsx0+jsxb).*bx+(jsy0+jsyb).*by+(jsz0+jszb).*bz)./((rho0+rhob).*bb);
eperpx=ex-epar.*bx./bb;
eperpy=ey-epar.*by./bb;
eperpz=ez-epar.*bz./bb;
eperp=sqrt(eperpx.^2+eperpy.^2+eperpz.^2);
uperpx=(jsx0+jsxb)./(rho0+rhob)-upar.*bx./bb;
uperpy=(jsy0+jsyb)./(rho0+rhob)-upar.*by./bb;
uperpz=(jsz0+jszb)./(rho0+rhob)-upar.*bz./bb;
uperp=sqrt(uperpx.^2+uperpy.^2+uperpz.^2);

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

coplot(axis1,axis2,jsx0,-ay,label1,label2,['J_{x0}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jsx0' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'jsx0' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,jsy0,-ay,label1,label2,['J_{y0}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jsy0' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'jsy0' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,jsz0,-ay,label1,label2,['J_{z0}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jsz0' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'jsz0' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,jsxb,-ay,label1,label2,['J_{xb}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jsxb' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'jsxb' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,jsyb,-ay,label1,label2,['J_{yb}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jsyb' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'jsyb' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,jszb,-ay,label1,label2,['J_{zb}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'jszb' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'jszb' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,ex,-ay,label1,label2,['E_x(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'ex' num2str(ig,'%8.8i')])
if(figs)  saveas(gcf,[film 'ex' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,ey,-ay,label1,label2,['E_y(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'ey' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'ey' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,ez,-ay,label1,label2,['E_z(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'ez' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'ez' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,eperpx,-ay,label1,label2,['E_{\perp,x}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'eperpx' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'eperpx' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,eperpy,-ay,label1,label2,['E_{\perp,y}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'eperpy' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'eperpy' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,eperpz,-ay,label1,label2,['E_{\perp,z}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'eperpz' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'epeprz' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,eperp,ay,label1,label2,['E_{\perp}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'eperp' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'eperp' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,epar,ay,label1,label2,['E_{||}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'epar' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'epar' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,(jsx0+jsxb)./(rho0+rhob),-ay,label1,label2,['u_{x}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'ux' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'ux' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,(jsy0+jsyb)./(rho0+rhob),-ay,label1,label2,['u_{y}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uy' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uy' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,(jsz0+jszb)./(rho0+rhob),-ay,label1,label2,['u_{z}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uz' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uz' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,uperpx,-ay,label1,label2,['u_{\perp,x}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uperpx' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uperpx' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,uperpy,-ay,label1,label2,['u_{\perp,y}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uperpy' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uperpy' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,uperpz,-ay,label1,label2,['u_{\perp,z}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uperpz' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uperpz' num2str(ig,'%8.8i') '.fig']); end
close all

coplot(axis1,axis2,uperp,ay,label1,label2,['u_{\perp}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'uperp' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'uperp' num2str(ig,'%8.8i') '.fig']); end
close all


coplot(axis1,axis2,upar,ay,label1,label2,['u_{||}(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'upar' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'upar' num2str(ig,'%8.8i') '.fig']); end
close all

if(background)
coplot(axis1,axis2,rhob*4*pi,ay,label1,label2,['\rho_b(\omega_{ci}t=' num2str(time) ')'])
set(gcf,'Renderer','zbuffer');
print('-dpng',[film 'rhob' num2str(ig,'%8.8i')])
if(figs) saveas(gcf,[film 'rhob' num2str(ig,'%8.8i') '.fig']); end
close all
end

