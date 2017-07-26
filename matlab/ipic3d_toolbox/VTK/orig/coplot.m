function [ym] = coplot(x,y,f,ay,xlab,ylab,tit,lmt,col)
global contours

f=double(f);
ay=double(ay);
x=double(x);
y=double(y);


%ym=max(y(:))/2;
nargin
hold off
if nargin<8
fmin1=min(f(:));
fmax1=max(f(:));
gg=pcolor(x,y,f);
else
fmin1=min(min((f(lmt(1):lmt(2),lmt(3):lmt(4)))));
fmax1=max(max((f(lmt(1):lmt(2),lmt(3):lmt(4)))));
gg=pcolor(x(lmt(1):lmt(2),lmt(3):lmt(4)),y(lmt(1):lmt(2),lmt(3):lmt(4)),f(lmt(1):lmt(2),lmt(3):lmt(4)));
end
title(tit,'fontsize',[18])
shading interp
colorbar
%set(gg,'fontsize',[18])
hold on
if nargin<8
aymin=min(ay(:));
aymax=max(ay(:));
else
aymin=min(min((ay(lmt(1):lmt(2),lmt(3):lmt(4)))));
aymax=max(max((ay(lmt(1):lmt(2),lmt(3):lmt(4)))));
end
fmax=max(fmax1,-fmin1);
fmin=0*min(fmin1,-fmax1);

caxis([fmin fmax])

ay2=((ay-aymin)/(aymax-aymin)*(fmax-fmin)+fmin);

if contours
if nargin<8
ff=contour(x,y,ay2,30,'w');
max(ay2(:))
else
ff=contour(x(lmt(1):lmt(2),lmt(3):lmt(4)),y(lmt(1):lmt(2),lmt(3):lmt(4)),ay2(lmt(1):lmt(2),lmt(3):lmt(4)),30,'w');
end
end

if nargin==9
caxis(col)
end

xlabel(xlab,'fontsize',[18])
ylabel(ylab,'fontsize',[18])
%set(ff,'fontsize',[18])
set(gca,'fontsize',[18])
%load cm_new
%colormap(cm_kbwrk)
axis equal
axis tight
set(gcf, 'Renderer', 'zbuffer');
print('-dpng',[tit '.png'])
%saveas(gcf,[tit '.fig'])
end
