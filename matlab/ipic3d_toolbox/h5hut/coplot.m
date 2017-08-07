function [ym] = coplot(it,x,y,f,ay,xlab,ylab,titolo,range1,range2)
global contours

f=double(f);
ay=double(ay);
x=double(x);
y=double(y);


hold off
if(range1(2)==range1(1))
    range1(1)=min(f(:));
    range1(2)=max(f(:));
end

if(range2(2)==range2(1))
    range2(1)=min(ay(:));
    range2(2)=max(ay(:));
end

if(range1(2)==range1(1))
    return
end
fmin=min(range1);
fmax=max(range1);
imagesc([min(x(:)) max(x(:))],[min(y(:)) max(y(:))],f);

title(titolo,'fontsize',[14])

colorbar
set(gca,'fontsize',[14])

hold on



caxis(range1)

%ay2=((ay-range2(1))/(range2(2)-range2(1))*(range1(2)-range1(2))+range1(1));


ff=contour(x,y,ay,linspace(range2(1),range2(2),30),'w');


xlabel(xlab,'fontsize',[18])
ylabel(ylab,'fontsize',[18])
colormap hsv

axis equal
axis tight
set(gcf, 'Renderer', 'zbuffer');
print('-dpng',[titolo '_cycle' num2str(it) '.png'])
print('-depsc','-painters',[titolo '_cycle' num2str(it) '.eps'])
%saveas(gcf,[tit '.fig'])

