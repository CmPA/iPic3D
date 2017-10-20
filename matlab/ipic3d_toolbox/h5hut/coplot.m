function [ym] = coplot(it,x,y,f,ay,xlab,ylab,titolo,range1,range2,figure_number)
global contours dx dy color_choice
if(nargin>10) 
    figure(figure_number) 
end;
    
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

load cm_new
load cm_multi4
        color = 'w'
        if(color_choice==0)
        colormap jet
        elseif (color_choice==1)
        colormap(cm_kbwrk)
        elseif (color_choice==2)
        colormap cool
        elseif (color_choice==3)
        colormap(cm_cool_hot_2); color = 'k'
        elseif (color_choice==4)
        colormap(flipud(hot))
        end
        

ff=contour(x,y,ay,linspace(range2(1),range2(2),60),color);


xlabel(xlab,'fontsize',[18])
ylabel(ylab,'fontsize',[18])
%colormap hsv % hsv prefers black contour lines, not white

axis equal
axis tight
set(gcf, 'Renderer', 'zbuffer');
print('-dpng',[titolo '_cycle' num2str(it) '.png'])
print('-depsc','-painters',[titolo '_cycle' num2str(it) '.eps'])
%saveas(gcf,[tit '.fig'])

