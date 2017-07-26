function vol3dtool(fig_target)
% Tool to edit color/alpha channels when 
% working with vol3d.m function.
%
% vol3dtool(FIG) Pass in figure handle
% vol3dtool  Use GCF.
%
% See example in vol3d.m

% Copyright Joe Conti, 2004

if nargin==0
    fig_target = gcf;
end

fig = hgload('vol3dtool.fig');
set(fig,'units','pixels');
set(fig,'Name','Color/Alpha Channel Edit');
red_ax = findall(fig,'type','axes','tag','red');
set(red_ax,'units','pixels');
green_ax = findall(fig,'type','axes','tag','green');
blue_ax = findall(fig,'type','axes','tag','blue');
alpha_ax = findall(fig,'type','axes','tag','alpha');

cmap = get(fig_target,'Colormap');
amap = get(fig_target,'Alphamap');

l_info.xdata = 1:length(cmap(:,1));
l_info.ydata = cmap(:,1);
l_info.parent = red_ax;
l_info.markersize = 5;
l_info.linewidth = 3;
l_info.color = 'r';
l_info.markerfacecolor = 'r';
l_info.buttondown = {@local_line,red_ax,fig,fig_target};
ylim(red_ax,[-.1,1.1]); 
xlim(red_ax,[1,length(cmap(:,1))]);
h = line(l_info);

l_info.ydata = cmap(:,2);
l_info.parent = green_ax;
l_info.color = 'g';
l_info.markerfacecolor = 'g';
l_info.buttondown = {@local_line,green_ax,fig,fig_target};
ylim(green_ax,[-.1,1.1]); 
xlim(green_ax,[1,length(cmap(:,1))]);
h = line(l_info);

l_info.ydata = cmap(:,3);
l_info.parent = blue_ax;
l_info.color = 'b';
l_info.markerfacecolor = 'b';
l_info.buttondown = {@local_line,blue_ax,fig,fig_target};
ylim(blue_ax,[-.1,1.1]); 
xlim(blue_ax,[1,length(cmap(:,1))]);
h = line(l_info);

l_info.xdata = 1:length(amap);
l_info.ydata = amap;
l_info.parent = alpha_ax;
l_info.color = 'k';
l_info.markerfacecolor = 'k';
l_info.buttondown = {@local_line,alpha_ax,fig,fig_target};
ylim(alpha_ax,[-.1,1.1]); 
xlim(alpha_ax,[1,length(amap)]);
h = line(l_info);

%--------------------------------------------------------%
function local_line(h,evd,ax,fig,fig_target)

cp = get(ax,'CurrentPoint');

% Find closest vertex
xdata = get(h,'xdata');
ydata = get(h,'ydata');
xp = cp(1,1);
yp = cp(2,2);
dist = (xp-xdata).^2 + (yp-ydata).^2;
[val,ind] = min(dist);

if ~isempty(ind)
  set(fig,'WindowButtonMotionFcn',{@local_mouse_motion,h,ax,fig,fig_target,ind});
end

%--------------------------------------------------------%
function local_mouse_motion(obj,evd,h,ax,fig,fig_target,ind)

set(fig,'WindowButtonUpFcn',{@local_mouse_up,fig});
cp = get(ax,'CurrentPoint');
yp = cp(2,2);
if yp<=0 
    yp = 0;
elseif yp >= 1
    yp = 1;
end

ydata = get(h,'ydata');
ydata(ind) = yp;

if ind > 0
  ydata(ind-1) = yp;
end

if ind < length(ydata)
   ydata(ind+1) = yp;
end

set(h,'ydata',ydata);

cmap = get(fig_target,'Colormap');
amap = get(fig_target,'Alphamap');

switch(get(ax,'tag'))
    case 'red'
        cmap(:,1) = ydata';
    case 'green'
        cmap(:,2) = ydata'; 
    case 'blue'
        cmap(:,3) = ydata';
    case 'alpha'
        amap = ydata;
end


set(fig_target,'Colormap',cmap);
set(fig_target,'AlphaMap',amap);


%--------------------------------------------------------%
function local_mouse_up(obj,evd,fig)

set(fig,'WindowButtonMotionFcn','');




