function [integral_phsp]=spaziofasi_int(xx,uu,qq,ndiv,xmin,xmax)

if(nargin<4)
ndiv=100;
end

if(nargin<6)
xmax=max(xx(:));
xmin=min(xx(:));
end

dx=(-xmin+xmax)/ndiv;
xrange=xmin:dx:xmax;

[nbinx,binx]=histc(xx,xrange);
ibin=max(size(nbinx));
integral_phsp = zeros(ibin,1);

for i=1:ibin
    utmp=uu(binx==i);
    qtmp=qq(binx==i);
    integral_phsp(i) = mean(utmp);
end

