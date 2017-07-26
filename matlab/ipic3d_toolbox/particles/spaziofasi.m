function [totnum,nbinu,xrange,urange]=spaziofasi(xx,uu,qq,cyl,umin,umax,ndiv)

if(nargin<7)
ndivv=100;
ndivx=50;
end

xmin=min(xx(:))
xmax=max(xx(:))
if(nargin<5)
umin=min(uu(:))
umax=max(uu(:))
end
% umin=-.4/4;
% umax=.4/4;
dx=(-xmin+xmax)/ndivx;
du=(-umin+umax)/ndivv;
xrange=xmin:dx:xmax;
urange=umin:du:umax;
[nbinx,binx]=histc(xx,xrange);
%ibin=max(size(nbinx))
nbinu=[];
totnum=0;
for i=1:ndivx+1
    utmp=uu(binx==i);
    qtmp=qq(binx==i);
    [nbin,bin]=histc(utmp,urange);
    totnum=totnum+sum(nbin);
    nbib=nbin*0;
    for j=1:ndivv+1;nbin(j)=sum(abs(qtmp(bin==j)));end
    nbinu=[nbinu;nbin.*(xrange(i)+dx/2).^cyl];
end

