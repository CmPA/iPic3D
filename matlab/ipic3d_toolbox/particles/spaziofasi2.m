function [totnum,nbinu,xrange,urange]=spaziofasi2(xx,uu,qq,cyl,xmin,xmax,umin,umax,ndiv)

if(nargin<9)
ndiv=60;
end


dx=(-xmin+xmax)/ndiv;
du=(-umin+umax)/ndiv;
xrange=xmin:dx:xmax;
urange=umin:du:umax;
[nbinx,binx]=histc(xx,xrange);
ibin=max(size(nbinx));
nbinu=[];
totnum=0;
for i=1:ibin
    utmp=uu(binx==i);
    qtmp=qq(binx==i);
    [nbin,bin]=histc(utmp,urange);
    totnum=totnum+sum(nbin);
    nbib=nbin*0;
    for j=1:ibin;nbin(j)=sum(abs(qtmp(bin==j)));end
    nbin_tmp=nbin.*(xrange(i)+dx/2).^cyl;
    %size(nbin_tmp)
    nbin_tmp=reshape(nbin_tmp,ndiv+1,1);
    nbinu=[nbinu nbin_tmp];
end

