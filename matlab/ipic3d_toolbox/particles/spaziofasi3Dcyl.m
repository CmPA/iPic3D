function [vdf]=spaziofasi3Dcyl(upar,uper,thp,q,u1min, u1max, u2min, u2max, u3min, u3max,ndiv)

dupar=(-u1min+u1max)/ndiv;
duper=(-u2min+u2max)/ndiv;
dth=(-u3min+u3max)/ndiv;
np=max(size(upar));
vdf=zeros(ndiv,ndiv,ndiv);
for ip=1:np

ivpar=floor((upar(ip)-u1min)/dupar);
ivper=floor((uper(ip)-u2min)/duper);
ivth=floor((thp(ip)-u3min)/dth);

if(ivpar>0 & ivpar<ndiv+1 & ivper>0 & ivper<ndiv+1 & ivth>0 & ivth<ndiv+1)

          vdf(ivpar,ivper,ivth)= vdf(ivpar,ivper,ivth)+abs(q(ip));

end
end
          
          %vdf=vdf./sum(vdf(:));

return
