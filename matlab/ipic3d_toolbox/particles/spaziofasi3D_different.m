function [vdf,umin,vmin,wmin,du,dv,dw]=spaziofasi3D_different(u,v,w,q,ndiv)

umin=min(u(:));umax=max(u(:));
vmin=min(v(:));vmax=max(v(:));
wmin=min(w(:));wmax=max(w(:));
du=(-umin+umax)/ndiv;
dv=(-vmin+vmax)/ndiv;
dw=(-wmin+wmax)/ndiv;
np=max(size(u));
vdf=zeros(ndiv,ndiv,ndiv);
for ip=1:np

ivx=floor((u(ip)-umin)/du);
ivy=floor((v(ip)-vmin)/dv);
ivz=floor((w(ip)-wmin)/dw);

if(ivx>0 & ivx<ndiv+1 & ivy>0 & ivy<ndiv+1 & ivz>0 & ivz<ndiv+1)

          vdf(ivx,ivy,ivz)= vdf(ivx,ivy,ivz)+abs(q(ip));

end
end
          
          %vdf=vdf./sum(vdf(:));

return
