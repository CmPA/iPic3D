% Function to convert an ensemble of particles to a 3D MATLAB velocity ditribution
% matrix
% By Giovanni Lapenta

function [vdf]=spaziofasi3D(u,v,w,q,umin,umax,ndiv)

du=(-umin+umax)/ndiv;
np=max(size(u));
vdf=zeros(ndiv,ndiv,ndiv);
for ip=1:np

ivx=floor((u(ip)-umin)/du);
ivy=floor((v(ip)-umin)/du);
ivz=floor((w(ip)-umin)/du);

if(ivx>0 & ivx<ndiv+1 & ivy>0 & ivy<ndiv+1 & ivz>0 & ivz<ndiv+1)

          vdf(ivx,ivy,ivz)= vdf(ivx,ivy,ivz)+abs(q(ip));

end
end
          
          %vdf=vdf./sum(vdf(:));

return
