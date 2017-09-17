function  udivP = compute_udivP(x,y,z,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz, N)

method='gaussian'
radius=5;
tmp = divergence(x,y,z,smooth3(permute(Pxx,[2 1 3]),method,radius), smooth3(permute(Pxy, [2 1 3]),method,radius), smooth3(permute(Pxz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = tmp.* Jx ./N;
tmp = divergence(x,y,z,smooth3(permute(Pxy,[2 1 3]),method,radius), smooth3(permute(Pyy, [2 1 3]),method,radius), smooth3(permute(Pyz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Jy ./N;
tmp = divergence(x,y,z,smooth3(permute(Pxz,[2 1 3]),method,radius), smooth3(permute(Pyz, [2 1 3]),method,radius), smooth3(permute(Pzz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Jz ./N;

end