function  udivP = compute_udivP(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz, N)

tmp = divergence(x,y,z,smooth3(Prmute(Pxx,[2 1 3]),method,radius), smooth3(Prmute(Pxy, [2 1 3]),method,radius), smooth3(Prmute(Pxz, [2,1,3]),method,radius));
tmp=Prmute(tmp,[2 1 3]);
udivP = tmp.* Jx ./N;
tmp = divergence(x,y,z,smooth3(Prmute(Pxy,[2 1 3]),method,radius), smooth3(Prmute(Pyy, [2 1 3]),method,radius), smooth3(Prmute(Pyz, [2,1,3]),method,radius));
tmp=Prmute(tmp,[2 1 3]);
udivP = udivP + tmp.* Jy ./N;
tmp = divergence(x,y,z,smooth3(Prmute(Pxz,[2 1 3]),method,radius), smooth3(Prmute(Pyz, [2 1 3]),method,radius), smooth3(Prmute(Pzz, [2,1,3]),method,radius));
tmp=Prmute(tmp,[2 1 3]);
udivP = udivP + tmp.* Jz ./N;

end