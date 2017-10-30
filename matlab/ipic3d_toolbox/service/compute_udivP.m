function  udivP = compute_udivP(x,y,z,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz, N)

small=1e-10;
method='gaussian'
radius=1;
%tmp = compute_div(x,y,z,smooth3(Pxx,method,radius),smooth3(Pxy,method,radius),smooth3(Pxz,method,radius));
tmp = compute_div(x,y,z,Pxx,Pxy,Pxz);
udivP = tmp.* Jx ./(N+small);
%tmp = compute_div(x,y,z,smooth3(Pxy,method,radius),smooth3(Pyy,method,radius),smooth3(Pyz,method,radius));
tmp = compute_div(x,y,z,Pxy,Pyy,Pyz);
udivP = udivP + tmp.* Jy ./(N+small);
%tmp = compute_div(x,y,z,smooth3(Pxz,method,radius),smooth3(Pyz,method,radius),smooth3(Pzz,method,radius));
tmp = compute_div(x,y,z,Pxz,Pyz,Pzz);
udivP = udivP + tmp.* Jz ./(N+small);

end