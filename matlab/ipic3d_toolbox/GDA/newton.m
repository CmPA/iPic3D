function dx = newton(t, x)
global Ex Ey Ez Bx By Bz xc yc qom Lx Ly
xp=x(1:3);
vp=x(4:6);
dx = zeros(6,1);

if(xp(1)<0 | xp(1)>Lx) 
return
end
if(xp(2)<0 | xp(2)>Ly)
return
end

forward=0;
if(forward)
tdir=1;
else 
tdir=-1;
end

dx(1) = tdir*vp(1);
dx(2) = tdir*vp(2);
dx(3) = tdir*vp(3);
Bp(1) = interp2(xc,yc,Bx',xp(1),xp(2));
Bp(2) = interp2(xc,yc,By',xp(1),xp(2));
Bp(3) = interp2(xc,yc,Bz',xp(1),xp(2));
Ep(1) = interp2(xc,yc,Ex',xp(1),xp(2));
Ep(2) = interp2(xc,yc,Ey',xp(1),xp(2));
Ep(3) = interp2(xc,yc,Ez',xp(1),xp(2));
Fp=cross(vp,Bp);

dx(4) = tdir*qom*(Ep(1) + Fp(1));
dx(5) = tdir*qom*(Ep(2) + Fp(2));
dx(6) = tdir*qom*(Ep(3) + Fp(3));
