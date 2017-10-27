function dx = newton_interp(t, x)
  global ex ey ez bx by bz xg yg  Lx Ly qom Rout
xp=x(1:3);
vp=x(4:6);
dx = zeros(6,1);

% if(xp(1)<0 | xp(1)>Lx) 
%     %disp('out in x')
% return
% end
% if(xp(2)<0 | xp(2)>Ly)
%     %disp('out in y')
% return
% end

r = sqrt(xp(3)^2+xp(1)^2);
theta = atan2(xp(3),xp(1));

% if (sqrt((xp(2)-Ly/2)^2+r^2)>Rout)
%      %disp('out in R')
% return
% end

forward=1;
if(forward)
tdir=1;
else 
tdir=-1;
end

dx(1) = tdir*vp(1);
dx(2) = tdir*vp(2);
dx(3) = tdir*vp(3);

methodflag=2;
Br = qinterp2(xg,yg,bx,r,xp(2),methodflag); % Code x is cylindrical coordiante r
Bp(2) = qinterp2(xg,yg,by,r,xp(2),methodflag); % Code y is cylindrical coordiante z
Btheta = qinterp2(xg,yg,bz,r,xp(2),methodflag); % Code z is cylindrical coordiante theta

Bp(1) = Br*cos(theta) - Btheta * sin(theta);
Bp(3) = Br*sin(theta) + Btheta * cos(theta);


Er = qinterp2(xg,yg,ex,r,xp(2),methodflag); % Code x is cylindrical coordiante r
Ep(2) = qinterp2(xg,yg,ey,r,xp(2),methodflag); % Code y is cylindrical coordiante z
Etheta = qinterp2(xg,yg,ez,r,xp(2),methodflag); % Code z is cylindrical coordiante theta

Ep(1) = Er*cos(theta) - Etheta * sin(theta);
Ep(3) = Er*sin(theta) + Etheta * cos(theta);

Fp=cross(vp,Bp);

dx(4) = tdir*qom*(Ep(1) + Fp(1));
dx(5) = tdir*qom*(Ep(2) + Fp(2));
dx(6) = tdir*qom*(Ep(3) + Fp(3));
