function dx = newton_interp3D(t, x)
  global ex ey ez bx by bz xn yn zn  Lx Ly Lz qom Rout
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

%r = sqrt(xp(3)^2+xp(1)^2);
%theta = atan2(xp(3),xp(1));

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

Bp(1) = interpn(xn,yn,zn,bx,xp(1),xp(2),xp(3));
Bp(2) = interpn(xn,yn,zn,by,xp(1),xp(2),xp(3));
Bp(3) = interpn(xn,yn,zn,bz,xp(1),xp(2),xp(3));

Ep(1) = interpn(xn,yn,zn,ex,xp(1),xp(2),xp(3));%*time_form(t);
Ep(2) = interpn(xn,yn,zn,ey,xp(1),xp(2),xp(3));%*time_form(t);
Ep(3) = interpn(xn,yn,zn,ez,xp(1),xp(2),xp(3));%*time_form(t);

% Bp(1) = lininterp3(xn,yn,zn,bx,xp(1),xp(2),xp(3));
% Bp(2) = lininterp3(xn,yn,zn,by,xp(1),xp(2),xp(3));
% Bp(3) = lininterp3(xn,yn,zn,bz,xp(1),xp(2),xp(3));
% 
% Ep(1) = lininterp3(xn,yn,zn,ex,xp(1),xp(2),xp(3));
% Ep(2) = lininterp3(xn,yn,zn,ey,xp(1),xp(2),xp(3));
% Ep(3) = lininterp3(xn,yn,zn,ez,xp(1),xp(2),xp(3));


Fp=cross(vp,Bp);

dx(4) = tdir*qom*(Ep(1) + Fp(1));
dx(5) = tdir*qom*(Ep(2) + Fp(2));
dx(6) = tdir*qom*(Ep(3) + Fp(3));
