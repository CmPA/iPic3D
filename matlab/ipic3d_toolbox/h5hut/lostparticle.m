
function [value,isterminal,direction] = lostparticle(t,y)
global ex ey ez bx by bz xg yg  Lx Ly qom Rout
xp=y(1:3);

r = sqrt(xp(3)^2+xp(1)^2);
value=(sqrt(0*(xp(2)-Ly/2)^2+r^2)- Rout);
isterminal=1;
direction=0;
end

