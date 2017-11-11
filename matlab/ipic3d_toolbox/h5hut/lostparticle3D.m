
function [value,isterminal,direction] = lostparticle3D(t,y)
global Lx Ly Lz Rout
xp=y(1:3);

r = sqrt((xp(3)-Lz/2)^2+(xp(2)-Ly/2)^2+(xp(1)-Lx/2)^2);
value=(r- Rout);
isterminal=1;
direction=0;
end

