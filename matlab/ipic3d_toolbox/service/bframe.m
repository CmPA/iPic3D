function [Apar Aperp1 Aperp2] = bframe(Bx, By, Bz, Ax, Ay, Az)

Bmod= sqrt(Bx.^2+By.^2+Bz.^2);
B2D= sqrt(Bx.^2+By.^2);
perp2x=Bz.*Bx./(Bmod.*B2D);
perp2y=Bz.*By./(Bmod.*B2D);
perp2z=-B2D./Bmod;

Apar= (Ax.*Bx+Ay.*By+Az.*Bz)./(Bmod+1e-10);

Aperp1=(By.*Ax-Bx.*Ay)./B2D;
Aperp2=perp2x.*Ax+perp2y.*Ay+perp2z.*Az;

end