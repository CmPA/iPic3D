function [Vexbx, Vexby, Vexbz] = cross_overmod2(Ex, Ey, Ez, Bx, By, Bz,B2)

Vexbx = (Ey.* Bz - Ez.* By)./B2;
Vexby = (Ez.* Bx - Ex.* Bz)./B2;
Vexbz = (Ex.* By - Ey.* Bx)./B2;
