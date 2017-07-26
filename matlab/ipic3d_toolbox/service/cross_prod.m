function  [Zx, Zy, Zz] = cross_prod(Vx, Vy, Vz, Bx, By, Bz);


Zx = (Vy.* Bz - Vz.* By);
Zy = (Vz.* Bx - Vx.* Bz);
Zz = (Vx.* By - Vy.* Bx);
