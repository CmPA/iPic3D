function [Pxx,Pyy,Pzz,Pxy,Pxz,Pyz, PPAR, PPER1, PPER2]=compute_pressure(Bx, By, Bz, Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz,N, qom)

small=1e-10;
Pxx = (Pxx - Jx.*Jx ./ (N+small) ) /qom;
Pyy = (Pyy - Jy.*Jy ./ (N+small) ) /qom;
Pzz = (Pzz - Jz.*Jz ./ (N+small) ) /qom;
Pxy = (Pxy - Jx.*Jy ./ (N+small) ) /qom;
Pxz = (Pxz - Jx.*Jz ./ (N+small) ) /qom;
Pyz = (Pyz - Jy.*Jz ./ (N+small) ) /qom;

b2D = 1e-10 + Bx.*Bx + By.*By;
b = b2D + Bz.*Bz;
PerP2x = Bz.*Bx./sqrt(b.*b2D);
PerP2y = Bz.*By./sqrt(b.*b2D);
PerP2z = -sqrt(b2D./b);
				
PPAR = Bx.*Pxx.*Bx + 2*Bx.*Pxy.*By + 2*Bx.*Pxz.*Bz;
PPAR = PPAR + By.*Pyy.*By + 2*By.*Pyz.*Bz;
PPAR = PPAR + Bz.*Pzz.*Bz;
				
PPAR = PPAR./b;
				
PPER1 = By.*Pxx.*By - 2*By.*Pxy.*Bx + Bx.*Pyy.*Bx;
				
PPER1 = PPER1./b2D;
				
PPER2 = PerP2x.*Pxx.*PerP2x + 2*PerP2x.*Pxy.*PerP2y + 2*PerP2x.*Pxz.*PerP2z;
PPER2 = PPER2 + PerP2y.*Pyy.*PerP2y + 2*PerP2y.*Pyz.*PerP2z;
PPER2 = PPER2 + PerP2z.*Pzz.*PerP2z;  
				

end