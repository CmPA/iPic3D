function [Pxx,Pyy,Pzz,Pxy,Pxz,Pyz]=compute_pressure(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz,N, qom)

Pxx = (Pxx - Jx.*Jx ./ N ) /qom;
Pyy = (Pxx - Jy.*Jy ./ N ) /qom;
Pzz = (Pzz - Jz.*Jz ./ N ) /qom;
Pxy = (Pxy - Jx.*Jy ./ N ) /qom;
Pxz = (Pxz - Jx.*Jz ./ N ) /qom;
Pyz = (Pyz - Jy.*Jz ./ N ) /qom;

end