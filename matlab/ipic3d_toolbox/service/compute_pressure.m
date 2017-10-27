function [Pxx,Pyy,Pzz,Pxy,Pxz,Pyz]=compute_pressure(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Jx,Jy,Jz,N, qom)

small=1e-10;
Pxx = (Pxx - Jx.*Jx ./ (N+small) ) /qom;
Pyy = (Pyy - Jy.*Jy ./ (N+small) ) /qom;
Pzz = (Pzz - Jz.*Jz ./ (N+small) ) /qom;
Pxy = (Pxy - Jx.*Jy ./ (N+small) ) /qom;
Pxz = (Pxz - Jx.*Jz ./ (N+small) ) /qom;
Pyz = (Pyz - Jy.*Jz ./ (N+small) ) /qom;

end