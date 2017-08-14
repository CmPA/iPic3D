function [Qenthx,Qenthy,Qenthz,Qbulkx,Qbulky,Qbulkz,Qhfx,Qhfy,Qhfz]=compute_energy_fluxes(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Qx,Qy,Qz,Jx,Jy,Jz,N, qom)

Ubulk = 0.5*(Jx.^2 + Jy.^2 + Jz.^2)./N/qom;
Qbulkx = Jx./N.*Ubulk;
Qbulky = Jy./N.*Ubulk;
Qbulkz = Jz./N.*Ubulk;

Uth = 0.5*(Pxx+Pyy+Pzz);

Qenthx = Jx./N.* Uth + (Jx .* Pxx + Jy.* Pxy + Jz.* Pxz) ./N;
Qenthy = Jy./N.* Uth + (Jx .* Pxy + Jy.* Pyy + Jz.* Pyz) ./N;
Qenthz = Jz./N.* Uth + (Jx .* Pxz + Jy.* Pyz + Jz.* Pzz) ./N;

Qhfx= Qx - Qbulkx - Qenthx;
Qhfy= Qy - Qbulky - Qenthy;
Qhfz= Qz - Qbulkz - Qenthz;
end