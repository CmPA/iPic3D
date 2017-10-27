function [Qenthx,Qenthy,Qenthz,Qbulkx,Qbulky,Qbulkz,Qhfx,Qhfy,Qhfz,Ubulk,Uth]=compute_energy_fluxes(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Qx,Qy,Qz,Jx,Jy,Jz,N, qom)

small=1e-10;
Ubulk = 0.5*(Jx.^2 + Jy.^2 + Jz.^2)./(N+small)/qom;
Qbulkx = Jx./(N+small).*Ubulk;
Qbulky = Jy./(N+small).*Ubulk;
Qbulkz = Jz./(N+small).*Ubulk;

Uth = 0.5*(Pxx+Pyy+Pzz);

Qenthx = Jx./(N+small).* Uth + (Jx .* Pxx + Jy.* Pxy + Jz.* Pxz) ./(N+small);
Qenthy = Jy./(N+small).* Uth + (Jx .* Pxy + Jy.* Pyy + Jz.* Pyz) ./(N+small);
Qenthz = Jz./(N+small).* Uth + (Jx .* Pxz + Jy.* Pyz + Jz.* Pzz) ./(N+small);

Qhfx= Qx - Qbulkx - Qenthx;
Qhfy= Qy - Qbulky - Qenthy;
Qhfz= Qz - Qbulkz - Qenthz;

end