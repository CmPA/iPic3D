function [Qenthx,Qenthy,Qenthz,Qbulkx,Qbulky,Qbulkz,Qhfx,Qhfy,Qhfz,Ubulk,Uth]=compute_energy_fluxes(Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,Qx,Qy,Qz,Jx,Jy,Jz,rho, qom)

small=1e-10;
Ubulk = 0.5*(Jx.^2 + Jy.^2 + Jz.^2)./(N+small)/qom;
Qbulkx = Jx./(rho+small).*Ubulk;
Qbulky = Jy./(rho+small).*Ubulk;
Qbulkz = Jz./(rho+small).*Ubulk;

Uth = 0.5*(Pxx+Pyy+Pzz);

Qenthx = Jx./(rho+small).* Uth + (Jx .* Pxx + Jy.* Pxy + Jz.* Pxz) ./(rho+small);
Qenthy = Jy./(rho+small).* Uth + (Jx .* Pxy + Jy.* Pyy + Jz.* Pyz) ./(rho+small);
Qenthz = Jz./(rho+small).* Uth + (Jx .* Pxz + Jy.* Pyz + Jz.* Pzz) ./(rho+small);

Qhfx= Qx - Qbulkx - Qenthx;
Qhfy= Qy - Qbulky - Qenthy;
Qhfz= Qz - Qbulkz - Qenthz;

end