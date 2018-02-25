function [Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV] = compute_energy_balance( ...
    rho, Jx, Jy, Jz,... 
    Qbulkx, Qbulky, Qbulkz, Qenthx, Qenthy, Qenthz, Qhfx, Qhfy, Qhfz, ...
    Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, x, y, z, dx, dy, dz, qom, radius)

%divQbulk = compute_div(x,y,z,smooth3(Qbulkx,'box',5),smooth3(Qbulky,'box',5),smooth3(Qbulkz,'box',5));    
divQbulk = compute_div(x,y,z,Qbulkx,Qbulky,Qbulkz, radius, 1);
divQenth = compute_div(x,y,z,Qenthx,Qenthy,Qenthz, radius, 1);
divQhf = compute_div(x,y,z,Qhfx,Qhfy,Qhfz, radius, 1);
Uth = (Pxx + Pyy + Pzz)/2;
Ubulk = ((Jx./rho).^2 + (Jy./rho).^2 + (Jz./rho).^2) .* rho /2/qom;  
    

Vx=imgaussfilt3(Jx./rho,radius);
Vy=imgaussfilt3(Jy./rho,radius);
Vz=imgaussfilt3(Jz./rho,radius);

tmp = compute_div(x,y,z,Pxx,Pxy,Pxz, radius, 1);
udivP = tmp.* Vx;
tmp = compute_div(x,y,z,Pxy,Pyy,Pyz, radius, 1);
udivP = udivP + tmp.* Vy;
tmp = compute_div(x,y,z,Pxz,Pyz,Pzz, radius, 1);
udivP = udivP + tmp.* Vz;

[tx, ty, tz] = gradient(imgaussfilt3(permute(Jx./rho,[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = tx.*Pxx + ty.*Pxy +tz.* Pxz;
[tx, ty, tz] = gradient(imgaussfilt3(permute(Jy./rho,[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*Pxy + ty.*Pyy +tz.* Pyz;
[tx, ty, tz] = gradient(imgaussfilt3(permute(Jz./rho,[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*Pxz + ty.*Pyz +tz.* Pzz;