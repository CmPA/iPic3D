function [Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV] = compute_energy_balance( ...
    rho, Jx, Jy, Jz,... 
    Qbulkx, Qbulky, Qbulkz, Qenthx, Qenthy, Qenthz, Qhfx, Qhfy, Qhfz, ...
    Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, x, y, z, dx, dy, dz, qom)

%divQbulk = compute_div(x,y,z,smooth3(Qbulkx,'box',5),smooth3(Qbulky,'box',5),smooth3(Qbulkz,'box',5));    
divQbulk = compute_div(x,y,z,Qbulkx,Qbulky,Qbulkz);
divQenth = compute_div(x,y,z,Qenthx,Qenthy,Qenthz);
divQhf = compute_div(x,y,z,Qhfx,Qhfy,Qhfz);
Uth = (Pxx + Pyy + Pzz)/2;
Ubulk = ((Jx./rho).^2 + (Jy./rho).^2 + (Jz./rho).^2) .* rho /2/qom;  
    
   
radius=5
method='gaussian';

Vx=smooth3(Jx./rho,method,radius);
Vy=smooth3(Jy./rho,method,radius);
Vz=smooth3(Jz./rho,method,radius);

tmp = divergence(x,y,z,smooth3(permute(Pxx,[2 1 3]),method,radius), smooth3(permute(Pxy, [2 1 3]),method,radius), smooth3(permute(Pxz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = tmp.* Vx;
tmp = divergence(x,y,z,smooth3(permute(Pxy,[2 1 3]),method,radius), smooth3(permute(Pyy, [2 1 3]),method,radius), smooth3(permute(Pyz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vy;
tmp = divergence(x,y,z,smooth3(permute(Pxz,[2 1 3]),method,radius), smooth3(permute(Pyz, [2 1 3]),method,radius), smooth3(permute(Pzz, [2,1,3]),method,radius));
tmp=permute(tmp,[2 1 3]);
udivP = udivP + tmp.* Vz;

[tx, ty, tz] = gradient(smooth3(permute(Jx./rho,[2 1 3]),method,radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = tx.*Pxx + ty.*Pxy +tz.* Pxz;
[tx, ty, tz] = gradient(smooth3(permute(Jy./rho,[2 1 3]),method,radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*Pxy + ty.*Pyy +tz.* Pyz;
[tx, ty, tz] = gradient(smooth3(permute(Jz./rho,[2 1 3]),method,radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*Pxz + ty.*Pyz +tz.* Pzz;