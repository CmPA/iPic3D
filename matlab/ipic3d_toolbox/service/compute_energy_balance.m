function [Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV, ugradp, pdivv, divUP] = compute_energy_balance( ...
    rho, Jx, Jy, Jz,... 
    Qbulkx, Qbulky, Qbulkz, Qenthx, Qenthy, Qenthz, Qhfx, Qhfy, Qhfz, ...
    Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, x, y, z, dx, dy, dz, qom, radius,cyl)

%divQbulk = compute_div(x,y,z,smooth3(Qbulkx,'box',5),smooth3(Qbulky,'box',5),smooth3(Qbulkz,'box',5));    
divQbulk = compute_div(x,y,z,Qbulkx,Qbulky,Qbulkz, radius, cyl);
divQenth = compute_div(x,y,z,Qenthx,Qenthy,Qenthz, radius, cyl);
divQhf = compute_div(x,y,z,Qhfx,Qhfy,Qhfz, radius, cyl);
Uth = (Pxx + Pyy + Pzz)/2;
Ubulk = ((Jx./rho).^2 + (Jy./rho).^2 + (Jz./rho).^2) .* rho /2/qom;  
    

Vx=imgaussfilt3(Jx./rho,radius);
Vy=imgaussfilt3(Jy./rho,radius);
Vz=imgaussfilt3(Jz./rho,radius);

tmp = compute_div(x,y,z,Pxx,Pxy,Pxz, radius, cyl);
udivP = tmp.* Vx;
tmp = compute_div(x,y,z,Pxy,Pyy,Pyz, radius, cyl);
udivP = udivP + tmp.* Vy;
tmp = compute_div(x,y,z,Pxz,Pyz,Pzz, radius, cyl);
udivP = udivP + tmp.* Vz;

p=(Pxx+Pyy+Pzz)/3;
pdivv = imgaussfilt3(p,radius).*compute_div(x,y,z,Vx,Vy,Vz, radius, cyl);

[tx, ty, tz] = gradient(imgaussfilt3(permute(p,[2 1 3]),radius), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
ugradp = tx.*Vx + ty.*Vy + tz.*Vz;

[tx, ty, tz] = gradient(permute(Vx,[2 1 3]), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = tx.*imgaussfilt3(Pxx,radius)+ ty.*imgaussfilt3(Pxy,radius) +tz.* imgaussfilt3(Pxz,radius);
[tx, ty, tz] = gradient(permute(Vy,[2 1 3]), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*imgaussfilt3(Pxy,radius) + ty.*imgaussfilt3(Pyy,radius) +tz.* imgaussfilt3(Pyz,radius);
[tx, ty, tz] = gradient(permute(Vz,[2 1 3]), dx, dy, dz);
tx=permute(tx,[2 1 3]);ty=permute(ty,[2 1 3]);tz=permute(tz,[2 1 3]);
PgradV = PgradV + tx.*imgaussfilt3(Pxz,radius) + ty.*imgaussfilt3(Pyz,radius) +tz.* imgaussfilt3(Pzz,radius);

divUP = compute_div(x,y,z,Pxx.*Vx,Pxy.*Vy,Pxz.*Vz, radius, cyl);
divUP = divUP + compute_div(x,y,z,Pxy.*Vx,Pyy.*Vy,Pyz.*Vz, radius, cyl);
divUP = divUP + compute_div(x,y,z,Pxz.*Vx,Pyz.*Vy,Pzz.*Vz, radius, cyl);
