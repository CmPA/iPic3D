cyl = 0

[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);

% Compute J dot E
JedotE=dot(Jex,Jey,Jez,Ex,Ey,Ez);


JedotEsm=dot(imgaussfilt3(Jex,radius),imgaussfilt3(Jey,radius),imgaussfilt3(Jez,radius), ...
    imgaussfilt3(Ex,radius),imgaussfilt3(Ey,radius),imgaussfilt3(Ez,radius));


JidotE=dot(Jix,Jiy,Jiz,Ex,Ey,Ez);

JdotE=JedotE+JidotE;

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx, By, Bz);
divS = compute_div(x,y,z,Sx,Sy,Sz,radius, cyl);


%
% Electrons
%


Vex=Jex./rhoe;Vey=Jey./rhoe;Vez=Jez./rhoe;
divVe = compute_div(x,y,z,Vex,Vey,Vez, radius, cyl);
Vix=Jix./rhoi;Viy=Jiy./rhoi;Viz=Jiz./rhoi;
divVi = compute_div(x,y,z,Vix,Viy,Viz, radius, cyl);

Vex=-imgaussfilt(squeeze(mean(Vex(:,jr,:),2)),radius);
Vey=imgaussfilt(squeeze(mean(Vey(:,jr,:),2)),radius);
Vez=imgaussfilt(squeeze(mean(Vez(:,jr,:),2)),radius);
AAze=vecpot(xc,zc,Vex,Vez);

Vix=-imgaussfilt(squeeze(mean(Vix(:,jr,:),2)),radius);
Viy=imgaussfilt(squeeze(mean(Viy(:,jr,:),2)),radius);
Viz=imgaussfilt(squeeze(mean(Viz(:,jr,:),2)),radius);
AAzi=vecpot(xc,zc,Vix,Viz);

AAz=vecpot(xc,yc,-mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));
Vx=-mean(Bx(:,:,kr),3);
Vy=mean(By(:,:,kr),3);

if(electrons)


[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV, ugradp, pdivv, divUP] = compute_energy_balance( ...
    rhoe, Jex, Jey, Jez,... 
    Qbulkex, Qbulkey, Qbulkez, Qenthex, Qenthey, Qenthez, Qhfex, Qhfey, Qhfez, ...
    Pexx, Peyy, Pezz, Pexy, Pexz, Peyz, x, y, z, dx, dy, dz, qom_ele, radius,cyl);
    
    
DUbulkDt = JedotE - Ubulk.*divVe - udivP;
DUthDt = -Uth.*divVe -PgradV;

end    


if(ions)
[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV, ugradp, pdivv, divUP] = compute_energy_balance( ...
    rhoi, Jix, Jiy, Jiz,... 
    Qbulkix, Qbulkiy, Qbulkiz, Qenthix, Qenthiy, Qenthiz, Qhfix, Qhfiy, Qhfiz, ...
    Pixx, Piyy, Pizz, Pixy, Pixz, Piyz, x, y, z, dx, dy, dz, 1.0, radius,cyl);

DUbulkDt = JidotE - Ubulk.*divVi - udivP;
DUthDt = -Uth.*divVi -PgradV;

end
