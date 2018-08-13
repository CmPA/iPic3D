cyl = 0;

[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);


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



AAz=vecpot(xc,yc,-mean(Bx(:,:,kr),3),mean(By(:,:,kr),3));
Vx=signx*mean(Bx(:,:,kr),3);
Vy=mean(By(:,:,kr),3);

if(electrons)
    
Vmhdx=(Vix+Vex/abs(qom))/(1+1/abs(qom));
Vmhdy=(Viy+Vey/abs(qom))/(1+1/abs(qom));
Vmhdz=(Viz+Vez/abs(qom))/(1+1/abs(qom));

Epx = Ex + (Vmhdy.*Bz - Vmhdz.*By);
Epy = Ey + (Vmhdz.*Bx - Vmhdx.*Bz);
Epz = Ez + (Vmhdx.*By - Vmhdy.*Bx);

JdotEp=(Jex+Jix).*Epx + (Jey+Jiy).*Epy + (Jez+Jiz).*Epz;

[Uth, Ubulk, divQbulk, divQenth, divQhf,  udivP, PgradV, ugradp, pdivv, divUP] = compute_energy_balance( ...
    rhoe, Jex, Jey, Jez,... 
    Qbulkex, Qbulkey, Qbulkez, Qenthex, Qenthey, Qenthez, Qhfex, Qhfey, Qhfez, ...
    Pexx, Peyy, Pezz, Pexy, Pexz, Peyz, x, y, z, dx, dy, dz, qom, radius,cyl);
    
    
DUbulkDt = JedotE - Ubulk.*divVe - udivP;
DUthDt = -Uth.*divVe -PgradV;

end    


