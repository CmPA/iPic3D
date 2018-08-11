function  [cAx, cAy, cAz, cAv] = compute_curl(x,y,z,Ax,Ay,Az,radius,cyl)
dx=x(1,2,1)-x(1,1,1);
dy=y(2,1,1)-y(1,1,1);
dz=y(1,1,2)-z(1,1,1);

small=1e-10;


    [cAx, cAy, cAz] = curl(x,y,z,permute(Ax,[2 1 3]), permute(Ay, [2 1 3]), permute(Az, [2,1,3]));

cAx = permute(cAx, [2 1 3]);
cAy = permute(cAy, [2 1 3]);
cAz = permute(cAz, [2 1 3]);
end