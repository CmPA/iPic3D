function  divA = compute_div(x,y,z,Ax,Ay,Az)
dx=x(1,2,1)-x(1,1,1);
dy=y(2,1,1)-y(1,1,1);
dz=y(1,1,2)-z(1,1,1);
global cyl 
small=1e-10
if(cyl)
    [divA1,tmp1,tmp2] = gradient(x.*permute(Ax,[2 1 3]),dx,dy,dz);
    [tmp1,divA2,tmp2] = gradient(permute(Ay,[2 1 3]),dx,dy,dz);
    divA=divA1./(x+small)+divA2;
    divA(:,1,:)=0;
else    
    divA = divergence(x,y,z,permute(Ax,[2 1 3]), permute(Ay, [2 1 3]), permute(Az, [2,1,3]));
end
divA = permute(divA, [2 1 3]);
end