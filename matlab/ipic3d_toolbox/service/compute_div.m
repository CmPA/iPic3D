function  divA = compute_div(x,y,z,Ax,Ay,Az)
dx=x(1,2,1)-x(1,1,1);
dy=y(2,1,1)-y(1,1,1);
dz=y(1,1,2)-z(1,1,1);
global cyl Nsm_div
small=1e-10
Nsm_div
if(cyl==1)
    method='gaussian'
radius=5;
    %divA = divergence(x,y,z,x.*permute(smooth3(Ax,method,radius),[2 1 3]), ...
    %    x.*permute(smooth3(Ay,method,radius),[2 1 3]), x.*permute(smooth3(Az,method,radius),[2 1 3]));
    divA = divergence(x,y,z,x.*permute(smooth3Dnew(Ax,Nsm_div),[2 1 3]), ...
        x.*permute(smooth3Dnew(Ay,Nsm_div),[2 1 3]), x.*permute(smooth3Dnew(Az,Nsm_div),[2 1 3]));
    %divA(:,1,:)=0;
elseif(cyl==2)
    [divA1,tmp1,tmp2] = gradient(x.*permute(smooth3Dnew(Ax,Nsm_div),[2 1 3]),dx,dy,dz);
    [tmp1,divA2,tmp2] = gradient(permute(smooth3Dnew(Ay,Nsm_div),[2 1 3]),dx,dy,dz);
    divA=divA1./(x+small)+divA2;
    divA(:,1,:)=0;
    divA = divergence(x,y,z,permute(Ax,[2 1 3]), permute(Ay, [2 1 3]), permute(Az, [2,1,3]));
end
divA = permute(divA, [2 1 3]);
end