function  [AgradBx, AgradBy, AgradBz, AgradBx1, AgradBx2] = compute_AgradB_cyl(x,y,z,Ax,Ay,Az,Bx,By,Bz)
dx=x(1,2,1)-x(1,1,1)
dy=y(2,1,1)-y(1,1,1)
dz=z(1,1,2)-z(1,1,1)
% x is r
% y is z
% z is theta
global cyl Nsm_div
small=1e-10
Nsm_div
if(cyl==1)
    method='gaussian'
radius=5;
    %divA = divergence(x,y,z,x.*permute(smooth3(Ax,method,radius),[2 1 3]), ...
    %    x.*permute(smooth3(Ay,method,radius),[2 1 3]), x.*permute(smooth3(Az,method,radius),[2 1 3]));
[tmpx,tmpy,tmpz] = gradient(permute(smooth3Dnew(Bx,Nsm_div),[2 1 3]),dx,dy,dz);
AgradBx1=permute(smooth3Dnew(Ax,Nsm_div),[2 1 3]).*tmpx+permute(smooth3Dnew(Ay,Nsm_div),[2 1 3]).*tmpy+permute(smooth3Dnew(Az,Nsm_div),[2 1 3])./x.*tmpz;
AgradBx2=-permute(smooth3Dnew(Az,Nsm_div),[2 1 3]).*permute(smooth3Dnew(Bz,Nsm_div),[2 1 3])./x;


[tmpx,tmpy,tmpz] = gradient(permute(smooth3Dnew(By,Nsm_div),[2 1 3]),dx,dy,dz);
AgradBy=permute(smooth3Dnew(Ax,Nsm_div),[2 1 3]).*tmpx+permute(smooth3Dnew(Ay,Nsm_div),[2 1 3]).*tmpy+permute(smooth3Dnew(Az,Nsm_div),[2 1 3])./x.*tmpz;

[tmpx,tmpy,tmpz] = gradient(permute(smooth3Dnew(Bz,Nsm_div),[2 1 3]),dx,dy,dz);
AgradBz=permute(smooth3Dnew(Ax,Nsm_div),[2 1 3]).*tmpx+permute(smooth3Dnew(Ay,Nsm_div),[2 1 3]).*tmpy+permute(smooth3Dnew(Az,Nsm_div),[2 1 3])./x.*tmpz;
AgradBz=AgradBz+permute(smooth3Dnew(Az,Nsm_div),[2 1 3]).*permute(smooth3Dnew(Bx,Nsm_div),[2 1 3])./x;
end
AgradBx = permute(AgradBx1+AgradBx2, [2 1 3]);
AgradBy = permute(AgradBy, [2 1 3]);
AgradBz = permute(AgradBz, [2 1 3]);
AgradBx1 = permute(AgradBx1, [2 1 3]);
AgradBx2 = permute(AgradBx2, [2 1 3]);
end