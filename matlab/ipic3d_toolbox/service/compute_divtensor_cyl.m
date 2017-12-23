function  [divAx, divAy, divAz, divAx1, divAx2] = compute_div_tensor_cyl(x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz)
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
    divAx1 = divergence(x,y,z,x.*permute(smooth3Dnew(Axx,Nsm_div),[2 1 3]), ...
        x.*permute(smooth3Dnew(Axy,Nsm_div),[2 1 3]), permute(smooth3Dnew(Axz,Nsm_div),[2 1 3]))./x;
    divAx2 = - permute(smooth3Dnew(Azz,Nsm_div),[2 1 3])./x;
    divAy = divergence(x,y,z,x.*permute(smooth3Dnew(Axy,Nsm_div),[2 1 3]), ...
        x.*permute(smooth3Dnew(Ayy,Nsm_div),[2 1 3]), permute(smooth3Dnew(Ayz,Nsm_div),[2 1 3]))./x;
    divAz = divergence(x,y,z,x.*permute(smooth3Dnew(Axz,Nsm_div),[2 1 3]), ...
        x.*permute(smooth3Dnew(Ayz,Nsm_div),[2 1 3]), permute(smooth3Dnew(Azz,Nsm_div),[2 1 3]))./x ...
        + permute(smooth3Dnew(Axz,Nsm_div),[2 1 3])./x;
  
divAx = permute(divAx1+divAx2, [2 1 3]);
divAy = permute(divAy, [2 1 3]);
divAz = permute(divAz, [2 1 3]);
divAx1 = permute(divAx1, [2 1 3]);
divAx2 = permute(divAx2, [2 1 3]);
end