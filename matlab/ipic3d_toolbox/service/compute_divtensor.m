function  [divAx, divAy, divAz] = compute_divtensor(x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz)
% x is r
% y is z
% z is theta

global cyl Nsm_div

    divAx = compute_div(x,y,z,Axx,Axy,Axz,1,0);
    divAy = compute_div(x,y,z,Axy,Ayy,Ayz,1,0);
    divAz = compute_div(x,y,z,Axz,Ayz,Azz,1,0);

end