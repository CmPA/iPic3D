function  [divAx, divAy, divAz] = compute_divtensor(x,y,z,Axx,Axy,Axz,Ayy,Ayz,Azz,radius,cyl)
% x is r
% y is z
% z is theta

    divAx = compute_div(x,y,z,Axx,Axy,Axz,radius,cyl);
    divAy = compute_div(x,y,z,Axy,Ayy,Ayz,radius,cyl);
    divAz = compute_div(x,y,z,Axz,Ayz,Azz,radius,cyl);

end