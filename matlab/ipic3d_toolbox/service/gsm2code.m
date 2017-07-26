function [x y z ix iy iz ipx ipy ipz ip] = gsm2code(xgsm, ygsm, zgsm)
global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN Nx Ny Nz
%
%	Method to convert GSM coordinates to code coordinates, 
%   grid position and processor id 
%


x= (Xgsmrange(2)-xgsm) * Lx /(Xgsmrange(2)-Xgsmrange(1));

y= (zgsm - Zgsmrange(1)) * Ly/ (Zgsmrange(2)-Zgsmrange(1));

z= (ygsm - Ygsmrange(1)) *Lz/ (Ygsmrange(2)-Ygsmrange(1));

ix= 0+ floor(x/dx); iy =0 +floor(y/dy);  iz = 0+ floor(z/dz);

ipx = 1+ floor(ix/(Nx/XLEN)); ipy = 1+ floor(iy/(Ny/YLEN)) ; ipz = 1+ floor(iz/(Nz/ZLEN))

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

return
