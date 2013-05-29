
function GEMparams()

TiTe = 5;
mime = 25.;
%;.0097;
B0 = .10;
lambda = .5;

[vti,vte,iondrift,elcdrift] = get_velocities(TiTe,mime,B0,lambda);
vte
vti
elcdrift
iondrift

dt_omegapi = .25;
duration_omegagi = 40.;
%omegapi_over_omegagi = 1/B0;
num_cycles = duration_omegagi/(B0*dt_omegapi)

end

%nbn0 = 0. %.2;

function [vti,vte,iondrift,elcdrift] = get_velocities(TiTe,mime,B0,lambda);

c = 1;
e = 1;
mi = 1.;

me = 1/mime;
n0 = 1/(4*pi);

T = B0^2/(8*pi)/n0; % (Ti+Te)

Ti = TiTe/(1+TiTe)*T;
Te = 1/(1+TiTe) * T;
vti = sqrt(Ti/mi);
vte = sqrt(Te/me);
vtecheck = vti*sqrt(mime/TiTe);

reldrift = -c*B0/(e*lambda*n0*4*pi);

iondrift = reldrift/(mime+1);
elcdrift = -reldrift*mime/(mime+1);

end
