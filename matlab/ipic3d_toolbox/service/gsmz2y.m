function ygsm = gsmz2y(z)
global Lz Ygsmrange

ygsm= z/Lz*(Ygsmrange(2)-Ygsmrange(1))+Ygsmrange(1);
%ygsm=-ygsm;
return
