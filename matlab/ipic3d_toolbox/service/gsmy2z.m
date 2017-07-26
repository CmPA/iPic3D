function zgsm = gsmy2z(y)
global Ly Zgsmrange
zgsm= y/Ly*(Zgsmrange(2)-Zgsmrange(1))+Zgsmrange(1);
return
