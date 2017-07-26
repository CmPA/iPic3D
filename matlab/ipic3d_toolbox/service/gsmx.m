function xgsm = gsmx(x)
global Lx Xgsmrange
xgsm=-x/Lx*(Xgsmrange(2)-Xgsmrange(1))+Xgsmrange(2);
return
