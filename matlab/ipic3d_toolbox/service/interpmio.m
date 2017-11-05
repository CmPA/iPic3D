function Tpic =  interpmio(x,y,z,T,Xpic,Ypic,Zpic);

Tpic=interpn(permute(x,[3 2 1]),permute(y,[3 2 1]), permute(z,[3 2 1]), permute(T,[3 2 1]), Xpic, Ypic, Zpic);
