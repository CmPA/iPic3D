function [J] = smooth_ipic2d(Jin)
J=Jin;
N=6;
for i=1:N 
  value=0.5;
  if(mod(i,2)==1) 
      value=0;
  end
  alfa=(1-value)/6;
  J(2:end-1,2:end-1)=(2*alfa+value)*J(2:end-1,2:end-1)+alfa*(J(1:end-2,2:end-1)+J(2:end-1,1:end-2)+J(3:end,2:end-1)+J(2:end-1,3:end)); 
end
