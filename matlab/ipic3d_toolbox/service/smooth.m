function [J] = smooth(Jin,N)
J=Jin;
for i=1:N 
J(2:end-1,2:end-1)=J(2:end-1,2:end-1)*.5+(J(1:end-2,2:end-1)+J(2:end-1,1:end-2)+J(3:end,2:end-1)+J(2:end-1,3:end))/8; 
end
