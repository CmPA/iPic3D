function [J] = smooth_1d(Jin,N)
J=Jin;
for i=1:N 
J(:,2:end-1)=J(:,2:end-1)*.5+.25*(J(:,1:end-2)+J(:,3:end)); 
end
