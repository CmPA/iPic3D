function [J] = smoothbc(Jin,N)
J=Jin;
for i=1:N 
J(2:end-1,2:end-1)=J(2:end-1,2:end-1)*.5+(J(1:end-2,2:end-1)+J(2:end-1,1:end-2)+J(3:end,2:end-1)+J(2:end-1,3:end))/8; 
J(1,2:end-1)=J(1,2:end-1)*.5+(J(2,2:end-1)+J(1,1:end-2)+J(1,3:end))/6;
J(end,2:end-1)=J(end,2:end-1)*.5+(J(end-1,2:end-1)+J(end,1:end-2)+J(end,3:end))/6;
J(2:end-1,1)=J(2:end-1,1)*.5+(J(2:end-1,2)+J(1:end-2,1)+J(3:end,1))/6;
J(2:end-1,end)=J(2:end-1,end)*.5+(J(2:end-1,end-1)+J(1:end-2,end)+J(3:end,end))/6;
%coreners missing
end
