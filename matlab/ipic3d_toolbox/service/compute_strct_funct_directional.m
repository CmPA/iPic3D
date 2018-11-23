function [autocorX,autocorY,autocorZ,volX,volY,volZ] = compute_strct_fucnt_directional(mask,var,var2,var3) 

if(nargin<3)
    var2=zeros(size(var));
end
if(nargin<4)
    var3=zeros(size(var));
end    
[nxn,nyn,nzn]=size(var);
autocorX = zeros(nxn,1);
               
    
volX=zeros(nxn,1);%sum(mask(:));

for j=1:nyn
for k=1:nzn
for i=1:nxn
      if(mask(i,j,k)~=0)
            for ir=0:nxn/2-i
                        if(mask(i+ir,j,k)==0) 
                            break
                        else    
                        autocorX(ir+1) = autocorX(ir+1) + (var(i,j,k) - var(i+ir,j,k)).^2 ...
                        + (var2(i,j,k) - var2(i+ir,j,k)).^2 + (var3(i,j,k) - var3(i+ir,j,k)).^2;
                        volX(ir+1)=volX(ir+1)+1;
                        end
            end
      end
end
end
end


autocorX = autocorX./(volX+1e-10);

    
autocorY = zeros(nyn,1);volY=zeros(nyn,1);
for i=1:nxn
for k=1:nzn
for j=1:nyn
      if(mask(i,j,k)~=0)
            for jr=0:nyn/2-j
                        if(mask(i,j+jr,k)==0) 
                            break
                        else    
                        autocorY(jr+1) = autocorY(jr+1) + (var(i,j,k) - var(i,j+jr,k)).^2 ...
                        + (var2(i,j,k) - var2(i,j+jr,k)).^2 + (var3(i,j,k) - var3(i,j+jr,k)).^2;
                        volY(jr+1)=volY(jr+1)+1;
                        end
            end
      end
end
end
end            




autocorY = autocorY./volY;

    
    
autocorZ = zeros(nzn,1);volZ=zeros(nzn,1);
               
    

for i=1:nxn
for j=1:nyn
for k=1:nzn
      if(mask(i,j,k)~=0)
            for kr=0:nzn/2-k
                        if(mask(i,j,k+kr)==0) 
                            break
                        else    
                        autocorZ(kr+1) = autocorZ(kr+1) + (var(i,j,k) - var(i,j,k+kr)).^2 ...
                         + (var2(i,j,k) - var2(i,j,k+kr)).^2  + (var3(i,j,k) - var3(i,j,k+kr)).^2  ;
                        volZ(kr+1)=volZ(kr+1)+1;
                        end
            end
      end
end
end
end            



autocorZ = autocorZ./volZ;

end
    