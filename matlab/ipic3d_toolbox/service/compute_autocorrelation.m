function [autocor] = compute_autocorrelation(var) 

[nxn,nyn,nzn]=size(var);
autocor = zeros(nxn,nyn,nzn);
               
    

for ir=1:nxn 
    ir
   for jr=1:nyn
      for kr=1:nzn 
                summa =0.0;
            for i=1:nxn-ir
                for j=1:nyn-jr
                    for k=1:nzn-kr
                        summa = summa + var(i,j,k) * var(i+ir,j+jr,k+kr);
                    end
                end
            end
            
        autocor(ir,jr,kr) =summa;
      end
   end
end


    for i=1:nxn 
        for j=1:nyn
            for k=1:nzn
       
            autocor(i,j,k) = autocor(i,j,k)/((nxn-i)*(nyn-j)*(nzn-k)+1.0e-10);
            end
        end
    end
    
end
    