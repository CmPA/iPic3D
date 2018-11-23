function [autocorX,autocorY,autocorZ] = compute_autocorrelation_directional(var) 

[nxn,nyn,nzn]=size(var);
autocorX = zeros(nxn,1);
               
    

for ir=0:nxn-1
                summa =0.0;
            for i=1:nxn-ir
                for j=1:nyn
                    for k=1:nzn
                        summa = summa + var(i,j,k) * var(i+ir,j,k);
                    end
                end
            end
            
        autocorX(ir+1) =summa;

end


    for i=1:nxn 
            autocorX(i) = autocorX(i)/((nxn-i+1)*(nyn)*(nzn));
    end
    
    
    autocorY = zeros(nyn,1);
               
    

for jr=0:nyn-1 
                summa =0.0;
            for i=1:nxn
                for j=1:nyn-jr
                    for k=1:nzn
                        summa = summa + var(i,j,k) * var(i,j+jr,k);
                    end
                end
            end
            
        autocorY(jr+1) =summa;

end


    for j=1:nyn 
            autocorY(j) = autocorY(j)/((nxn)*(nyn-j+1)*(nzn));
    end
    
    
    autocorZ = zeros(nzn,1);
               
    

for kr=0:nzn-1 
                summa =0.0;
            for i=1:nxn
                for j=1:nyn
                    for k=1:nzn-kr
                        summa = summa + var(i,j,k) * var(i,j,k+kr);
                    end
                end
            end
            
        autocorZ(kr+1) =summa;

end


    for k=1:nzn 
            autocorZ(k) = autocorZ(k)/((nxn)*(nyn)*(nzn-k+1));
    end
end
    