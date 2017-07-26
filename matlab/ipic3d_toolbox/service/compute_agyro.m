function [Agyro,Agyro_aunai,Nongyro_swisdak]=compute_Agyro(Bx, By, Bz, ...
 PeXX, PeXY, PeXZ, PeYY, PeYZ, PeZZ ) 


[nx ny nz]= size(PeXX);

for i=1:nx
for iy=1:ny
for k=1:nz

p(1,1)=(PeXX(i,iy,k));
p(1,2)=(PeXY(i,iy,k));
p(1,3)=(PeXZ(i,iy,k));
p(2,2)=(PeYY(i,iy,k));
p(2,3)=(PeYZ(i,iy,k));
p(3,3)=(PeZZ(i,iy,k));
p(2,1)=p(1,2);
p(3,1)=p(1,3);
p(3,2)=p(2,3);

b(1)=(Bx(i,iy,k));
b(2)=(By(i,iy,k));
b(3)=(Bz(i,iy,k));

b=reshape(b,3,1);
b=b./sqrt(sum(b.^2));

%%%%%%%%%%
% Scudder
%%%%%%%%%%

for l=1:3
N1(l,:)=cross(b,p(l,:));
end

for l=1:3
N(:,l)=cross(b,N1(:,l));
end

lamb=sort(eig(N));
lambda1(i,iy,k)=lamb(1);
lambda2(i,iy,k)=lamb(2);
lambda3(i,iy,k)=lamb(3);

Agyro(i,iy,k)= 2*(lamb(3)-lamb(2))/(lamb(3)+lamb(2));


%%%%%%%%%%%
% Aunai
%%%%%%%%%%%

Tr=(p(1,1)+p(2,2)+p(3,3));
Ppar=b'*p*b;
Pper=(Tr-Ppar)/2;
G=eye(3,3)*Pper+(Ppar-Pper)*kron(b,b');
N=p-G;
Agyro_aunai(i,iy,k)=sqrt(sum(N(:).^2))./Tr;

%%%%%%%%%%%
% Swisdak
%%%%%%%%%%%

I2=p(1,1)*p(2,2)+p(1,1)*p(3,3)+p(2,2)*p(3,3);
I2=I2-(p(1,2).^2+p(1,3).^2+p(2,3).^2);
Q=1-4*I2./((Tr-Ppar).*(Tr+3*Ppar));
Nongyro_swisdak(i,iy,k)=sqrt(Q);
% The following formula form Swisdak paper is actually wrong
%Nongyro_aunai(i,k)=sqrt(8*(p(1,2).^2+p(1,3).^2+p(2,3).^2))./(Ppar+2*Pper);



end
end
end


