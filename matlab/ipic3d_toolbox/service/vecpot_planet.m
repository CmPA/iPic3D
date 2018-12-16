function ay=vecpot_planet(bx,by,dx,dy)


[nx,nz]=size(bx);

nzmezzo=nz%ceil(nz/2);

ay=zeros(size(bx));
for ind=2:nx
	ay(ind,nzmezzo)=ay(ind-1,nzmezzo)+trapz(by(ind-1:ind,nzmezzo),1)*dx;
end
for ind=nzmezzo+1:nz
%   ind;
   ay(:,ind)=ay(:,ind-1)-trapz(bx(:,ind-1:ind),2)*dy;
end
for ind=nzmezzo-1:-1:1
   ay(:,ind)=ay(:,ind+1)+trapz(bx(:,ind:ind+1),2)*dy;
end

% nxmezzo=nx%ceil(nx/2);
% 
% ay=zeros(size(bx));
% for ind=2:nz
% 	ay(nxmezzo,ind)=ay(nxmezzo,ind-1)-(bx(nxmezzo,ind-1)+bx(nxmezzo,ind))*dy/2;
% end
% for ind=nxmezzo+1:nx
% %   ind;
%    %ay(ind,:)=ay(ind-1,:)+(by(ind-1,:)+by(ind-1,:))*dx/2;
% end
% for ind=nxmezzo-1:-1:1
%    ay(ind,:)=ay(ind+1,:)-(by(ind,:)+by(ind+1,:))*dx/2;
% end
% 
% return
% for i in range(1,nx):
%       az[i][nymezzo] = az[i-1][nymezzo]- (bx[i-1][nymezzo]+bx[i][nymezzo])*dy/2.0
% 
%     for ind in range(nymezzo+1,ny):
%         for j in range(0,nx):
%             az[j,ind]=az[j][ind-1]+ (by[j][ind-1] + by[j][ind-1])*dx/2.0
%     
%     for ind in range(nymezzo-1,-1,-1):
%         for j in range(0,nx):
%             az[j,ind]=az[j][ind+1]- (by[j][ind+1] + by[j][ind])*dx/2.0
