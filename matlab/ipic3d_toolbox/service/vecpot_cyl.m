function ath=vecpot_cyl(r,z,br,bz)

[nz,nr]=size(br);
nzmezzo=ceil(nz/2);
nrmezzo=ceil(nr/10);

bzr=bz(nzmezzo,:).*r;
ath=zeros(size(br));
rath2=zeros(size(br));
rath=zeros(1,nr);

for ind=2:nr;
	rath(ind)=rath(ind-1)+trapz(r(ind-1:ind),bzr(ind-1:ind));
end
ath(nzmezzo,:)=rath./r;
for ind=nzmezzo+1:nz
%   ind;
   ath(ind,:)=ath(ind-1,:)-trapz(z(ind-1:ind),br(ind-1:ind,:));
end
for ind=nzmezzo-1:-1:1
   ath(ind,:)=ath(ind+1,:)+trapz(z(ind:ind+1),br(ind:ind+1,:));
end

% for iz=1:nz
%     bzr(iz,:)=bz(iz,:).*r;
% end 
% [size(bz) size(bzr)]
% rath2(:,nrmezzo)=ath(:,nrmezzo).*r(nrmezzo);
% for ir=nrmezzo+1:nr
%     for iz=1:nz
%         rath2(iz,ir)=rath2(iz,ir-1)+sum(bzr(iz,ir-1:ir))/2*(r(ir)-r(ir-1));
%     end
% end    
% for ir=nrmezzo-1:-1:2
%     for iz=1:nz
%         rath2(iz,ir)=rath2(iz,ir+1)-sum(bzr(iz,ir:ir+1))/2*(r(ir+1)-r(ir));
%     end
% end
% for iz=1:nz
%     ath2(iz,2:nr)=rath2(iz,2:nr)./r(2:nr);
% end   
% %ath=ath2;
% ath=.5*(2*ath+0*ath2);
