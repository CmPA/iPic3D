function ath=vecpot_cyl(r,z,br,bz)

[nz,nr]=size(br);
nzmezzo=ceil(nz/2);

bzr=bz(nzmezzo,:).*r
ath=zeros(size(br));
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
