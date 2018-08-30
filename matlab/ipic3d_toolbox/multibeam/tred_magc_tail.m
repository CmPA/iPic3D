clear all
close all



%In nT:
Bx=-0.4189;
By=1.5791;
Bz=12.7405; 

n0=1; % particle per cc


load('uvw-tail.mat')
load('tail-vdf.mat')
fcutoff=permute(fcutoff,[2 1 3]);
fcutoff=smooth3(fcutoff,'box',5);
Np=max(size(up));
X2=zeros(Np,3);
X2(:,1)=up;
X2(:,2)=vp;
X2(:,3)=wp;

Nc=4
opts = statset('Display','final');
[idx,C] = kmeans(X2,Nc,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

figure(1)
c=['r' 'g' 'y' 'm' 'r']
hold on
for ic=1:Nc
%plot(X(idx==ic,1),X(idx==ic,2),[c(ic) '.'],'MarkerSize',12)
plot3(X2(idx==ic,1),X2(idx==ic,2),X2(idx==ic,3),'.','color',rand(1,3),'MarkerSize',1)
cluster=[C(ic,1) C(ic,2) C(ic,3)]
cluster_avg=[mean(X2(idx==ic,1)) mean(X2(idx==ic,2)) mean(X2(idx==ic,3))]
ncluster(ic)=max(size(X2(idx==ic,1)))/max(size(X2))
end
plot3(C(:,1),C(:,2),C(:,3),'kx',...
     'MarkerSize',15,'LineWidth',3) 
%legend('Cluster 1','Cluster 2','Centroids',...
%       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off


U=mean(X2);

M=[C'; ones(1,Nc)];
b=[U'; 1];
n=M\b

conservation=n'*C-U


%C=C/(Nx-1)*560*2-560

%U=U/(Nx-1)*560*2-560


figure(2)
subplot(2,2,1)
aa=squeeze(mean(fcutoff,3));
imagesc([-.06 .06],[-.06 .06],aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(mean(fcutoff,2));
imagesc([-.06 .06],[-.06 .06],aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(mean(fcutoff,1));
imagesc([-.06 .06],[-.06 .06],aa')
title('|| - \perp_2')

hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

print('figure_linear','-dpng')

figure(3)
subplot(2,2,1)
aa=squeeze(log(mean(fcutoff,3)));
imagesc([-.06 .06],[-.06 .06],aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(log(mean(fcutoff,2)));
imagesc([-.06 .06],[-.06 .06],aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(log(mean(fcutoff,1)));
imagesc([-.06 .06],[-.06 .06],aa')
title('|| - \perp_2')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy


print('figure_log','-dpng')


[n ncluster']



figure(4)


dist=[]
for icluster=1:Nc
    dist=[dist sum((XHE2(:,:)-C(icluster,:)).^2,2)];
end
 [dm,im]=min(dist,[],2);
 
for ic=1:Nc
 plot3(XHE2(im==ic,1),XHE2(im==ic,2),XHE2(im==ic,3),'.','color',rand(1,3),'MarkerSize',1)
 hold on


x_mean(ic,:) = (sum(X2(idx==ic,:).*Q(idx==ic),1)+sum(XHE2(im==ic,:).*QHE(im==ic),1)) ...
    ./(sum(Q(idx==ic))+sum(QHE(im==ic)))

x2_mean(ic,:) = (sum((X2(idx==ic,:)-x_mean(ic,:)).^2.*Q(idx==ic),1)+sum((XHE2(im==ic,:)-x_mean(ic,:)).^2.*QHE(im==ic),1)) ...
    ./(sum(Q(idx==ic))+sum(QHE(im==ic)))

x3_mean(ic,:) = (sum((X2(idx==ic,:)-x_mean(ic,:)).^3.*Q(idx==ic),1)+sum((XHE2(im==ic,:)-x_mean(ic,:)).^3.*QHE(im==ic),1)) ...
    ./(sum(Q(idx==ic))+sum(QHE(im==ic)))

Q_mean(ic,:) = x_mean(ic,:).^2*(sum(Q(idx==ic))+sum(QHE(im==ic)));
    
end 
 
x_mean_t = (sum(X2(:,:).*Q,1)+sum(XHE2(:,:).*QHE,1)) ...
    ./(sum(Q)+sum(QHE))

x2_mean_t = (sum((X2(:,:)-x_mean_t).^2.*Q,1)+sum((XHE2(:,:)-x_mean_t).^2.*QHE,1)) ...
    ./(sum(Q)+sum(QHE))

x3_mean_t = (sum((X2(:,:)-x_mean_t).^3.*Q,1)+sum((XHE2(:,:)-x_mean_t).^3.*QHE,1)) ...
    ./(sum(Q)+sum(QHE))

Q_mean_t = (sum((X2(:,:)).^2.*Q,1)+sum((XHE2(:,:)).^2.*QHE,1)) ...
    ./(sum(Q)+sum(QHE))


Q_non_bulk_1beam=(Q_mean_t-x_mean_t.^2)./sum(x_mean_t.^2)
Q_non_bulk_multibeam=(Q_mean_t.*(sum(Q)+sum(QHE))-sum(Q_mean))./(sum(Q)+sum(QHE))./sum(x_mean_t.^2)
