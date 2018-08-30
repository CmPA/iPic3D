clear all
generateX=0

load('bimodal_ions.mat');
[Nx,Ny,Nz]=size(fcutoff);
if(generateX)
a=floor(10*fcutoff./max(fcutoff(:)));



X=[];
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
        for ip=1:a(i,j,k)
            X=[X;[i j k]];
        end  
                [ i j k max(size(X))]
        end
    end
end
X2=X+rand(size(X));
plot3(X2(:,1),X2(:,2),X2(:,3)'.','MarkerSize',[1])
save('X3d.mat','X','X2','-mat')
else
    load('X3d.mat')
end


close all
Nc=4
opts = statset('Display','final');
[idx,C] = kmeans(X2,Nc,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

figure;
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

C=C/(Nx-1)*560*2-560

figure;
subplot(2,2,1)
aa=squeeze(mean(fcutoff,3));
imagesc([-560 560],[-560 560],aa')
title('xy')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(mean(fcutoff,2));
imagesc([-560 560],[-560 560],aa')
title('xz')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(mean(fcutoff,1));
imagesc([-560 560],[-560 560],aa')
title('yz')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

print('figure_linear','-dpng')
figure;
subplot(2,2,1)
aa=squeeze(log(mean(fcutoff,3)));
imagesc([-560 560],[-560 560],aa')
title('xy')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(log(mean(fcutoff,2)));
imagesc([-560 560],[-560 560],aa')
title('xz')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(log(mean(fcutoff,1)));
imagesc([-560 560],[-560 560],aa')
title('yz')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy


print('figure_log','-dpng')

[n ncluster']


