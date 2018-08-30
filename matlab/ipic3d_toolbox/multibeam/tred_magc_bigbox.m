clear all
close all

generateX=1

%In nT:
Bx=-0.4189;
By=1.5791;
Bz=12.7405; 

n0=1; % particle per cc

load('bimodal_ions_FAC_bigbox.mat');
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
        end
    end
            [ i  max(size(X))]
end
X2=X+rand(size(X));
Q=ones(size(X2(:,1)));
plot3(X2(:,1),X2(:,2),X2(:,3),'.','MarkerSize',[1])
save('X3d_FAC.mat','X2','Q','-mat')


a2=(10*fcutoff./max(fcutoff(:)));
ii=a2>=1;
a2(ii)=0;
ii=a2>0.1;
QHE=a2(ii);
NpHE=max(size(QHE));
[x1,x2,x3]=ndgrid(1:Nx,1:Ny,1:Nz);
XHE=[x1(ii) x2(ii) x3(ii)];

XHE2=XHE+rand(size(XHE));
plot3(XHE2(:,1),XHE2(:,2),XHE2(:,3)'.','MarkerSize',[1])
save('X3dHE_FAC.mat','XHE2','QHE','-mat')
else
    load('X3d_FAC.mat')
    load('X3dHE_FAC.mat')
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


C=C/(Nx-1)*560*2-560;
C=C*3.75

U= U/(Nx-1)*560*2-560;
U=U*3.75

close all
figure;
subplot(2,2,1)
aa=squeeze(mean(fcutoff,3));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(mean(fcutoff,2));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(mean(fcutoff,1));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('|| - \perp_2')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

print('figure_linear','-dpng')

figure;
subplot(2,2,1)
aa=squeeze(log(mean(fcutoff,3)));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(log(mean(fcutoff,2)));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(log(mean(fcutoff,1)));
imagesc([-560 560]*3.75,[-560 560]*3.75,aa')
title('|| - \perp_2')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy


print('figure_log','-dpng')


[n ncluster']


