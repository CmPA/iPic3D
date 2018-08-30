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



 Q=ones(size(X2));
for ic=1:Nc


x_mean(ic,:) = (sum(X2(idx==ic,:).*Q(idx==ic),1)) ...
    ./(sum(Q(idx==ic)));

x2_mean(ic,:) = (sum((X2(idx==ic,:)-x_mean(ic,:)).^2.*Q(idx==ic),1)) ...
    ./(sum(Q(idx==ic)));

x3_mean(ic,:) = (sum((X2(idx==ic,:)-x_mean(ic,:)).^3.*Q(idx==ic),1)) ...
    ./(sum(Q(idx==ic)));

Q_mean(ic,:) = x_mean(ic,:).^2*(sum(Q(idx==ic)));
    
end 

fid=fopen('MyLatex.tex','w');
input.dataFormat = {'%e'};

x_mean_t = (sum(X2(:,:).*Q,1)) ...
    ./(sum(Q));

Ubar=[x_mean_t ;(x_mean'*n)']
input.data = Ubar

% LaTex table caption:
input.tableCaption = 'Mean Velocity ';

% LaTex table label:
input.tableLabel = 'mean_velocity';

latex = latexTable(input)


[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fprintf(fid,'\n');


x2_mean_t = (sum((X2(:,:)-x_mean_t).^2.*Q,1)) ...
    ./(sum(Q));
Uth=[x2_mean_t/2 ; (x2_mean'*n)'/2]

% LaTex table caption:
input.tableCaption = 'Thermal Energy';

% LaTex table label:
input.tableLabel = 'thermal_energy';



input.data = Uth
latex = latexTable(input)

[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fprintf(fid,'\n');

x3_mean_t = (sum((X2(:,:)-x_mean_t).^3.*Q,1)) ...
    ./(sum(Q));
HF=[x3_mean_t; (x3_mean'*n)']

% LaTex table caption:
input.tableCaption = 'Heat Flux';

% LaTex table label:
input.tableLabel = 'heat_flux';

input.data = HF
latex = latexTable(input)

[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end



fclose(fid);
fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');



