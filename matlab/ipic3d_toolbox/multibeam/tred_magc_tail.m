clear all
close all



%In nT:
Bx=-0.4189;
By=1.5791;
Bz=12.7405; 

n0=1; % particle per cc


load('uvw-tail.mat')
load('tail-vdf.mat')
%fcutoff=permute(fcutoff,[2 1 3]);
fcutoff=smooth3(fcutoff,'box',5);
Np=max(size(up));
X2=zeros(Np,3);
X2(:,1)=vp;
X2(:,2)=up;
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

load ../service/gist_ncar.mat


figure(2)
aa=squeeze(mean(fcutoff,3));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_x','fontsize',[20])
ylabel('v_y','fontsize',[20])
set(gca,'fontsize',[20])
hold on
plot(C(:,1),C(:,2),'mx',...
     'MarkerSize',15,'LineWidth',1) 
axis xy
axis equal
axis tight
colormap(parula)
print('figure_linear_1','-dpng','-r300')

figure(3)
aa=squeeze(mean(fcutoff,2));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_x','fontsize',[20])
ylabel('v_z','fontsize',[20])
set(gca,'fontsize',[20])
hold on
plot(C(:,1),C(:,3),'mx',...
     'MarkerSize',15,'LineWidth',1) 
axis xy
axis equal
axis tight
colormap(parula)
print('figure_linear_2','-dpng','-r300')

figure(4)
aa=squeeze(mean(fcutoff,1));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_y','fontsize',[20])
ylabel('v_z','fontsize',[20])
set(gca,'fontsize',[20])

hold on
plot(C(:,2),C(:,3),'mx',...
     'MarkerSize',15,'LineWidth',1) 
axis xy
axis equal
axis tight
colormap(parula)
print('figure_linear_3','-dpng','-r300')

figure(5)
subplot(2,2,1)
aa=squeeze(log(mean(fcutoff,3)));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_x','fontsize',[20])
ylabel('v_y','fontsize',[20])
hold on
plot(C(:,1),C(:,2),'wx',...
     'MarkerSize',6,'LineWidth',1) 
axis xy
colormap(parula)

subplot(2,2,2)
aa=squeeze(log(mean(fcutoff,2)));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_x','fontsize',[20])
ylabel('v_z','fontsize',[20])
hold on
plot(C(:,1),C(:,3),'wx',...
     'MarkerSize',6,'LineWidth',1) 
axis xy
colormap(parula)

subplot(2,2,3)
aa=squeeze(log(mean(fcutoff,1)));
imagesc([-.06 .06],[-.06 .06],aa')
xlabel('v_y','fontsize',[20])
ylabel('v_z','fontsize',[20])
hold on
plot(C(:,2),C(:,3),'wx',...
     'MarkerSize',6,'LineWidth',1) 
axis xy
colormap(parula)

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

input.tableColLabels = {'x','y','z'};
input.tableRowLabels = {'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'};


% LaTex table caption:
input.tableCaption = 'Cluster Centers';
% LaTex table label:
input.tableLabel = 'cluster_center';
input.data = C
latex = latexTable(input)

input.tableColLabels = {'x','y','z','mod'};
input.tableRowLabels = {'1 cluster','4 clusters'};


[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fprintf(fid,'\n');

x_mean_t = (sum(X2(:,:).*Q,1)) ...
    ./(sum(Q));

Ubar=[x_mean_t ;(x_mean'*n)']
Ubar=[Ubar sqrt(sum(Ubar.^2,2))]
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
Uth=[Uth sum(Uth,2)]

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
HF=[x3_mean_t/2; (x3_mean'*n)'/2]
HF=[HF sqrt(sum(HF.^2,2))]

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


fid=fopen('MyLatex_SI.tex','w');
input.dataFormat = {'%e'};
input.tableColLabels = {'x','y','z','mod'};
input.tableRowLabels = {'1 cluster','4 clusters'};


%Convert to SI
% 0.5*n*me*256*U^2*c^2 (*1e9 to have it in nW/m^3)
c=299792458;
n_mi_c2=0.01*1e6*9.1e-31*(c)^2*256;


input.data = Ubar.^2/2*n_mi_c2*1e9

% LaTex table caption:
input.tableCaption = 'Bulk Energy $nJ/m^3$ ';
% LaTex table label:
input.tableLabel = 'bulk_energy_si';

latex = latexTable(input)

[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fprintf(fid,'\n');

% Uth*c^2*n*me*256


% LaTex table caption:
input.tableCaption = 'Thermal Energy $nJ/m^3$';

% LaTex table label:
input.tableLabel = 'thermal_energy_si';



input.data = Uth*n_mi_c2*1e9
latex = latexTable(input)

[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fprintf(fid,'\n');

% HF*c^2*n*me*256


% LaTex table caption:
input.tableCaption = 'Heat Flux $mW/m^2$';

% LaTex table label:
input.tableLabel = 'heat_flux_si';

input.data = HF*n_mi_c2*1e3*c
latex = latexTable(input)

[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved as ''MyLatex.tex'' in your working directory\n');

!cp *.tex ~/Dropbox/Science/tex/marty-multibeam
!/usr/local/bin/mogrify -trim *.png
!cp *.png ~/Dropbox/Science/tex/marty-multibeam