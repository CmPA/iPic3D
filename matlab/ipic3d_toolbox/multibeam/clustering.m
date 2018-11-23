function [hist_NC,hist_frac,x2_mean_t,Q_tot_t] = clustering()
clear all
close all

generateX=true
separate_he=false
graphicose=true


hist_NC=[];
hist_frac=[];

%In nT:
Bx=-0.4189;
By=1.5791;
Bz=12.7405; 

n0=1; % particle per cc


% reading in the FPI distribution

load('bimodal_ions_FAC.mat');
[Nx,Ny,Nz]=size(fcutoff);
[xg,yg,zg]=ndgrid(0:Nx-1,0:Ny-1,0:Nz-1);
xg=xg/(Nx-1)*560*2-560;
yg=yg/(Ny-1)*560*2-560;
zg=zg/(Nz-1)*560*2-560;
Uxdata=mean(mean(mean(fcutoff.*xg)))./mean(mean(mean(fcutoff)));
Uydata=mean(mean(mean(fcutoff.*yg)))./mean(mean(mean(fcutoff)));
Uzdata=mean(mean(mean(fcutoff.*zg)))./mean(mean(mean(fcutoff)));
Udata=[Uxdata Uydata Uzdata];

a=fcutoff;

if(generateX)
    

if(separate_he)
    [XHE2 QHE] = generate_he(fcutoff);
    a=(10*fcutoff./max(fcutoff(:)));
    ii=a<1;
    a(ii)=0;
end

% generating particles with it - rejection method
[X2,Q] = generate_rejection(a);

% generating particles with it - rejection method
%[X2,Q] = generate_cdf(a);

% generating particles with it - integer steps method
%[X2,Q] = generate_cut(a);

if(separate_he)
plot3(XHE2(:,1),XHE2(:,2),XHE2(:,3)'.','MarkerSize',[1])
save('X3dHE_FAC.mat','XHE2','QHE','-mat')
else
        XHE2=[];
        QHE=[];
end

save('X3d_FAC.mat','X2','Q','-mat')

else
    load('X3d_FAC.mat')
    if(separate_he)
    load('X3dHE_FAC.mat')
    else
        XHE2=[];
        QHE=[];
    end
end

for Nc=2:10
opts = statset('Display','final');
[idx,C] = kmeans(X2,Nc,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);



%c=['r' 'g' 'y' 'm' 'r']
for ic=1:Nc
%plot(X(idx==ic,1),X(idx==ic,2),[c(ic) '.'],'MarkerSize',12)

if(graphicose)
figure(1)
plot3(X2(idx==ic,1),X2(idx==ic,2),X2(idx==ic,3),'.','color',rand(1,3),'MarkerSize',1)
hold on
end
cluster=[C(ic,1) C(ic,2) C(ic,3)]
cluster_avg=[mean(X2(idx==ic,1)) mean(X2(idx==ic,2)) mean(X2(idx==ic,3))]
ncluster(ic)=max(size(X2(idx==ic,1)))/max(size(X2))
end
if(graphicose)
plot3(C(:,1),C(:,2),C(:,3),'kx',...
     'MarkerSize',15,'LineWidth',3) 
%legend('Cluster 1','Cluster 2','Centroids',...
%       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off
end

U=mean(X2);

M=[C'; ones(1,Nc)];
b=[U'; 1];
n=M\b

conservation=n'*C-U


%C=C/(Nx-1)*560*2-560

%U=U/(Nx-1)*560*2-560


if(graphicose)
figure(2)
subplot(2,2,1)
aa=squeeze(mean(fcutoff,3));
imagesc([-560 560],[-560 560],aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(mean(fcutoff,2));
imagesc([-560 560],[-560 560],aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(mean(fcutoff,1));
imagesc([-560 560],[-560 560],aa')
title('|| - \perp_2')

hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

print([num2str(Nc) 'figure_linear'],'-dpng')

figure(3)
subplot(2,2,1)
aa=squeeze(log(mean(fcutoff,3)));
imagesc([-560 560],[-560 560],aa')
title('\perp_1 - ||')
hold on
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,2)
aa=squeeze(log(mean(fcutoff,2)));
imagesc([-560 560],[-560 560],aa')
title('\perp_1 - \perp_2')
hold on
plot(C(:,1),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy

subplot(2,2,3)
aa=squeeze(log(mean(fcutoff,1)));
imagesc([-560 560],[-560 560],aa')
title('|| - \perp_2')
hold on
plot(C(:,2),C(:,3),'kx',...
     'MarkerSize',3,'LineWidth',1) 
axis xy


print([num2str(Nc) 'figure_log'],'-dpng')
end

[n ncluster']


for ic=1:Nc
% Moment one of each cluster
x_mean(ic,:) = sum(X2(idx==ic,:).*Q(idx==ic),1)./sum(Q(idx==ic));

% Moment two of each cluster
x2_mean(ic,:) = sum((X2(idx==ic,:)-x_mean(ic,:)).^2.*Q(idx==ic),1);

% Bulk Energy  of each cluster
Q_mean(ic,:) = x_mean(ic,:).^2*sum(Q(idx==ic));
    
end 
 


if(separate_he)
    x_mean_t = (sum(X2.*Q,1)+sum(XHE2.*QHE,1)) ...
    ./(sum(Q)+sum(QHE));

x2_mean_t = sum((X2(:,:)-x_mean_t).^2.*Q,1)+sum((XHE2-x_mean_t).^2.*QHE,1);


Q_mean_t = x_mean_t.^2.*(sum(Q)+sum(QHE));
Q_tot_t = sum(X2.^2.*Q,1)+sum(XHE2.^2.*QHE,1);
conservation_energy=sum(x2_mean+Q_mean)+sum(XHE2.^2.*QHE,1)-Q_tot_t

else
    size(Q)
    size(X2)
    x_mean_t = sum(X2.*Q,1)./sum(Q);

x2_mean_t = sum((X2-x_mean_t).^2.*Q,1);


Q_mean_t = x_mean_t.^2.*sum(Q);
Q_tot_t= sum(X2.^2.*Q);
conservation_energy=sum(x2_mean+Q_mean)-Q_tot_t

end


hist_NC=[hist_NC;Nc];
hist_frac=[hist_frac;sum(sum(x2_mean))./sum(Q_tot_t)];
end
figure(10)
plot(hist_NC,hist_frac)
fit((hist_NC),(hist_frac),'power1')
end

function [X2,Q] = generate_rejection(fcutoff)
[Nx,Ny,Nz]=size(fcutoff);


a=(10*fcutoff./max(fcutoff(:)));


X=zeros(1000000,3);ip=1
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            Nv=floor(a(i,j,k));
            if(a(i,j,k)-Nv>rand(1,1))
                Nv=Nv+1;
            end    
            %[a(i,j,k)-Nv, Nv]
            %X=[X;repmat([i j k],Nv,1)];
        for ip2=1:Nv
            X(ip,:)=[i j k];
            ip=ip+1;
        end  
        end
    end
            [ i  ip]
end
X2=X+rand(size(X));
X2=X2(1:ip-1,:)/(Nx-1)*560*2-560;
Q=ones(size(X2(:,1)));
end


function [X,Q] = generate_cdf(fcutoff)
% from David Newman
[nx,ny,nz]=size(fcutoff);
vmax=560                             %for small box
%vmax=560*3.75;                       %for big box

dv=2*vmax/(nx-1);
v1D=-vmax:dv:vmax;
[vvx,vvy,vvz]=meshgrid(v1D,v1D,v1D);
f1D=fcutoff(:);
fmin=min(f1D(f1D>0));

%Distribution cannot equal zero for method to work
%so set zeros (if the exist) to small fraction of 
%smallest nonzero value.
f1D(f1D==0)=fmin/1000;      
vx1D=vvx(:)'; vy1D=vvy(:)'; vz1D=vvz(:)';

Np=1e5;          %Number of randomly distributed particles

Ng=size(f1D,1);
ranarr=rand(4,Np);
fcum=cumsum(f1D);
fcum=Ng*fcum/fcum(Ng);
Pg=interp1(fcum',1:Ng,Ng*ranarr(1,:));
Pg=1+floor(Pg);
xp=vx1D(Pg)+dv*ranarr(2,:)-dv/2;
yp=vy1D(Pg)+dv*ranarr(3,:)-dv/2;
zp=vz1D(Pg)+dv*ranarr(4,:)-dv/2;

means=[mean(xp) mean(yp) mean(zp)]
%figure(1)
%scatter3(xp(1:1000:Np),yp(1:1000:Np),zp(1:1000:Np));
%daspect([1 1 1])

X= [yp' xp' zp'];
size(X)
Q=ones(size(X(:,1)));

end

function [X2,Q] = generate_cut(fcutoff)
[Nx,Ny,Nz]=size(fcutoff);
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
X2=X2/(Nx-1)*560*2-560;
Q=ones(size(X2(:,1)));
end


function [XHE2 QHE] = generate_he(fcutoff)
[Nx,Ny,Nz]=size(fcutoff);
a2=(10*fcutoff./max(fcutoff(:)));
ii=a2>=1;
a2(ii)=0;
ii=a2>0.1;
QHE=a2(ii);
NpHE=max(size(QHE));
[x1,x2,x3]=ndgrid(1:Nx,1:Ny,1:Nz);
XHE=[x1(ii) x2(ii) x3(ii)];

XHE2=XHE+rand(size(XHE));
XHE2=XHE2/(Nx-1)*560*2-560;

end