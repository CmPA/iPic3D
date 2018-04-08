function [h1,h2] = pdf_fluctuations( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all
%B=A-smooth3(A,'box',5);
data=A(:)-mean(A(:)); %data=randn(size(A));
sigma=std(data(:));
ii=abs(data)>sigma/3;
data=data(ii);

%data=data(:)-mean(data(:));
data=data./sigma;
std(data)
edges=-20:.01:20;
h1=hist(data,edges);
h3=exp(-edges.^2/2);h3=h3./sum(h3).*sum(h1)./max(h1);
%[maxval, ii]=max(h1);
%h1(ii)=0;
semilogy(edges,h1./max(h1),edges,h3)
ylim([1./max(h1) 10])
print -dpng figura_pdf
end

