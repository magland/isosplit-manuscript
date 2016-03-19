function FIG_decision_boundaries

mfile_path=fileparts(mfilename('fullpath'));

close all;

rng(1);

fontsize=16;

simoptions.centers={[0,0],[3.5,4.2],[-1,4.8]};
simoptions.pops={1000,400,400};
simoptions.shapes={[2,2,0],[0.5,0.5,0],[1.5,0.4,0]};

[X,labels]=generate_samples_decision_boundaries(simoptions);

labels1=isosplit2(X);
labels2=local_kmeans_sorber(X,3);

figure;
set(gcf,'Position',[50,50,1600,500]);
subplot_opts.yoffset=-0.02;

subplot_jfm(1,3,1,subplot_opts);
plot_samples_2d(X,labels);
title('Truth');

subplot_jfm(1,3,2,subplot_opts);
plot_samples_2d(X,labels1); hold on;
pt1=[1.1,4.0];
pt2=[-5.2,2.9]; pt2=pt2+(pt2-pt1)*10;
pt3=[2,5.1]; pt3=pt3+(pt3-pt1)*10;
pt1b=[1.7,3.6];
pt1c=[1.1,4.0];
pt4=[4.2,1.9]; pt4=pt4+(pt4-pt1b)*10;
plot([pt1c(1),pt2(1)],[pt1c(2),pt2(2)],'k--');
plot([pt1(1),pt3(1)],[pt1(2),pt3(2)],'k--');
plot([pt1b(1),pt4(1)],[pt1b(2),pt4(2)],'k--');
title('ISO-SPLIT');

subplot_jfm(1,3,3,subplot_opts);
plot_samples_2d(X,labels2); hold on; 
pt1=[0.8,2.1];
pt2=[-5.5,0.5]; pt2=pt2+(pt2-pt1)*10;
pt3=[1.0,5.7]; pt3=pt3+(pt3-pt1)*10;
pt4=[5.5,-1.2]; pt4=pt4+(pt4-pt1)*10;
plot([pt1(1),pt2(1)],[pt1(2),pt2(2)],'k--');
plot([pt1(1),pt3(1)],[pt1(2),pt3(2)],'k--');
plot([pt1(1),pt4(1)],[pt1(2),pt4(2)],'k--');
title('K-means');

set(gcf,'position',[0,0,14*100,4*100]);
set(gcf,'paperpositionmode','auto');
%set(gcf,'paperposition',[0,0,14,4]);
print([mfile_path,'/../images/decision_boundaries.eps'],'-depsc2');

end

function [samples,labels]=generate_samples_decision_boundaries(opts)

centers=opts.centers;
pops=opts.pops;
shapes=opts.shapes;

samples=zeros(2,0);
labels=zeros(1,0);

for j=1:length(centers)
	xx=randn(1,pops{j});
	yy=randn(1,pops{j});
	shape=shapes{j};
	xx2=xx*shape(1)+yy*shape(3);
	yy2=-xx*shape(3)+yy*shape(2);
	center=centers{j};
	xx2=xx2+center(1);
	yy2=yy2+center(2);
	tmp=zeros(2,pops{j});
	tmp(1,:)=xx2; tmp(2,:)=yy2;
	samples=[samples,tmp];
	labels=[labels,ones(1,pops{j})*j];
end;

end

function plot_samples_2d(samples,labels)

if (nargin<4) subtitle=''; end;

colors='rgbkymc';
for j=1:max(labels)
	xx=samples(1,find(labels==j));
	yy=samples(2,find(labels==j));
	col=colors(mod(j-1,length(colors))+1);
	plot(xx,yy,['.',col]);
	if (j==1) hold on; end;
end;
set(gca,'xtick',[],'ytick',[])

xmin=min(samples(1,:)); xmax=max(samples(1,:));
ymin=min(samples(2,:)); ymax=max(samples(2,:));
xlim([floor(xmin-2),ceil(xmax+2)]);
ylim([floor(ymin-1),ceil(ymax+2)]);

end

function [L,C]=local_kmeans_sorber(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
%        D = cumsum(sqrt(dot(D,D,1)));  % orig, seems to be dist (l=1)
        D = cumsum(dot(D,D,1));  % Arthur-Vassilvitskii use dist^2 (l=2)
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end

end
