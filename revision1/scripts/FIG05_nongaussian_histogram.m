function FIG05_nongaussian_histogram

Z=randn(1,5e7);
Z2=log(abs(Z+3));
figure; hist(Z2,10000);
xlim([-2.5,2.5]);
set(gcf,'Color','w');
set(gca,'xtick',[])
set(gca,'ytick',[])

set(gcf,'paperposition',[0,0,5,3.5]);
print('../images/nongaussian_histogram.eps','-depsc2');

end
