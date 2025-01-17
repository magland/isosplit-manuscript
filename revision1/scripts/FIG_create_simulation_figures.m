function FIG_create_simulation_figures

close all;

hgload('../images/simulation1.fig');
set_props;
set(gcf,'paperposition',[0,0,9,6]);
print('../images/simulation1.eps','-depsc2');

hgload('../images/simulation2.fig');
set_props;
set(gcf,'paperposition',[0,0,9,6]);
print('../images/simulation2.eps','-depsc2');

hgload('../images/simulation3.fig');
set_props;
set(gcf,'paperposition',[0,0,9,6]);
print('../images/simulation3.eps','-depsc2');

hgload('../images/simulation4.fig');
set_props;
set(gcf,'paperposition',[0,0,9,6]);
print('../images/simulation4.eps','-depsc2');

hgload('../simulation5.fig');
set_props;
set(gcf,'paperposition',[0,0,9,6]);
print('../images/simulation5.eps','-depsc2');

end

function set_props

marker_size=5;

ch1=get(gcf,'children');
for j=1:length(ch1)
	ch2=get(ch1(j),'children');
	for k=1:length(ch2)
		set(ch2(k),'markersize',marker_size);
	end;
end;

end
