function prepare_figure_for_isosplit_paper

close all;

clip_size=200;
labels_to_view=[4,6,8,10,11];

mfile_path=fileparts(mfilename('fullpath'));
%path0=[mfile_path,'/sort_dl12_20151208_NNF_r1_tet16_17/output_tet16'];
path0=['/home/magland/dev/ms_franklab/experiments/2016_03_15/sort_dl12_20151208_NNF_r1_tet16_17/output_tet16'];
% firings=readmda([path0,'/firings.mda']);
% outlier_scores=firings(5,:);
% inds=find(outlier_scores<=3);
% firings=firings(:,inds);
% K=max(firings(3,:));
% labels_map=zeros(1,K);
% labels_map(labels_to_view)=1:length(labels_to_view);
% labels=firings(3,:);
% labels=labels_map(labels);
% firings(3,:)=labels;
% inds=find(labels>0);
% firings=firings(:,inds);
% times=firings(2,:);
% labels=firings(3,:);
% 
% writemda64(firings,[mfile_path,'/firings_for_isosplit_paper.mda']);

%Let's read this particular firings so we know the label ordering is the
%same
firings=readmda([mfile_path,'/firings_for_isosplit_paper.mda']);
times=firings(2,:);
labels=firings(3,:);

fprintf('Reading raw...\n');
pre0=readmda([path0,'/pre0.mda']);
pre2=readmda([path0,'/pre2.mda']);
fprintf('Extracting clips...\n');
clips0=ms_extract_clips2(pre0,times,clip_size);
clips2=ms_extract_clips2(pre2,times,clip_size);
fprintf('ms_templates...\n');
templates=ms_templates(clips0,labels);

aa=randsample(length(times),1000);
FF=ms_event_features(clips2,3);

figure;
subplot(1,2,1);
ms_view_clusters(FF,labels);
xlabel(''); ylabel(''); zlabel('');
set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[]);
subplot(1,2,2);
ms_view_templates(templates);
set(gcf,'position',[450,1300,1300,450]);

end