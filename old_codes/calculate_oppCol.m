% calculate_oppcol space
% Inspired by Thesis by Elmaraghi 2003

%%
nstep = 1000; % for plotting 
options = struct('Normalize','on');
%options = struct('Normalize','off');
mu_s = 4; sigma_s = 2; % values for normalization

[ oppCol_image, indx_sat] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 

%% Plot the color opponency space (cite SIC)
%if strcmpi(plotflag,'on')
%% get the image name for plotting
split_string = regexp(imname,'\.','split');
savename = fullfile(resultdir,split_string{1});
%h1=plot(pink_stain_man_sic(1),pink_stain_man_sic(2),'bs','MarkerSize',10,...%'MarkerFaceColor',pink_manual_rgb./255,
%    'MarkerEdgeColor','k','LineWidth',3);
%h2=plot(purple_stain_man_sic(1),purple_stain_man_sic(2),'bo','MarkerSize',10,...%'MarkerFaceColor',purple_manual_rgb./255,
%    'MarkerEdgeColor','r','LineWidth',3);
%legend([h1, h2],'manual pink stain', 'manual purple stain');

X = oppCol_image';
opts = statset('Display','final');
%dist_type = 'city'; 
dist_type ='sqEuclidean';
[idx,ctrs] = kmeans(X,2,'Distance',dist_type,...
    'Replicates',5,'Options',opts, 'start','sample');

%%
h1 = figure; 
scatter(oppCol_image(1,1:nstep:end),oppCol_image(2,1:nstep:end),20,rgb_image(:,1:nstep:end)'./255,'filled');
%plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
%plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
hold off
axis([-1 1 -1 1])
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h1,'-dtiff', [savename '_oppCol_rep.tiff']);
end

%%
h2 = figure; 
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
hold off
axis([-1 1 -1 1])
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h2,'-dtiff', [savename '_oppCol_rep.tiff']);
end