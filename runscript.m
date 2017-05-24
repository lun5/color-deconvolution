% script to calculate the stain density using Ruifrok and Johnston
% for all normalized images by Tammy Ma
% code written by  Luong Nguyen July 25, 2016
%{
data_dir = 'D:\Documents\Tiles_Norm';
im_dir = fullfile(data_dir,'normalized_images');
im_dir = fullfile(data_dir,'Target15Norm');
stain_dir = fullfile(data_dir,'stain_densities');

im_dir = fullfile(data_dir,'Target15ReNorm');
stain_dir = fullfile(im_dir,'stain_densities');

method_names = {'Khan','Macenko','Reinhard','Vahadane', 'VahadaneFast','Luong'};

for mm = 1:length(method_names)
   curr_dir = fullfile(im_dir,method_names{mm});
   output_dir = fullfile(stain_dir,method_names{mm});
   if ~exist(output_dir,'dir')
       mkdir(output_dir);
   end
   imlist = dir(fullfile(curr_dir,'*.tif'));
   imlist = {imlist.name}';
   tic;
   parfor i = 1:length(imlist)
       %imname = 'klmqi6sq7wl6-mbdqhorkuxs';
       imname = imlist{i}(1:end-4);
       output_name = fullfile(output_dir,[imname '.mat']);
       if ~exist(output_name,'file')
           im_rgb = imread(fullfile(curr_dir,[imname '.tif']));
           [stain_density,~, ~] = calc_stain_density(im_rgb);
           parsave(output_name,stain_density);
       end
       %csvwrite(fullfile(output_dir,[imname '.csv']),stain_density);
   end
   T = toc;
   fprintf('Done with method %s in %.2f seconds\n',method_names{mm},T);
end
%}
%% Target
%{
% first test if the stain vector can actually separate the stains
target_list = {'4d0ylpdlwf','13nedzdzfh','9uixinhtjjis','anaggwovpxanwq0',...
    '95f7k8loesyevi','dj0ebjbuxshxz','jbakl4tseqt','k2yxq1tbr6kpny0',...
    'sfzywlg4291ch9','ocmmhhrtzz5','2ale5ngryfnpo',...
    'ffwtgxylhyna','ibhyyugefpbn','klmqi6sq7wl6','r4oq0dcxzlbdi'};
%target_list = dir(fullfile(data_dir,'Target','*.tif'));
%target_list = {target_list.name}';

for i = 1:length(target_list)
    %imname = 'klmqi6sq7wl6';
    imname = target_list{i}; %(1:end-4);
    %if ~exist(fullfile(data_dir,'Target',[imname '.mat']),'file')
        im_rgb = imread(fullfile(data_dir,'Target',[imname '.tif']));
        [stain_density, hem_im_rgb, eosin_im_rgb] = calc_stain_density(im_rgb,1);
        parsave(fullfile(data_dir,'Target',[imname '.mat']),stain_density);
    %end
    %figure; subplot(1,2,1); imshow(hem_im_rgb); subplot(1,2,2); imshow(eosin_im_rgb);
    % reconstruct the image    
end
%}
%% Original Images
%{
first test if the stain vector can actually separate the stains
target_dir = fullfile(data_dir,'Tiles_512');
target_list = dir(fullfile(target_dir,'*.tif'));
target_list = {target_list.name}';

for i = 1:length(target_list)
    imname = target_list{i}(1:end-4);
    %if ~exist(fullfile(data_dir,'Target',[imname '.mat']),'file')
        im_rgb = imread(fullfile(target_dir,[imname '.tif']));
        [stain_density, hem_im_rgb, eosin_im_rgb] = calc_stain_density(im_rgb,1);
        parsave(fullfile(target_dir,[imname '.mat']),stain_density);
    %end
    %figure; subplot(1,2,1); imshow(hem_im_rgb); subplot(1,2,2); imshow(eosin_im_rgb);
    % reconstruct the image    
end
%}
%%
% Compare between target and source
%{
im_dir = fullfile(data_dir,'Target15ReNorm');
stain_dir = fullfile(im_dir,'stain_densities');
target_dir = fullfile(data_dir,'Tiles_512');
if ~exist(fullfile(im_dir,'HistDist'),'dir')
    mkdir(fullfile(im_dir,'HistDist'));
end
PDollar_toolbox_path = 'C:\Users\luong_nguyen\Documents\GitHub\toolbox';
addpath(genpath(PDollar_toolbox_path));
for mm = 1:length(method_names)
   tic;
   curr_dir = fullfile(im_dir,method_names{mm});
   output_dir = fullfile(stain_dir,method_names{mm});
   imlist = dir(fullfile(curr_dir,'*.tif'));
   imlist = {imlist.name}';
   source_images = cell(length(imlist),1);
   target_images = cell(length(imlist),1);
   ks_dist_hem = cell(length(imlist),1);
   ks_dist_eos = cell(length(imlist),1);
   ks_dist = cell(length(imlist),1); 
   chisq_dist_hem = cell(length(imlist),1);
   chisq_dist_eos = cell(length(imlist),1);
   chisq_dist = cell(length(imlist),1);
   emd_dist_hem = cell(length(imlist),1);
   emd_dist_eos = cell(length(imlist),1);
   emd_dist =  cell(length(imlist),1);
   kl_dist_hem = cell(length(imlist),1);
   kl_dist_eos = cell(length(imlist),1);
   kl_dist =  cell(length(imlist),1);
   
   parfor i = 1:length(imlist)
      imname =  imlist{i}(1:end-4);
      imname_split = strsplit(imname,'-');
      source_images{i} = imname_split{2};
      target_images{i} = imname_split{1};
      
      tmp = load(fullfile(output_dir,[imname '.mat']));
      source_stain_densities = tmp.data; % tmp.stain_density;
      if sum(isnan(source_stain_densities(:)) > 0)
         fprintf('Method %s has image %s with NAN stain densities\n',...
             method_names{mm}, imname);
         chisq_dist_hem{i} = 15;
         chisq_dist_eos{i} = 15;
         emd_dist_hem{i} = 15;
         emd_dist_eos{i} = 15;
         ks_dist_hem{i} = 15;
         ks_dist_eos{i} = 15;
         continue;
      end
      
      [cs_hem,es_hem] = histcounts(source_stain_densities(1,:),100,'Normalization', 'probability');
      [cs_eos,~] = histcounts(source_stain_densities(2,:),'BinEdges',es_hem,'Normalization', 'probability');
      
      %tmp = load(fullfile(data_dir,'Target',[target_images{i} '.mat']),'stain_density');
      tmp = load(fullfile(target_dir,[source_images{i} '.mat']));
      target_stain_densities = tmp.data;
      [ct_hem,~] = histcounts(target_stain_densities(1,:),'BinEdges',es_hem,'Normalization', 'probability');
      [ct_eos,~] = histcounts(target_stain_densities(2,:),'BinEdges',es_hem,'Normalization', 'probability');
      
      [~,~,ks_dist_hem{i}] = kstest2(source_stain_densities(1,:),target_stain_densities(1,:));
      [~,~,ks_dist_eos{i}] = kstest2(source_stain_densities(2,:),target_stain_densities(2,:));
      ks_dist{i} = ks_dist_hem{i} + ks_dist_eos{i};
      
      chisq_dist_hem{i} = pdist2(cs_hem, ct_hem,'chisq');
      chisq_dist_eos{i} = pdist2(cs_eos,ct_eos,'chisq');
      chisq_dist{i} = chisq_dist_hem{i} + chisq_dist_eos{i};
      
      emd_dist_hem{i} = pdist2(cs_hem, ct_hem,'emd');
      emd_dist_eos{i} = pdist2(cs_eos,ct_eos,'emd');
      emd_dist{i} = emd_dist_hem{i} + emd_dist_eos{i};
      
      tmp = cs_hem.*log((cs_hem+eps)./(ct_hem+eps)); tmp(isnan(tmp)) = 0;
      kl_dist_hem{i} = sum(tmp,2);
      tmp = cs_eos.*log((cs_eos+eps)./(ct_eos+eps)); tmp(isnan(tmp)) = 0;
      kl_dist_eos{i} = sum(tmp,2);
      kl_dist{i} = kl_dist_hem{i} + kl_dist_eos{i};
      
      %chisq_dist_hem{i} = sum((cs_hem - ct_hem).^2./(cs_hem+ct_hem+eps))./2;
      %emd_dist_hem{i} = sum(abs(cumsum(cs_hem) - cumsum(ct_hem))); 
   end
   T = table(source_images,target_images,ks_dist_hem,ks_dist_eos, ks_dist, ...
       chisq_dist_hem,chisq_dist_eos,chisq_dist,...
       emd_dist_hem, emd_dist_eos, emd_dist, kl_dist_hem, kl_dist_eos,kl_dist);
   T.Properties.VariableNames = {'SourceImages','TargetImages',...
       'ks_hem','ks_oes','ks','chisq_hem','chisq_eos','chisq',...
       'emd_hem','emd_eos','emd','kl_hem','kl_eos','kl'};
   writetable(T,fullfile(im_dir,'HistDist',[method_names{mm} '.txt']),'Delimiter',',');
   timepassed = toc;
   fprintf('Done with method %s in %.2f seconds\n',method_names{mm},timepassed);
end
%}

%% box plot for HSV space
%{
Compare between target and source
if ~exist(fullfile(im_dir,'HistDist_HSV'),'dir')
    mkdir(fullfile(im_dir,'HistDist_HSV'));
end
PDollar_toolbox_path = 'C:\Users\luong_nguyen\Documents\GitHub\toolbox';
addpath(genpath(PDollar_toolbox_path));
metric_names = {'ks','chisq','emd','kl'};
%im_dir = fullfile(data_dir,'Target15Norm');
%target_dir = fullfile(data_dir,'Target');

im_dir = fullfile(data_dir,'Target15ReNorm');
target_dir = fullfile(data_dir,'Tiles_512');
for mm = 1:length(method_names)
   tic;
   curr_dir = fullfile(im_dir,method_names{mm});
   output_dir = fullfile(stain_dir,method_names{mm});
   imlist = dir(fullfile(curr_dir,'*.tif'));
   imlist = {imlist.name}';
   source_images = cell(length(imlist),1);
   target_images = cell(length(imlist),1);
   ks_dist = cell(length(imlist),1); 
   chisq_dist = cell(length(imlist),1);
   emd_dist =  cell(length(imlist),1);
   kl_dist =  cell(length(imlist),1);
   
   measures = cell(length(imlist),1);
      
   parfor i = 1:length(imlist)
      imname =  imlist{i}(1:end-4);
      imname_split = strsplit(imname,'-');
      source_images{i} = imname_split{2};
      target_images{i} = imname_split{1};
           
      source_im = imread(fullfile(curr_dir,[imname '.tif']));
      %target_im = imread(fullfile(target_dir,[target_images{i} '.tif']));
      target_im = imread(fullfile(target_dir,[source_images{i} '.tif']));
      source_hsv = rgb2hsv(source_im);
      target_hsv = rgb2hsv(target_im);
      
      metrics = zeros(1,length(metric_names)*(size(target_hsv,3)+1));
      channel_hist_counts = cell(2,3);
      for cc = 1:3
         channel = source_hsv(:,:,cc); 
         [channel_hist_counts{1,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
         channel = target_hsv(:,:,cc); 
         [channel_hist_counts{2,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
      end   
          
      for metr = 1:length(metric_names)
          metric = metric_names{metr};
          start_indx = (metr-1)*(size(target_hsv,3)+1) + 1;
          if strcmp(metric,'ks')
              for cc = 1:3
                ss = source_hsv(:,:,cc); tt = target_hsv(:,:,cc);
                [~,~,metrics(start_indx+cc)] = kstest2(ss(:),tt(:));
              end
          else
              for cc = 1:3
                 metrics(start_indx+cc) = pdist2(channel_hist_counts{1,cc}, ...
                     channel_hist_counts{2,cc},metric);
              end
          end
          metrics(start_indx) = sum(metrics((start_indx+1):(start_indx+3)));
      end
      measures{i} = metrics;
   end
   measures = num2cell(cat(1,measures{:}));
   T = cell2table(cat(2,source_images,target_images,measures));
   T.Properties.VariableNames = {'SourceImages','TargetImages',...
       'ks','ks_h','ks_s','ks_v','chisq','chisq_h','chisq_s','chisq_v',...
       'emd','emd_h','emd_s','emd_v','kl','kl_h','kl_s','kl_v'};
   writetable(T,fullfile(im_dir,'HistDist_HSV',[method_names{mm} '.txt']),'Delimiter',',');
   timepassed = toc;
   fprintf('Done with method %s in %.2f seconds\n',method_names{mm},timepassed);
end
%}
%% boxplot for the methods
%{
ks_scores = cell(length(method_names),1);
chisq_scores = cell(length(method_names),1);
emd_scores = cell(length(method_names),1);
kl_scores = cell(length(method_names),1); 
group_names = cell(length(method_names),1);
plot_names =  {'Khan','Macenko','Reinhard','Vahadane', 'VahadaneFast','HEColStats'};
for mm = 1:length(method_names)
   %T = readtable(fullfile(im_dir,'HistDist',[method_names{mm} '.txt']),'Delimiter',',');
   T = readtable(fullfile(im_dir,'HistDist_HSV',[method_names{mm} '.txt']),'Delimiter',',');
   ks_scores{mm} = T.ks;
   chisq_scores{mm} = T.chisq;
   emd_scores{mm} = T.emd;
   kl_scores{mm} = T.kl;
   gp = cell(length(ks_scores{mm}),1);
   gp(:) = {plot_names{mm}};
   group_names{mm} = gp;
end

t test
calculate the median
ks_means = cell(length(method_names),1);
emd_means = cell(length(method_names),1);
chisq_means = cell(length(method_names),1);
kl_means = cell(length(method_names),1);
for mm = 1:length(method_names)
%     ks_means{mm} = median(ks_scores{mm});%
%     emd_means{mm} = median(emd_scores{mm}); %
%     chisq_means{mm} = median(chisq_scores{mm});%
%     kl_means{mm} = median(kl_scores{mm});
    ks_means{mm} = mean(ks_scores{mm});
    emd_means{mm} = mean(emd_scores{mm});
    chisq_means{mm} = mean(chisq_scores{mm});
    kl_means{mm} = mean(kl_scores{mm});
end

metric_names = {'ks_scores','emd_scores','chisq_scores','kl_scores'};
mean_metrics = cell(length(metric_names),1);
mean_metrics{1} = cat(1,ks_means{:});
mean_metrics{2} = cat(1,emd_means{:});
mean_metrics{3} = cat(1,chisq_means{:});
mean_metrics{4} = cat(1,kl_means{:});

for met = 1:4
    [sort_metric,sort_id] = sort(mean_metrics{met},'ascend');
    fprintf('\n\nRanking for metric %s is \n',metric_names{met});
    for mm = 1:length(sort_id)
        %fprintf ('\t Method %s with median %.4f\n',method_names{sort_id(mm)},sort_metric(mm));
        fprintf ('\t Method %s with mean %.4f\n',method_names{sort_id(mm)},sort_metric(mm));
    end
    for mm = 1:(length(sort_id)-1)
       method1 = method_names(sort_id(mm));
       method2 = method_names(sort_id(mm+1));    
       [h,p] = eval(['ttest2(',metric_names{met},'{sort_id(mm)},',...
          metric_names{met},'{sort_id(mm+1)})']);
       fprintf('ttest p-value for metric %s of method %s (mean %.2f) and method %s (mean %.2f) is %.4f\n',...
          metric_names{met}, method1{1}, sort_metric(mm), method2{1}, sort_metric(mm+1), p);
%        p = eval(['signrank(',metric_names{met},'{sort_id(mm)},',...
%            metric_names{met},'{sort_id(mm+1)})']);
%        fprintf('sign rank test p-value for metric %s of method %s (median %.2f) and method %s (median %.2f) is %.4f\n',...
%            metric_names{met}, method1{1}, sort_metric(mm), method2{1}, sort_metric(mm+1), p);
    end
end

% box plot
ks_scores = cat(1,ks_scores{:});
chisq_scores = cat(1,chisq_scores{:});
emd_scores = cat(1,emd_scores{:});
kl_scores = cat(1,kl_scores{:});
group_names = cat(1,group_names{:});

figure; boxplot(chisq_scores,group_names);
set(gca,'FontSize',16);
ylabel('\chi^2 differences'); %xlabel('Methods');

figure; boxplot(emd_scores,group_names);
set(gca,'FontSize',16);
ylabel('Earth Mover distances'); %xlabel('Methods');

figure; boxplot(ks_scores,group_names);
set(gca,'FontSize',16);
ylabel('KS statistics'); %xlabel('Methods');
%

figure; boxplot(kl_scores,group_names);
set(gca,'FontSize',16);
ylabel('KL divergences'); 
%}
%% box plot for user ranking
human_ranking_dir = 'D:\Documents\Tiles_Norm\Human_ranking';
file_list = dir(fullfile(human_ranking_dir,'user_results*.txt'));
file_list = {file_list.name}';
num_sources = 10;
num_targets = 2;
method_names = {'Luong', 'Macenko', 'Reinhard','Khan', 'Vahadane', 'Vahadane_Fast'};

group_names = cell(6,1);
for mm = 1:6
    gp = cell(num_targets*num_sources,1);
    gp(:) = {method_names{mm}};
    group_names{mm} = gp;
end
group_names = cat(1,group_names{:});

for i = 1: length(file_list)
   fname = fullfile(human_ranking_dir,file_list{i});
   T = readtable(fname);
   % test the logic
   
   summary_table = zeros(num_targets*num_sources,2+length(method_names));
   % col 1: source, col 2: target, col 3 ->8: method  1 ==> 6   
   count = 0;
   for tt = 1:num_targets
       for ss = 1:num_sources
           count = count + 1;
           summary_table(count,1) = ss;
           summary_table(count,2) = tt;
           indx_pair = (T.source_num == ss) & (T.target_num == tt);
           %scores = zeros(length(method_names),1);
           
           for mm = 1:length(method_names)
               %scores(mm) = sum(pairwise_comp_mat(:,mm) == 1) + 0.5*sum(pairwise_comp_mat(:,mm) == 0);
               summary_table(count,mm+2) = sum(T.m1_num(indx_pair) == mm & T.Results(indx_pair) == 1) + ...
                   sum(T.m2_num(indx_pair) == mm & T.Results(indx_pair) == 0) + ...
                   0.5*sum((T.m1_num(indx_pair) == mm |T.m2_num(indx_pair) == mm) & T.Results(indx_pair) == -1);
           end
       end
   end
   
   mean_scores = mean(summary_table(:,3:8),1); 
   med_scores = median(summary_table(:,3:8),1); 
   sum_scores = sum(summary_table(:,3:8), 1);
   [sort_scores, sort_id] = sort(mean_scores,'descend');
   [~, sort_id_med] = sort(med_scores,'descend');
   
   count_inconst = 0;
   %[sort_scores, sorted_id] = sort(med_scores,'descend');
   for tt = 1:num_targets
       for ss = 1:num_sources
           indx_pair = (T.source_num == ss) & (T.target_num == tt);
           %pairwise_comp_mat = sparse([T.m1_num(indx_pair) T.m2_num(indx_pair)],...
           %    [T.m2_num(indx_pair) T.m1_num(indx_pair)], [T.Results(indx_pair) T.Results(indx_pair)]);
           pairwise_comp_mat = sparse(T.m1_num(indx_pair),...
               T.m2_num(indx_pair), T.Results(indx_pair),length(method_names),length(method_names));
           pairwise_comp_mat = full(pairwise_comp_mat);
           for m1 = 1:5
               for m2 = (m1+1):6
                   m1_num = sort_id(m1);
                   m2_num = sort_id(m2);
                   if m1_num > max(size(pairwise_comp_mat)) || m2_num > max(size(pairwise_comp_mat))
                       %disp('something wrong');
                       %uiwait;
                       break;
                   end
                   comp_result = pairwise_comp_mat(m1_num,m2_num);
                   %if (comp_result == 1 && sort_scores(m1_num) <= sort_scores(m2_num))...
                   %        || (comp_result == 0 && sort_scores(m1_num) >= sort_scores(m2_num))...
                   %        || (comp_result == -1 && sort_scores(m1_num) ~= sort_scores(m2_num))
                       %fprintf('method %d score is %.2f, method %d score is %.2f\n', ...
                       %    m1_num, sorted_scores(m1_num),m2_num,  sorted_scores(m2_num));
                   if (comp_result == 0)% || (comp_result == -1 && mean_scores(m1_num)>mean_scores(m2_num))
                       count_inconst = count_inconst + 1;
                   end
               end
           end
       end
   end
   
   avg_inconst = count_inconst/300;
   fprintf('\n\nUser %s ',upper(file_list{i}(14:end-4)));
   fprintf('%.2f %0.2f %d  %0.2f %d %.2f\n', sum_scores(1), mean_scores(1),...
       sort_id(1), med_scores(1), sort_id_med(1),avg_inconst);
   
   for mm = 1:(length(sort_id)-1)
       method1 = method_names(sort_id(mm));
       method2 = method_names(sort_id(mm+1));
       [~,p] = ttest2(summary_table(:,sort_id(mm) +2), summary_table(:,sort_id(mm+1) +2));
       fprintf('ttest p-value for user %s of method %s (mean %.2f) and method %s (mean %.2f) is %.4f\n',...
           file_list{i}(14:end-4), method1{1}, sort_scores(mm), method2{1}, sort_scores(mm+1), p);
       p = signrank(summary_table(:,sort_id(mm) + 2), summary_table(:,sort_id(mm+1) + 2));
       fprintf('sign rank test p-value for  user %s of method %s (median %.2f) and method %s (median %.2f) is %.4f\n',...
           file_list{i}(14:end-4), method1{1}, median(summary_table(:,sort_id(mm) + 2)),...
           method2{1}, median(summary_table(:,sort_id(mm+1) + 2)), p);
   end
   
   scores = cat(1,summary_table(:,3),summary_table(:,4),...
       summary_table(:,5),summary_table(:,6),summary_table(:,7),summary_table(:,8));
   %figure; boxplot(scores,group_names);
   %title(upper(file_list{i}(14:end-4)));   
end

%{
% %% generate latex table
% data_dir = '/Users/lun5/Box Sync/ColorNormalizationPaper/Tiles_512_Validation_Data';
% %eval_dir = fullfile(data_dir,'Norm_Images','HistDist_HSV');
% eval_dir = fullfile(data_dir,'Norm_Images','HistDist');
% %eval_dir = fullfile(data_dir,'Renorm_Images','HistDist_HSV');
% %eval_dir = fullfile(data_dir,'Renorm_Images','HistDist');
% 
% method_names = {'Macenko','Reinhard','Khan','Vahadane', 'VahadaneFast','Luong'};
% metrics_gp = cell(length(method_names), 1);
% %ss_names = cell(length(method_names), 1);   
% %tt_names = cell(length(method_names), 1);   
% group_names = cell(length(method_names),1);
% for mm = 1:length(method_names)
%    T = readtable(fullfile(eval_dir,[method_names{mm} '.txt']),'Delimiter',',');
% %    T.Properties.VariableNames = {'SourceImages','TargetImages',...
% %        'ks','ks_h','ks_s','ks_v','chisq','chisq_h','chisq_s','chisq_v',...
% %        'emd','emd_h','emd_s','emd_v','kl','kl_h','kl_s','kl_v'};
% %    writetable(T,fullfile(eval_dir,[method_names{mm} '.txt']),'Delimiter',',');
%    %metrics_gp{mm} = table2array(T(:,3:4:end));
%    metrics_gp{mm} = table2array(T(:,3:3:end));
%    %ss_names{mm} = T.Source;
%    %tt_names{mm} = T.Target;
%    gp = cell(size(metrics_gp{mm},1),1);
%    gp(:) = {method_names{mm}};
%    group_names{mm} = gp;
% end
% 
% mean_metrics = cell(length(method_names),1);
% median_metrics = cell(length(method_names), 1);   
% 
% for mm = 1:length(method_names)
%    mean_metrics{mm} = mean(metrics_gp{mm});
%    median_metrics{mm} = median(metrics_gp{mm});
% end
% 
% mean_metrics = cat(1,mean_metrics{:});
% rank_means = zeros(size(mean_metrics));
% p_means = zeros(size(mean_metrics));
% median_metrics = cat(1,median_metrics{:});
% rank_medians = zeros(size(median_metrics));
% p_medians  = zeros(size(mean_metrics));
% metric_names = {'ks_scores','chisq_scores','emd_scores','kl_scores'};
% 
% for met = 1:length(metric_names)
%    [sort_metrics, sort_id] = sort(mean_metrics(:,met),'ascend');
%    [~,ii] = sort(sort_id);
%    rank_means(:,met) = ii;
%    fprintf('\n\nRanking for mean metric %s is \n',metric_names{met});
%    for mm = 1:length(sort_id)
%        fprintf('\t Method %s with mean %.4f\n',method_names{sort_id(mm)},sort_metrics(mm));
%    end
%    
%    for mm = 1:(length(sort_id) -1)
%        method1 = method_names(sort_id(mm));
%        method2 = method_names(sort_id(mm+1));
%        metrics_m1 = metrics_gp{sort_id(mm)}(:,met);
%        metrics_m2 = metrics_gp{sort_id(mm+1)}(:,met);    
%        [h,p] = ttest2(metrics_m1,metrics_m2);
%        p_means(sort_id(mm),met) = p;
%        %fprintf('ttest p-value for metric %s of method %s (mean %.2f) and method %s (mean %.2f) is %.4f\n',...
%        %     metric_names{met}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
%    end
%  
%    fprintf('\n\nRanking for median metric %s is \n',metric_names{met});
%    [sort_metrics, sort_id] = sort(median_metrics(:,met),'ascend');
%    [~,ii] = sort(sort_id);
%    rank_medians(:,met) = ii;
%    for mm = 1:length(sort_id)
%        fprintf('\t Method %s with median %.4f\n',method_names{sort_id(mm)},sort_metrics(mm)); 
%    end
%    
%    for mm = 1:(length(sort_id)-1)
%        method1 = method_names(sort_id(mm));
%        method2 = method_names(sort_id(mm+1));  
%        metrics_m1 = metrics_gp{sort_id(mm)}(:,met);
%        metrics_m2 = metrics_gp{sort_id(mm+1)}(:,met);    
%        p = signrank(metrics_m1,metrics_m2);
%        p_medians(sort_id(mm),met) = p;
%        %fprintf('sign rank test p-value for metric %s of method %s (median %.2f) and method %s (median %.2f) is %.4f\n',...
%        %   metric_names{met}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
%    end
% end
% 
% mean_score_rank_p = zeros(size(mean_metrics,1), 3*size(mean_metrics,2));
% med_score_rank_p = zeros(size(mean_metrics,1), 3*size(mean_metrics,2));
% for i = 1:length(metric_names)
%    mean_score_rank_p(:,(i-1)*3+1) = mean_metrics(:,i);
%    mean_score_rank_p(:,(i-1)*3+2) = rank_means(:,i);
%    mean_score_rank_p(:,(i-1)*3+3) = p_means(:,i);
%    med_score_rank_p(:,(i-1)*3+1) = median_metrics(:,i);
%    med_score_rank_p(:,(i-1)*3+2) = rank_medians(:,i);
%    med_score_rank_p(:,(i-1)*3+3) = p_medians(:,i);
% end
% 
% table_method_names = {'MK','RH','Khan','VH','VHF','SCAN'};
% for mm = 1:length(method_names)
%    fprintf('%s & %.3f & %d & %.2f & %.3f & %d & %.2f & %.3f & %d & %.2f & %.3f & %d & %.2f \\\\ \n',...
%        table_method_names{mm},mean_score_rank_p(mm,1),uint8(mean_score_rank_p(mm,2)),mean_score_rank_p(mm,3),...
%        mean_score_rank_p(mm,4),uint8(mean_score_rank_p(mm,5)),mean_score_rank_p(mm,6),...
%        mean_score_rank_p(mm,7),uint8(mean_score_rank_p(mm,8)),mean_score_rank_p(mm,9),...
%        mean_score_rank_p(mm,10),uint8(mean_score_rank_p(mm,11)),mean_score_rank_p(mm,12)); 
% end
% 
% fprintf('\n\n');
% for mm = 1:length(method_names)
%    fprintf('%s & %.3f & %d & %.2f & %.3f & %d & %.2f & %.3f & %d & %.2f & %.3f & %d & %.2f \\\\\n',...
%        table_method_names{mm},med_score_rank_p(mm,1),uint8(med_score_rank_p(mm,2)),med_score_rank_p(mm,3),...
%        med_score_rank_p(mm,4),uint8(med_score_rank_p(mm,5)),med_score_rank_p(mm,6),...
%        med_score_rank_p(mm,7),uint8(med_score_rank_p(mm,8)),med_score_rank_p(mm,9),...
%        med_score_rank_p(mm,10),uint8(med_score_rank_p(mm,11)),med_score_rank_p(mm,12)); 
% end
% 
%}


