%% Script for evaluating segmentation results from color segmentation
% Luong Nguyen 7/27/16
% segmentation results computed using interNucDist methods
% segmentation results calculated by Burak Tosun
% color normalization run by Tammy Ma

%function script_eval_seg_col_norm()
 method_names = {'Luong','Khan','Macenko','Reinhard','Vahadane', 'VahadaneFast','non_normalized'};
 data_dir = '/home/lun5/ColorNorm';
% %im_dir = fullfile(data_dir,'Tiles_Norm');
% %seg_dir = fullfile(data_dir,'seg_results_15Norm');
% %seg_dir = fullfile(data_dir,'seg_results_15Norm_aug19');
 seg_dir = fullfile(data_dir,'JSEG_results','mat_files');
% %seg_dir = fullfile(data_dir,'EGB_results','segmented_images_seism');
% 
 seism_dir = '/home/lun5/github/seism';
 addpath(genpath(seism_dir));
% % List of measures to compute
 measures = {%
     'fb'  ,... % Precision-recall for boundaries
     'fop' ,... % Precision-recall for objects and parts
     'fr'  ,... % Precision-recall for regions
    'voi' ,... % Variation of information
     'nvoi',... % Normalized variation of information
     'pri' ,... % Probabilistic Rand index
     'sc'  ,'ssc' ,... % Segmentation covering (two directions)
     'dhd' ,'sdhd',... % Directional Hamming distance (two directions)
     'bgm' ,... % Bipartite graph matching
    'vd'  ,... % Van Dongen
    'bce' ,... % Bidirectional consistency error
    'gce' ,... % Global consistency error
    'lce' ,... % Local consistency error
    };

 % table of evaluation results are saved here
 eval_dir = fullfile(data_dir,'JSEG_results');%seg_dir; %fullfile(data_dir,'EGB_eval_results');
 if ~exist(eval_dir,'dir'); mkdir(eval_dir); end
 gt_set = 'all_files';
 gt_dir = fullfile('/home/lun5/HEproject/groundTruth/coarse_fine_GT_512_512/',gt_set);
 gt_list = dir(fullfile(gt_dir,'*.mat'));
 gt_list = {gt_list.name}';

%{
for mm = 1:length(method_names)
    tic;
    curr_seg_dir = fullfile(seg_dir,method_names{mm});
    seg_list = dir(fullfile(curr_seg_dir,'*.mat'));
    seg_list = {seg_list.name}';
    num_im = length(seg_list);
    fprintf('%s has %d images\n',curr_seg_dir,num_im);
    %for tt = 1:length(measures)
    %    eval([measures{tt} '_results = cell(num_im,1);']);
    %end
    
    %eval_results = zeros(num_im,length(measures));
    eval_results = cell(num_im,1);
    source_names = cell(num_im,1);
    target_names = cell(num_im,1);
    alpha = 0.75; 
    for ii = 1:num_im
        seg_name = seg_list{ii};
        split_name = strsplit(seg_name,'-');
        %target_names{ii} = split_name{1};
        %split_name = strsplit(split_name{2},'_');
        target_names{ii} = 'NA';
        source_names{ii} = split_name{1}(1:end-4);
        
        tmp = load(fullfile(curr_seg_dir,seg_name));
        seg = uint16(tmp.data{1});
        tmp = load(fullfile(gt_dir,[source_names{ii} '.mat']));
        gt = tmp.groundTruth{1}.Segmentation;
	groundTruth = tmp.groundTruth;
        %seg = gt; %tmp.seg;
	eval_im_results = zeros(1,length(measures)+1);
	if ~isempty(seg)
          for tt = 1:length(measures)
            %eval([measures{tt} '_results{ii} = eval_segm(seg,gt,measures{tt});']);
            %eval_results(ii,tt) = eval_segm(seg,gt,measures{tt});
	    result = eval_segm(seg,gt,measures{tt});
            eval_im_results(tt) = result(1); 
          end
	  [ff_score, bb_score] = evalRegions(groundTruth,seg);
	  if bb_score == -1
	    eval_im_results(tt+1) = ff_score;
	  else
	    eval_im_results(tt+1) = alpha*ff_score + (1-alpha)*bb_score;
	  end
	end
	eval_results{ii} = eval_im_results;
    end
    
    result_table = table(source_names,target_names);
    eval_results = num2cell(cat(1,eval_results{:}));
    result_table = [result_table eval_results]; 
    result_table.Properties.VariableNames = cat(2,{'Source','Target'},measures,{'f_overlap'});
    %result_table.Properties.VariableNames = cat(2,{'Source','Target'},{'regions'});
    fprintf('done with method %s in %.2f seconds\n',method_names{mm},toc);
    writetable(result_table,fullfile(eval_dir,[method_names{mm} '_' gt_set '.txt']),'Delimiter',',');
end


avg_metrics = zeros(length(method_names), length(measures)+1);

for mm = 1:length(method_names)
    T = readtable(fullfile(eval_dir,[method_names{mm} '_all_files.txt']),'Delimiter',','); 
    %for jj = 1:length(measures)+1   
    avg_metrics(mm,:) = mean(table2array(T(:,3:end))); 
    %end
end

% calculate the scores for well defined set
avg_metrics_wd =  zeros(length(method_names), length(measures)+1);
for mm = 1:length(method_names)
    T = readtable(fullfile(eval_dir,[method_names{mm} '_all_files.txt']),'Delimiter',',');
    source_names = T.Source;
    gt_dir_wd = '/home/lun5/HEproject/groundTruth/coarse_fine_GT_512_512/well_defined';
    gt_files_wd = dir(fullfile(gt_dir_wd,'*.mat'));
    gt_files_wd = {gt_files_wd.name}';
    gt_files_wd = cellfun(@(x) x(1:end-4), gt_files_wd, 'UniformOutput',false);
    [Lia, Lob] = ismember(source_names, gt_files_wd);
    T_wd = T(Lia,:);
    avg_metrics_wd(mm,:) = mean(table2array(T_wd(:,3:end)));
end
%}
% for the target ocm
%data_dir = '/Users/lun5/Box Sync/ColorNormalizationPaper/Tiles_512_Validation_Data';
eval_dir = fullfile(data_dir,'JSEG_results');
method_names = {'Macenko','Reinhard','Khan','Vahadane', 'VahadaneFast','Luong','non_normalized'};
metrics_gp = cell(length(method_names), 1);
ss_names = cell(length(method_names), 1);   
tt_names = cell(length(method_names), 1);   
group_names = cell(length(method_names),1);
for mm = 1:length(method_names)
   T = readtable(fullfile(eval_dir,[method_names{mm} '_all_files.txt']),'Delimiter',',');
   %indx = ismember(T.Target,{'jbakl4tseqt'}); 
   indx = 1:1:length(T.Source); 
   metrics_gp{mm} = table2array(T(indx,3:end));
   ss_names{mm} = T.Source;
   tt_names{mm} = T.Target;
   gp = cell(size(metrics_gp{mm},1),1);
   gp(:) = {method_names{mm}};
   group_names{mm} = gp;
end

mean_metrics = cell(length(method_names),1);
median_metrics = cell(length(method_names), 1);   

for mm = 1:length(method_names)
   mean_metrics{mm} = mean(metrics_gp{mm});
   median_metrics{mm} = median(metrics_gp{mm});
end

%metric_names = cat(2,measures,{'f_overlap'});
metric_names = {'fb','f_overlap'};
mean_metrics = cat(1,mean_metrics{:});
median_metrics = cat(1,median_metrics{:});

big_mean_metrics = mean_metrics;
big_med_metrics = median_metrics;
mean_metrics = mean_metrics(1:end-1,[1, 16]);
median_metrics = median_metrics(1:end-1,[1, 16]);

rank_means = zeros(size(mean_metrics));
p_means = zeros(size(mean_metrics));
rank_medians = zeros(size(median_metrics));
p_medians  = zeros(size(mean_metrics));

for met = 1:length(metric_names)
   [sort_metrics, sort_id] = sort(mean_metrics(:,met),'descend');
   [~,ii] = sort(sort_id);
   rank_means(:,met) = ii;

   fprintf('\n\nRanking for mean metric %s is \n',metric_names{met});
   for mm = 1:length(sort_id)
       fprintf('\t Method %s with mean %.4f\n',method_names{sort_id(mm)},sort_metrics(mm));
   end
   
   for mm = 1:(length(sort_id) -1)
       method1 = method_names(sort_id(mm));
       method2 = method_names(sort_id(mm+1));
       metrics_m1 = metrics_gp{sort_id(mm)}(:,met);
       metrics_m2 = metrics_gp{sort_id(mm+1)}(:,met);    
       [h,p] = ttest2(metrics_m1,metrics_m2);
       p_means(sort_id(mm),met) = p;
       fprintf('ttest p-value for metric %s of method %s (mean %.2f) and method %s (mean %.2f) is %.4f\n',...
            metric_names{met}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
   end
 
   fprintf('\n\nRanking for mean metric %s is \n',metric_names{met});
   [sort_metrics, sort_id] = sort(median_metrics(:,met),'descend');
   [~,ii] = sort(sort_id);
   rank_medians(:,met) = ii;
   for mm = 1:length(sort_id)
       fprintf('\t Method %s with median %.4f\n',method_names{sort_id(mm)},sort_metrics(mm)); 
   end
   
   for mm = 1:(length(sort_id)-1)
       method1 = method_names(sort_id(mm));
       method2 = method_names(sort_id(mm+1));  
       metrics_m1 = metrics_gp{sort_id(mm)}(:,met);
       metrics_m2 = metrics_gp{sort_id(mm+1)}(:,met);    
       p = signrank(metrics_m1,metrics_m2);
       p_medians(sort_id(mm),met) = p;
       fprintf('sign rank test p-value for metric %s of method %s (median %.2f) and method %s (median %.2f) is %.4f\n',...
          metric_names{met}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
   end
end

mean_score_rank_p = zeros(size(mean_metrics,1), 3*size(mean_metrics,2));
med_score_rank_p = zeros(size(mean_metrics,1), 3*size(mean_metrics,2));
for i = 1:length(metric_names)
   mean_score_rank_p(:,(i-1)*3+1) = mean_metrics(:,i);
   mean_score_rank_p(:,(i-1)*3+2) = rank_means(:,i);
   mean_score_rank_p(:,(i-1)*3+3) = p_means(:,i);
   med_score_rank_p(:,(i-1)*3+1) = median_metrics(:,i);
   med_score_rank_p(:,(i-1)*3+2) = rank_medians(:,i);
   med_score_rank_p(:,(i-1)*3+3) = p_medians(:,i);
end

table_method_names = {'MK','RH','Khan','VH','VHF','SCAN'};
fprintf('MEAN\n');
for mm = 1:length(method_names)-1
   fprintf('%s & %.4f & %d & %.2f & %.4f & %d & %.2f \\\\ \n',...
       table_method_names{mm},mean_score_rank_p(mm,1),uint8(mean_score_rank_p(mm,2)),mean_score_rank_p(mm,3),...
       mean_score_rank_p(mm,4),uint8(mean_score_rank_p(mm,5)),mean_score_rank_p(mm,6)); 
end

fprintf('\n\nMEDIAN\n');
for mm = 1:length(method_names)-1
   fprintf('%s & %.4f & %d & %.2f & %.4f & %d & %.2f \\\\\n',...
       table_method_names{mm},med_score_rank_p(mm,1),uint8(med_score_rank_p(mm,2)),med_score_rank_p(mm,3),...
       med_score_rank_p(mm,4),uint8(med_score_rank_p(mm,5)),med_score_rank_p(mm,6)); 
end

% fb_results = cell(num_im,1);
%    fop_results = cell(num_im,1);
%    fr_results = cell(num_im,1);
%    voi_results = cell(num_im,1);
%    nvoi_results= cell(num_im,1);
%    pri_results = cell(num_im,1);
%    sc_results = cell(num_im,1);
%    ssc_results = cell(num_im,1);
%    dhd_results = cell(num_im,1);
%    sdhd_results = cell(num_im,1);
%    bgm_results = cell(num_im,1);
%    vd_results = cell(num_im,1);
%    bce_results = cell(num_im,1);
%    gce_results = cell(num_im,1);
%    lce_results = cell(num_im,1);
