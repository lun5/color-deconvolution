%% run script for MITOS

%% STEP 1: HSV calculation
data_dir = 'M:\MITOS';
dir_list = {'training_aperio','testing_aperio','training_hamamatsu','testing_hamamatsu'};
method_names = {'Luong','Macenko','Khan','Reinhard','Vahadane','VahadaneFast'};
PDollar_toolbox_path = 'C:\Users\luong_nguyen\Documents\GitHub\toolbox';
addpath(genpath(PDollar_toolbox_path));
metric_names = {'ks','chisq','emd','kl'};

variable_names = {'SourceImages','TargetImages',...
                's2t_ks','s2t_ks_h','s2t_ks_s','s2t_ks_v',...
                's2t_chisq','s2t_chisq_h','s2t_chisq_s','s2t_chisq_v',...
                's2t_emd','s2t_emd_h','s2t_emd_s','s2t_emd_v',...
                's2t_kl','s2t_kl_h','s2t_kl_s','s2t_kl_v',...
                't2s_ks','t2s_ks_h','t2s_ks_s','t2s_ks_v',...
                't2s_chisq','t2s_chisq_h','t2s_chisq_s','t2s_chisq_v',...
                't2s_emd','t2s_emd_h','t2s_emd_s','t2s_emd_v',...
                't2s_kl','t2s_kl_h','t2s_kl_s','t2s_kl_v'};
% for mm = 1:length(method_names)
%     met = method_names{mm};    
%     t1 = tic;
%     fprintf('Start with method %s\n',met);
%     measures_over_traintest = cell(2,1);
%     source_over_traintest = cell(2,1);
%     target_over_traintest = cell(2,1);
%     for dd = 1:2
%         source_list = dir(fullfile(data_dir,dir_list{dd},'*.tar'));
%         t2 = tic;
%         if ~exist(fullfile(data_dir,dir_list{dd},'HistDist_HSV'),'dir')
%             mkdir(fullfile(data_dir,dir_list{dd},'HistDist_HSV'));
%         end
%         source_list = {source_list.name}';
%         measures_over_sourcelist = cell(length(source_list),1); 
%         source_over_sourcelist = cell(length(source_list),1);
%         target_over_sourcelist = cell(length(source_list),1);
%         for ss = 1:length(source_list)
%             source_name = source_list{ss}(1:end-4);
%             t3 = tic;
%             resol_levels = dir(fullfile(data_dir,dir_list{dd},source_name,'frames'));
%             resol_levels = {resol_levels.name}';
%             resol_levels(strcmp(resol_levels,'.')) = [];
%             resol_levels(strcmp(resol_levels,'..')) = [];
%             indx = strfind(resol_levels,'stain');
%             resol_levels(~cellfun(@isempty,indx)) = [];
%             if exist(fullfile(data_dir,dir_list{dd},'HistDist_HSV',...
%                     [method_names{mm} '_' source_name '.txt']),'file')
%                 fprintf('Already calculated for dir %s\n', source_name);
%                 T1 = readtable(fullfile(data_dir,dir_list{dd},'HistDist_HSV',...
%                     [method_names{mm} '_' source_name '.txt']), 'Delimiter',',');
%                 measures_over_resol = table2cell(T1(:,3:end));
%                 source_over_resol = table2cell(T1(:,1));
%                 target_over_resol = table2cell(T1(:,2));
%                 measures_over_sourcelist{ss} = measures_over_resol;
%                 source_over_sourcelist{ss} = source_over_resol;
%                 target_over_sourcelist{ss} = target_over_resol;
%                 continue;
%             end
%             measures_over_resol = cell(length(resol_levels),1);
%             source_over_resol = cell(length(resol_levels),1);
%             target_over_resol = cell(length(resol_levels),1);
% 
%             for rr = 1:length(resol_levels)
%                 if strcmp(resol_levels{rr},'.') || strcmp(resol_levels{rr},'..')...
%                         || ~isempty(strfind(resol_levels{rr},'stain'))
%                     continue;
%                 end
%                 t4 = tic;
%                 source_dir = fullfile(data_dir,dir_list{dd},source_name,...
%                     'frames',resol_levels{rr});
%                 target_dir = fullfile(data_dir,dir_list{dd+2},...
%                     ['H' source_name(2:end)],'frames',resol_levels{rr});
%                 s2t_dir = fullfile(data_dir, dir_list{dd},source_name,...
%                     [met '_normalized_frames'],resol_levels{rr});
%                 
%                 if ~exist(s2t_dir,'dir'); mkdir(s2t_dir); end;
%                 t2s_dir = fullfile(data_dir, dir_list{dd+2},...
%                     ['H' source_name(2:end)],[met '_normalized_frames'],resol_levels{rr});
%                 if ~exist(t2s_dir,'dir'); mkdir(t2s_dir); end;
%                 
%             
%                 imlist = dir(fullfile(source_dir,'*.tiff'));
%                 imlist = {imlist.name}';
%                 measures_over_imlist = cell(length(imlist),1);
%                 source_over_imlist = cell(length(imlist),1);
%                 target_over_imlist = cell(length(imlist),1);
%                 parfor ii = 1:length(imlist)
%                     imname = imlist{ii};
%                     source_im = imread(fullfile(source_dir,imname));
%                     target_im = imread(fullfile(target_dir,['H' imname(2:end)]));
%                     source_over_imlist{ii} = cellstr(fullfile(source_dir,imname));
%                     target_over_imlist{ii} = cellstr(fullfile(target_dir,['H' imname(2:end)]));
%                     
%                     s2t_im = imread(fullfile(s2t_dir,imname));
%                     s2t_hsv = rgb2hsv(s2t_im); % normalized hsv
%                     t2s_im = imread(fullfile(t2s_dir,['H' imname(2:end)]));
%                     t2s_hsv = rgb2hsv(t2s_im);
%                     source_hsv = rgb2hsv(source_im);
%                     target_hsv = rgb2hsv(target_im);
%                     metrics = zeros(2,length(metric_names)*(size(target_hsv,3)+1));% row 1: s2t, row 2: t2s
%                     channel_hist_counts = cell(4,3);
%                     for cc = 1:3
%                         channel = source_hsv(:,:,cc);
%                         [channel_hist_counts{1,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
%                         channel = target_hsv(:,:,cc);
%                         [channel_hist_counts{2,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
%                         channel = s2t_hsv(:,:,cc);
%                         [channel_hist_counts{3,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
%                         channel = t2s_hsv(:,:,cc);
%                         [channel_hist_counts{4,cc},~] =  histcounts(channel(:),100,'Normalization', 'probability');
%                     end
%                     
%                     for metr = 1:length(metric_names)
%                         metric = metric_names{metr};
%                         start_indx = (metr-1)*(size(target_hsv,3)+1) + 1;
%                         if strcmp(metric,'ks')
%                             for cc = 1:3
%                                 sc = source_hsv(:,:,cc); tc = target_hsv(:,:,cc);
%                                 s2t = s2t_hsv(:,:,cc); t2s = t2s_hsv(:,:,cc);
%                                 [~,~,metrics(1,start_indx + cc)] = kstest2(s2t(:), tc(:));
%                                 [~,~,metrics(2,start_indx + cc)] = kstest2(t2s(:), sc(:));
%                             end
%                         else
%                             for cc = 1:3
%                                 % target vs s2t
%                                 metrics(1,start_indx+cc) = pdist2(channel_hist_counts{2,cc}, ...
%                                     channel_hist_counts{3,cc},metric);
%                                 % source vs t2s
%                                 metrics(2,start_indx+cc) = pdist2(channel_hist_counts{1,cc}, ...
%                                     channel_hist_counts{4,cc},metric);
%                             end
%                         end
%                         metrics(1,start_indx) = sum(metrics(1,(start_indx+1):(start_indx+3)));
%                         metrics(2,start_indx) = sum(metrics(2,(start_indx+1):(start_indx+3)));
%                     end
%                     measures_over_imlist{ii} = cat(2,metrics(1,:), metrics(2,:));
%                 end
%                 
%                 measures_over_imlist = num2cell(cat(1,measures_over_imlist{:}));
%                 source_over_imlist = cat(1,source_over_imlist{:});
%                 target_over_imlist = cat(1,target_over_imlist{:});
%                 clear t2s_im; clear s2t_im;
%                 measures_over_resol{rr} = measures_over_imlist;
%                 source_over_resol{rr} = source_over_imlist;
%                 target_over_resol{rr} = target_over_imlist;
%                 tt4 = toc(t4);
%                 fprintf('Done with resolution %s in %.2f s\n',resol_levels{rr},tt4);
%             end
%             measures_over_resol = num2cell(cat(1,measures_over_resol{:}));
%             source_over_resol = cat(1,source_over_resol{:});
%             target_over_resol = cat(1,target_over_resol{:});
%             measures_over_sourcelist{ss} = measures_over_resol;
%             source_over_sourcelist{ss} = source_over_resol;
%             target_over_sourcelist{ss} = target_over_resol;
%   
%             T1 = cell2table(cat(2,source_over_resol,target_over_resol,measures_over_resol));            
%             T1.Properties.VariableNames = variable_names;
%             writetable(T1,fullfile(data_dir,dir_list{dd},'HistDist_HSV',...
%                 [method_names{mm} '_' source_name '.txt']),'Delimiter',',');
%             tt3 = toc(t3);
%             fprintf('Done with dir %s in %.2f seconds\n',source_name,tt3);
%         end
%         measures_over_sourcelist = cat(1,measures_over_sourcelist{:});
%         source_over_sourcelist = cat(1,source_over_sourcelist{:});
%         target_over_sourcelist = cat(1,target_over_sourcelist{:});
%         
%         measures_over_traintest{dd} = measures_over_sourcelist;
%         source_over_traintest{dd} = source_over_sourcelist;
%         target_over_traintest{dd} = target_over_sourcelist;
%         tt2 = toc(t2);
%         fprintf('Done with dir %s in %.2f seconds\n',dir_list{dd},tt2);
%     end
%     measures_over_traintest = cat(1,measures_over_traintest{:});
%     source_over_traintest = cat(1,source_over_traintest{:});
%     target_over_traintest = cat(1,target_over_traintest{:});
%     T = cell2table(cat(2,source_over_traintest,target_over_traintest,measures_over_traintest));
%     
%     T.Properties.VariableNames = variable_names;
%     writetable(T,fullfile(data_dir,dir_list{dd},'HistDist_HSV',[method_names{mm} '.txt']),'Delimiter',',');
%     tt1 = toc(t1);
%     fprintf('Done with method %s in %.2f seconds\n',met,tt1);
% end

%% Statistics
dd = 2;
metrics_gp = cell(length(method_names), 1);
group_names = cell(length(method_names),1);
indx = 1:2058;%randperm(2058,300);
%method_names = {'Macenko','Reinhard','Khan','Vahadane', 'VahadaneFast','Luong'};
method_names = {'Vahadane','VahadaneFast','Luong'};
for mm = 1:length(method_names)
   T = readtable(fullfile(data_dir,dir_list{dd},'HistDist_HSV',[method_names{mm} '.txt']),'Delimiter',',');
   %metrics = sortrows(table2array(T(:,3:4:end)))./3;
   %indx = 1:300;
   metrics = table2array(T(:,3:4:end));
   metrics_gp{mm} = metrics(indx,:);
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

mean_metrics = cat(1,mean_metrics{:});
rank_means = zeros(size(mean_metrics));
p_means = zeros(size(mean_metrics));
median_metrics = cat(1,median_metrics{:});
rank_medians = zeros(size(median_metrics));
p_medians  = zeros(size(mean_metrics));

for met = 1:length(metric_names)*2
   [sort_metrics, sort_id] = sort(mean_metrics(:,met),'ascend');
   [~,ii] = sort(sort_id);
   rank_means(:,met) = ii;
   if met == 4 || met == 8
       met_indx = 4; 
   else
       met_indx = mod(met,4);
   end
   fprintf('\n\nRanking for mean metric %s is \n',metric_names{met_indx});
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
            metric_names{met_indx}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
   end
 
   fprintf('\n\nRanking for median metric %s is \n',metric_names{met_indx});
   [sort_metrics, sort_id] = sort(median_metrics(:,met),'ascend');
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
          metric_names{met_indx}, method1{1}, sort_metrics(mm), method2{1}, sort_metrics(mm+1), p);
   end
end

for met = 1:length(metric_names)*2
    fprintf('& Mean & %.2f & %.2f & %.2f \\\\ \n',mean_metrics(1,met), mean_metrics(2,met), mean_metrics(3,met));
    fprintf('&& Rank & %d & %d & %d \\\\ \n', rank_means(1,met), rank_means(2,met), rank_means(3,met));
    fprintf('&& p-val & %.2f & %.2f & %.2f \\\\ \n',p_means(1,met),p_means(2,met), p_means(3,met));
end

fprintf('\n\n');

for met = 1:length(metric_names)*2
    fprintf('& Med & %.2f & %.2f & %.2f \\\\ \n',median_metrics(1,met), median_metrics(2,met), median_metrics(3,met));
    fprintf('&& Rank & %d & %d & %d \\\\ \n', rank_medians(1,met), rank_medians(2,met), rank_medians(3,met));
    fprintf('&& p-val & %.2f & %.2f & %.2f \\\\ \n',p_medians(1,met),p_medians(2,met), p_medians(3,met));
end





