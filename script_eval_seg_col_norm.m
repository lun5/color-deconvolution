%% Script for evaluating segmentation results from color segmentation
% Luong Nguyen 7/27/16
% segmentation results computed using interNucDist methods
% segmentation results calculated by Burak Tosun
% color normalization run by Tammy Ma

function script_eval_seg_col_norm()
method_names = {'Khan','Luong','Macenko','Reinhard','Vahadane', 'Vahadane Fast'};
data_dir = '/home/lun5/ColorNorm';
%im_dir = fullfile(data_dir,'Tiles_Norm');
seg_dir = fullfile(data_dir,'Tiles_seg_results');
gt_dir = '/home/lun5/HEproject/groundTruth/groundTruth_512_512_fine_coarse';

% List of measures to compute
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
eval_dir = fullfile('data_dir','eval_results');
if ~exist(eval_dir,'dir'); mkdir(eval_dir); end

for mm = 1%:length(method_names)
    tic;
    curr_seg_dir = fullfile(seg_dir,method_names{mm});
    seg_list = dir(fullfile(curr_seg_dir,'*.mat'));
    seg_list = {seg_list.name}';
    num_im = length(seg_list);
    
    %for tt = 1:length(measures)
    %    eval([measures{tt} '_results = cell(num_im,1);']);
    %end
    
    eval_results = zeros(num_im,length(measures));
    
    source_names = cell(num_im,1);
    target_names = cell(num_im,1);
    
    parfor ii = 1:num_im
        seg_name = seg_list{ii};
        split_name = strsplit(seg_name,'-');
        target_names{ii} = split_name{1};
        split_name = strsplit(split_name,'_');
        source_names{ii} = split_name{1};
        
        %tmp = load(fullfile(seg_dir,seg_name));
        %seg = tmp.seg;
        tmp = load(fullfile(gt_dir,[source_names{ii} '.mat']));
        gt = tmp.groundTruth{1}.Segmentation;
        seg = tmp.seg;
        for tt = 1:length(measures)
            %eval([measures{tt} '_results{ii} = eval_segm(seg,gt,measures{tt});']);
            eval_results(ii,tt) = eval_segm(seg,gt,measures{tt});
        end
    end
    
    result_table = table(source_names,target_names);
    %for tt = 1:length(measures)
    %   eval(['result_table.' measures{tt} ' =' measures{tt} '_results;']);
    %end
    result_table(:,3:(2+length(measures))) = eval_results;
    result_table.Properties.VariableNames = cat(2,{'Source','Target'},measures);
    writetable(result_table,fullfile(eval_dir,[method_names{mm} '.mat']));
    fprintf('done with method %s in %.2f seconds\n',method_names{mm},toc);
end
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
