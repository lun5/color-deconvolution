% write code to summary the results for each of the methods:
% each image has N number of settings
% each setting --> score board
% save the score for each image 
% calculate the mean at each setting
% choose the setting that have the best values

%% Windows
%githubdir ='C:\Users\tma\Documents\MATLAB\HE-tumor-object-segmentation-master'; 
%cd(githubdir)
%seismdir = 'C:\Users\tma\Documents\MATLAB\seism-master';  addpath(genpath(seismdir));
%DATA_DIR = 'D:\ChakraLab\Normalization\'; %%idli
%GT_DIR = 'D:\ChakraLab\Normalization\ground_truth\all_files';
%IMG_DIR = 'C:\Users\tma\Documents\MATLAB\Target15Norm_newresults\Khan';
%githubdir = '/Users/lun5/Research/github/HE-tumor-object-segmentation'; addpath(genpath(githubdir));
%seismdir = '/Users/lun5/Research/github/seism';  addpath(genpath(seismdir));
%DATA_DIR = '/Users/lun5/Research/HE_Segmentation/';
%GT_DIR = fullfile(DATA_DIR,'groundTruth','coarse_fine_GT_512_512','all_files');
%IMG_DIR = '/Users/lun5/Research/HE_color_normalization/interNucDistSeg';


githubdir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation'; addpath(genpath(githubdir));
seismdir = 'C:\Users\luong_nguyen\Documents\GitHub\seism';  addpath(genpath(seismdir));
DATA_DIR = 'D:\Documents\HE_Segmentation\data';
GT_DIR = fullfile(DATA_DIR,'GroundTruth','coarse_fine_GT_512_512','all_files');
IMG_DIR = 'D:\Documents\Tiles_Norm\Target15Norm_newresults_march';
%addpath(genpath(githubdir));
addpath(genpath(seismdir));

%all_methods = {'Luong', 'Macenko', 'Reinhard', 'Vahadane', 'VahadaneFast'};
all_methods = {'Luong','Khan','Macenko', 'Reinhard','Vahadane', 'VahadaneFast','NoNorm'};
%all_methods = {'NoNorm'};

train_fname = fullfile(githubdir,'otherMethods','train_tiles.txt');
train_table = readtable(train_fname,'ReadVariableNames',false, 'Delimiter',',');
train_table.Properties.VariableNames = {'tile_names','wsi_names'};
im_list = train_table.tile_names;

measures = {
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
        
test_fname = fullfile(githubdir,'otherMethods','test_tiles.txt');
test_table = readtable(test_fname,'ReadVariableNames',false, 'Delimiter',',');
test_table.Properties.VariableNames = {'tile_names','wsi_names'};
test_im_list = test_table.tile_names;

RESULTS_DIR = cell(length(all_methods),1);

for i = 1:length(all_methods)
   RESULTS_DIR{i} = fullfile(IMG_DIR,all_methods{i});
end

maxDist = 0.02;
maxDist_vec = 0.02;
target_names = {'ocmmhhrtzz5'};
%target_names = {'ffwtgxylhyna','ocmmhhrtzz5'};


for i= 1:length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
   fprintf('Start with method %s\n',all_methods{i});
   if ~exist(param_scan_dir,'dir')
       mkdir(param_scan_dir)
   end
   tic;
   parfor j = 1:length(im_list) %%load sources
      % load the mat file
      im_name = lower(im_list{j});
      % load ground truth
      tmp = load(fullfile(GT_DIR,[im_name '.mat']));
      gt = tmp.groundTruth{1}.Segmentation;
      groundTruth = tmp.groundTruth;      
      
      for tt = 1:length(target_names)
          target_name = target_names{tt};
          if strcmp(target_name, im_name)
              continue;
          end
          
          if i < 7
              outname = fullfile(param_scan_dir,[target_name '-' im_name '.mat']);
          else
              outname = fullfile(param_scan_dir,[im_name '.mat']);
          end
          %if exist(fullfile(param_scan_dir,[target_name '-' im_name '.mat']),'file')
          if exist(outname,'file')
              continue;
          end
          
          if i < 7
              fname = fullfile(RESULTS_DIR{i},[target_name '-' im_name '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']);
          else
              fname = fullfile(RESULTS_DIR{i},[im_name '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']);
          end
          
          %tmp = load(fullfile(RESULTS_DIR{i},[target_name '-' im_name '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']));
          if ~exist(fname,'file')
             fprintf('File %s does not exist\n',im_name);
             continue;
          end
          
          tmp = load(fname);
          segs = tmp.segs;
          nSegments = length(segs); % segments 2:2:200
          param_scan_results = cell(nSegments,1);
          t1 = tic;
          for k = 1:nSegments
              seg = segs{k,1};
              if size(seg,1) > size(gt,1)
                  seg = seg(1:4:end,1:4:end);
              end
              if min(seg(:)) ~= 1
                  diff = min(seg(:)) - 1;
                  seg = seg - diff;
              end
              seg = uint16(seg);
              curr_results = zeros(1,length(measures)+1);
              for m = 1:length(measures)
                  result = eval_segm(seg, gt, measures{m},maxDist);
                  curr_results(m) = result(1);
              end
              [ff_score, bb_score] = evalRegions(groundTruth,seg);
              if bb_score == -1 % no background in gt
                  curr_results(end) = ff_score;
              else
                  alpha = 0.75;
                  curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
              end
              param_scan_results{k} = curr_results;
          end
          fprintf('Finish with image %s in %.2f seconds\n',im_name,toc(t1));
          param_scan_results = cat(1,param_scan_results{:});
          
          %save(fullfile(param_scan_dir,[target_name '-' im_name '.mat']),'param_scan_results');
          parsave(outname,param_scan_results);
      end
   end
   fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);
end

%% Put together the results across multiple maxDist
% don't need this step because there is only one maxDist
%{
for i = length(all_methods)
    sum_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_summary');
    if ~exist(sum_scan_dir,'dir')
        mkdir(sum_scan_dir);
    end
    for j = 1 :length(im_list)
        im_name = im_list{j};
        if exist(fullfile(sum_scan_dir,[im_name '.mat']),'file')
            continue;
        end
        nSegments = 40;
        param_scan_results = zeros(nSegments, length(measures) + 4);
        for mm = 1:length(maxDist_vec)
            maxDist = maxDist_vec(mm);
            param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan' '_' num2str(maxDist)]);
            tmp = load(fullfile(param_scan_dir,['0rpetmn8ks.mat']));
            if mm == 1
                param_scan_results(:,1) = tmp.param_scan_results(:,1);
                param_scan_results(:,5:end) = tmp.param_scan_results(:,2:end);
            else
               param_scan_results(:,mm) = tmp.param_scan_results(:,1);
            end
        end
        save(fullfile(sum_scan_dir,[target_name '-' im_name '.mat']),'param_scan_results');
    end
end
%}

all_test_measures = cell(length(all_methods),1);

for i = 1:length(all_methods)
    fprintf('Start with method %s\n',all_methods{i});
    %segmented_image_dir = fullfile(RESULTS_DIR{i},'segmented_images');
    param_scan_dir = fullfile(RESULTS_DIR{i},['param_scan_' num2str(maxDist)]);
    nSegments = 40;
    mean_measures_thres = zeros(nSegments,length(measures)+1);
    count = 0;
    count_none = 0;
    tic;
    for j = 1:length(im_list)
        im_name = lower(im_list{j});
        for k=1:length(target_names)
            target_name = target_names{k};
            if strcmp(target_name, im_name)
                continue;
            end
            %tmp = load(fullfile(param_scan_dir,[target_name '-' im_name '.mat']));
            if i < 7
                fname = fullfile(param_scan_dir,[target_name '-' im_name '.mat']);
            else
                fname = fullfile(param_scan_dir,[im_name '.mat']);
            end
            
            if ~exist(fname,'file')
                fprintf('There is no %s\n',im_name)
                continue
            end
            tmp = load(fname);
            if isfield(tmp,'data')
               data = tmp.data; 
            else
               data = tmp.param_scan_results;
            end
            if ~isempty(data)
                %mean_measures_thres = mean_measures_thres + tmp.param_scan_results;
                mean_measures_thres = mean_measures_thres + data;
            else
                fprintf('File %s has no results\n',im_name);
                count_none = count_none+1;
            end
            count = count+1;           
        end
    end
    % find the best thres
    fprintf('There are %d images with no resutls\n',count_none);
    mean_measures_thres = mean_measures_thres./count;
    [max_val, max_id] = max(mean_measures_thres,[],1);
    save(fullfile(param_scan_dir, 'average_measures.mat'),...
        'mean_measures_thres','max_val','max_id');
    disp(max_val)
    %test_measures = zeros(length(test_im_list),length(measures)+1);
    test_measures = cell(length(test_im_list),1);
    if exist(fullfile(param_scan_dir, 'test_measures.mat'),'file')
        tmp = load(fullfile(param_scan_dir, 'test_measures.mat'));
        all_test_measures{i} = tmp.test_measures;
        mean_test_measures = tmp.mean_test_measures;
        fprintf('Method %s\n',all_methods{i});
        disp(mean_test_measures)
    else
        parfor j = 1:length(test_im_list)
            %load the mat file
            im_name = lower(test_im_list{j});
            for k=1:length(target_names)
                target_name = target_names{k};
                if strcmp(target_name, im_name)
                    continue;
                end
                %if exist(fullfile(param_scan_dir,[target_name '-' im_name '.mat']),'file')
                if exist(fullfile(param_scan_dir,[im_name '.mat']),'file')
                    continue;
                end
                
                tmp = load(fullfile(GT_DIR,[im_name '.mat']));
                gt = tmp.groundTruth{1}.Segmentation;
                groundTruth = tmp.groundTruth;
                
                %tmp = load(fullfile(RESULTS_DIR{i},[target_name '-' im_name ...
                %    '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']));
                
                if i < 7
                   fname =  fullfile(RESULTS_DIR{i},[target_name '-' im_name ...
                    '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']);
                else
                   fname = fullfile(RESULTS_DIR{i},[im_name ...
                    '_se1_minNuc3_minStr5_minLum5_segResults.mat_segmentationTilesResult.mat']);
                end
                tmp = load(fname);
                segs = tmp.segs;
                seg = segs{max_id,1};
                
                if size(seg,1) > size(gt,1)
                    seg = seg(1:4:end,1:4:end);
                end
                if min(seg(:)) ~= 1
                    diff = min(seg(:)) - 1;
                    seg = seg - diff;
                end
                seg = uint16(seg);
                
                curr_results = zeros(1,length(measures)+1);
                for m = 1:length(measures)
                    result = eval_segm(seg, gt, measures{m},maxDist);
                    curr_results(m) = result(1);
                end
                [ff_score, bb_score] = evalRegions(groundTruth,seg);
                if bb_score == -1 % no background in gt
                    curr_results(end) = ff_score;
                else
                    alpha = 0.75;
                    curr_results(end) = alpha*ff_score + (1-alpha)*bb_score;
                end
                %test_measures(j,:) = curr_results;
                test_measures{j} = curr_results;
            end
        end
        test_measures = cat(1,test_measures{:});
        mean_test_measures = mean(test_measures,1);
        fprintf('Method %s\n',all_methods{i});
        disp(mean_test_measures)
        save(fullfile(param_scan_dir, 'test_measures.mat'),...
            'test_measures','mean_test_measures');
    end
    fprintf('Finish with method %s in %.2f seconds\n',all_methods{i},toc);    
end


for i = 1:2
    for j = (i+1):length(all_methods)
        fprintf('Compare method %s to method %s\n', all_methods{i}, all_methods{j});
        fprintf('ttest\n');
        fprintf('Fb: ');
        [h,p] = ttest2(all_test_measures{i}(:,1),all_test_measures{j}(:,1));
        fprintf('Reject = %d, p_val = %.5f\n',h, p);
        fprintf('Dice: ');
        [h,p] = ttest2(all_test_measures{i}(:,3),all_test_measures{j}(:,3));
        fprintf('Reject = %d, p_val = %.5f\n',h, p);
        fprintf('Fr: ');
        [h,p] = ttest2(all_test_measures{i}(:,end),all_test_measures{j}(:,end));
        fprintf('Reject = %d, p_val = %.5f\n',h, p);
        
%         fprintf('sign rank\n');
%         fprintf('Fb\n');
%         [p,h] = signrank(all_test_measures{i}(:,1),all_test_measures{j}(:,1));
%         fprintf('Reject = %d, p_val = %.5f\n',h, p);
%         fprintf('Dice\n');
%         [p,h] = signrank(all_test_measures{i}(:,3),all_test_measures{j}(:,3));
%         fprintf('Reject = %d, p_val = %.5f\n',h, p);
%         fprintf('Fr\n');
%         [p,h] = signrank(all_test_measures{i}(:,end),all_test_measures{j}(:,end));
%         fprintf('Reject = %d, p_val = %.5f\n',h, p);
    end
end

%{
for i = length(all_methods)
   param_scan_dir = fullfile(RESULTS_DIR{i},'param_scan_summary');
   fprintf('Start with method %s\n',all_methods{i});
   tmp = load(fullfile(param_scan_dir, 'average_measures.mat'));
   max_val = tmp.max_val; max_id = tmp.max_id;

   tmp = load(fullfile(param_scan_dir,'average_measures.mat'));
   mean_measures_thres = tmp.mean_measures_thres;
   
  for mm =1:length(maxDist_vec)
      maxDist = maxDist_vec(mm);
      curr_max_id = max_id(mm);
      tmp = load(fullfile(param_scan_dir, ['test_measures_' num2str(maxDist) '.mat']));
      mean_test_measures = tmp.mean_test_measures;
      fprintf('%.4f, %.4f, %.4f,%.4f, %.4f, %.4f,%.4f, %.4f,%.2f\n',...
            mean_measures_thres(curr_max_id,[mm,5:6,end]),...
            mean_test_measures(1:3),mean_test_measures(end),thres(curr_max_id));
   end
end
%}



