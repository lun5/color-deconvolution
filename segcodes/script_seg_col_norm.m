%% Script for evaluating segmentation results from color segmentation
% Luong Nguyen 7/27/16
% segmentation results computed using interNucDist methods
% segmentation results calculated by Burak Tosun
% color normalization run by Tammy Ma

%function script_eval_seg_col_norm()
method_names = {'Khan','Luong','Macenko','Reinhard','Vahadane', 'VahadaneFast'};
data_dir = '/home/lun5/ColorNorm';
%im_dir = fullfile(data_dir,'Tiles_Norm');
im_dir = fullfile(data_dir,'Target15Norm');
seg_dir = fullfile(data_dir,'seg_results_15Norm_aug19');
mask_dir = fullfile(data_dir,'masks');
if ~exist(seg_dir,'dir'); mkdir(seg_dir); end
adj_dir = fullfile(data_dir,'Target15Norm_withAdjDelas');
param_str = '_se1_minNuc2_minStr3_minLum3';

for mm = 1:length(method_names)
    tic;
    curr_im_dir = fullfile(im_dir,method_names{mm});
    curr_adj_dir = fullfile(adj_dir,method_names{mm});
    curr_seg_dir =fullfile(seg_dir,method_names{mm});
    if ~exist(curr_seg_dir,'dir'); mkdir(curr_seg_dir); end
 
    im_list = dir(fullfile(curr_im_dir,'*.tif'));
    im_list = {im_list.name}';
    num_im = length(im_list);
    fprintf('number of images is %d\n',num_im);
    
    parfor ii = 1:num_im
        imname = im_list{ii}(1:end-4);
    	segs = object_proposal_all_types(curr_im_dir, curr_adj_dir, mask_dir, imname, param_str);
        parsave(fullfile(curr_seg_dir,[imname '.mat']),{segs});  
        %fprintf('Done with image %s in method %s\n',imname, method_names{mm});
    end
    
    fprintf('done with method %s in %.2f seconds\n',method_names{mm},toc);
end

