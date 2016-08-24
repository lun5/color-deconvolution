%% Run JSEG algorithm on the data
% This actually behaves quite well
% Luong Nguyen 10/05/2015

github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
addpath(genpath(github_dir));
jseg_dir = fullfile(github_dir,'otherMethods','JSEG');
cd(jseg_dir);
%im_dir = 'Z:\HEproject\data\normalization_512_jpg\';
im_dir = 'Z:\HEproject\data\normalization_512_jpg_gray\';
%im_dir = 'Z:\Tiles_512_jpg\';
%output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','multi_scale');
%output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','one_scale');
output_dir = 'Z:\HEproject\normalized_evaluation_results\JSEG_grayscale';
%output_dir = fullfile('Z:\HEproject\normalized_evaluation_results\JSEG','new_params');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

if ~exist(fullfile(output_dir,'seg_im'),'dir')
    mkdir(fullfile(output_dir,'seg_im'));
end

if ~exist(fullfile(output_dir,'bdry_im'),'dir')
    mkdir(fullfile(output_dir,'bdry_im'));
end

if ~exist(fullfile(output_dir,'mat_files'),'dir')
    mkdir(fullfile(output_dir,'mat_files'));
end

if ~exist(fullfile(output_dir,'gif_files'),'dir')
    mkdir(fullfile(output_dir,'gif_files'));
end

qthresh_vec = [50 150 250 350 450 550 600];
merge_vec = [0.05 .1 .2 .3 .5 .7 .9 1];
scale_vec = [1 2 3];
% params = combvec(qthresh_vec, merge_vec, scale_vec);
load('Jsegparams.mat','params');
params_scale1 = params(:, params(3,:) == 1);
params_scale2 = params(:, params(3,:) == 2);
params_scale3 = params(:, params(3,:) == 3);
save('Jsegparams_scale1.mat','params_scale1');
save('Jsegparams_scale2.mat','params_scale2');
save('Jsegparams_scale3.mat','params_scale3');
%params = params(:,1:10:end);
%quantize_threshold = 25:25:600;
%num_thres = length(quantize_threshold);
num_thres = size(params,2);
im_list = dir(fullfile(im_dir,'*.jpg'));
im_list = {im_list.name}';
num_images = length(im_list);
quote = '''';
nrow = 512; ncol = 512;
parfor i = 1:num_images
    T = tic; 
    im_name = im_list{i}(1:end-4);
    fprintf('Start with image %s...',im_name);
    for j = 1: num_thres
        q_thresh = params(1,j); m_thresh = params(2,j); scale = params(3,j);
        gif_file = fullfile(output_dir,'gif_files',[im_name '_qthr'  num2str(q_thresh)...
            '_mthresh' num2str(m_thresh) '_scale' num2str(scale) '.gif']);
        if ~exist(gif_file,'file')
%            expr = ['segwin -i ', fullfile(im_dir,im_list{i}), ' -t 6 -r9 ', ...
%              gif_file,' -q ' num2str(q_thresh),' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
            expr = ['segwin -i ', fullfile(im_dir,im_list{i}), ' -t 3 -s ',...
                num2str(nrow),' ', num2str(ncol),' -r9 ', gif_file,...
                ' -q ' num2str(q_thresh),' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
            %expr = ['segwin -i ', fullfile(im_dir,im_list{i}), ' -t 6 -r9 ', ...
            %    gif_file,' -l 1 -q ' num2str(q_thresh)];
            %out_expr = evalc(['system(' quote expr quote ')']);
            out_expr = evalc_parfor(expr);
        end
    end
    t = toc(T); fprintf(' Done in %.2f seconds\n',t);
end

%% get segs by scales
for scale = 1:3
params_scale = params(:,params(3,:) == scale);
num_thres = size(params_scale,2);    
parfor i = 1:num_images
    T = tic; 
    im_name = im_list{i}(1:end-4);
    fprintf('Start with image %s...',im_name);
    I = imread( fullfile(im_dir,im_list{i}));
    segs = cell(num_thres,1);
    for j = 1:num_thres
        %q_thresh = quantize_threshold(j);
        %gif_file = fullfile(output_dir,'gif_files',[im_name '_'  num2str(q_thresh) '.gif']);
        q_thresh = params(1,j); m_thresh = params(2,j); %scale = params(3,j);
        gif_file = fullfile(output_dir,'gif_files',[im_name, '_qthr',  num2str(q_thresh),...
            '_mthresh', num2str(m_thresh), '_scale', num2str(scale), '.gif']);
        if exist(gif_file,'file')
            bdry_im_fname = fullfile(output_dir,'bdry_im',[im_name, '_' num2str(q_thresh),...
                '_mthresh', num2str(m_thresh), '_scale', num2str(scale), '_bdry.jpg']);
            seg_im_fname = fullfile(output_dir,'bdry_im',[im_name, '_' num2str(q_thresh),...
                '_mthresh', num2str(m_thresh), '_scale', num2str(scale), '_seg.jpg']);    
            labels = imread(gif_file,1);
            segs{j} = labels+1;
            if ~exist(bdry_im_fname,'file')
                edge_map = seg2bdry(labels,'imageSize');
                % change thickness of edges
                edge_map = imdilate(edge_map, strel('disk',1));
                edge_map_im = I.*uint8(repmat(~edge_map,[1 1 1]));
                imwrite(edge_map_im,bdry_im_fname);
                %imwrite(label2rgb(labels),seg_im_fname);
            end
        end
    end
    parsave(fullfile(output_dir,'mat_files',[im_name '_scale' num2str(scale) '.mat']),segs);
    t = toc(T); fprintf(' Done in %.2f seconds\n',t);
end
end

disp('done');
%color quantization threshold - specify  
%values 0-600, leave blank for automatic determination. 
%The higher the value, the less number of quantized colors in the image. 
%For color images, try 250. If you are unsatisfied with the result 
%because two neighboring regions with similar colors are not getting separated, 
%please try a smaller value say 150.