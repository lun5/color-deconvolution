%% Run JSEG algorithm on the data
% This actually behaves quite well
% Luong Nguyen 10/05/2015

github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
addpath(genpath(github_dir));
jseg_dir = fullfile(github_dir,'otherMethods','JSEG');
cd(jseg_dir);

%% scan parameters
% load('Jsegparams.mat','params');
% params_scale1 = params(:, params(3,:) == 1);
% params_scale2 = params(:, params(3,:) == 2);
% params_scale3 = params(:, params(3,:) == 3);
% save('Jsegparams_scale1.mat','params_scale1');
% save('Jsegparams_scale2.mat','params_scale2');
% save('Jsegparams_scale3.mat','params_scale3');
% 
% %num_thres = size(params,2);
% num_thres = size(params_scale1,2);
% 
% im_dir = 'D:\Documents\Tiles_Norm\Target15Norm_jpg';
% result_dir = 'D:\Documents\Tiles_Norm\JSEG_results';
% mat_dir = fullfile(result_dir,'mat_files');
% gif_dir = fullfile(result_dir,'gif_files');
% bdry_dir = fullfile(result_dir,'bdry_files'); 
% 
% if ~exist(result_dir,'dir');mkdir(result_dir);end
% if ~exist(mat_dir,'dir'); mkdir(mat_dir); end
% if ~exist(gif_dir,'dir'); mkdir(gif_dir);end
% if ~exist(bdry_dir,'dir'); mkdir(bdry_dir);end
% 
% method = 'Luong';%method_names{mm};
% curr_dir = fullfile(im_dir,method);
% 
% curr_matfile_dir = fullfile(mat_dir,method);
% if ~exist(curr_matfile_dir,'dir'); mkdir(curr_matfile_dir); end;
% 
% curr_gif_dir = fullfile(gif_dir,method);
% if ~exist(curr_gif_dir,'dir'); mkdir(curr_gif_dir); end;
% 
% curr_bdry_dir = fullfile(bdry_dir,method);
% if ~exist(curr_bdry_dir,'dir'); mkdir(curr_bdry_dir); end;
% 
% im_list = dir(fullfile(curr_dir,'ocmmhhrtzz5-*.jpg'));
% im_list = {im_list.name}';
% num_images = length(im_list);
% nrow = 512; ncol = 512;
% parfor i = 1:num_images
%     T = tic; 
%     im_name = im_list{i}(1:end-4);
%     fprintf('Start with image %s...',im_name);
%     gif_file = fullfile(curr_gif_dir,[im_name '.gif']);
%     expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 6 -r9 ',gif_file, ...
%              ' -s ', num2str(nrow),' ', num2str(ncol),...
%              ' -m 0.3',' -l 3' ];
%     out_expr = evalc_parfor(expr);
% %     for j = 1: num_thres
% %         %q_thresh = params(1,j); m_thresh = params(2,j); scale = params(3,j);
% %         q_thresh = params_scale1(1,j); m_thresh = params_scale1(2,j); scale = params_scale1(3,j);
% %         gif_file = fullfile(curr_gif_dir,[im_name '_qthr'  num2str(q_thresh)...
% %             '_mthresh' num2str(m_thresh) '_scale' num2str(scale) '.gif']);
% %         if ~exist(gif_file,'file')
% %            % for color image
% %            expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 6 -r9 ',gif_file, ...
% %              ' -s ', num2str(nrow),' ', num2str(ncol),' -q ' num2str(q_thresh),...
% %              ' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
% %            % for gray image
% %            %expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 3 -s ',...
% %            %     num2str(nrow),' ', num2str(ncol),' -r9 ', gif_file,...
% %            %     ' -q ' num2str(q_thresh),' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
% %             out_expr = evalc_parfor(expr);
% %         end
% %     end
%     t = toc(T); fprintf(' Done in %.2f seconds\n',t);
% end
% 
% 
% parfor i = 1:num_images    
%     T = tic;
%     im_name = im_list{i}(1:end-4);
%     fprintf('Start with image %s...',im_name);
%     I = imread( fullfile(curr_dir,im_list{i}));
%     gif_file = fullfile(curr_gif_dir,[im_name '.gif']);
%     labels = imread(gif_file,1);
%     edge_map = seg2bdry(labels,'imageSize');
%     edge_map = imdilate(edge_map, strel('disk',1));
%     edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
%     bdry_im_fname = fullfile(curr_bdry_dir,[im_name, '_bdry.jpg']);
%     imwrite(edge_map_im,bdry_im_fname);
%     t = toc(T); fprintf(' Done in %.2f seconds\n',t);
% end

% parfor i = 1:8%num_images    
%     T = tic;
%     im_name = im_list{i}(1:end-4);
%     fprintf('Start with image %s...',im_name);
%     I = imread( fullfile(curr_dir,im_list{i}));
%     segs = cell(num_thres,1);
%     for j = 1:num_thres
%         %q_thresh = quantize_threshold(j);
%         %gif_file = fullfile(output_dir,'gif_files',[im_name '_'  num2str(q_thresh) '.gif']);
%         q_thresh = params_scale1(1,j); m_thresh = params_scale1(2,j); scale = params_scale1(3,j);
%         gif_file = fullfile(curr_gif_dir,[im_name '_qthr'  num2str(q_thresh)...
%             '_mthresh' num2str(m_thresh) '_scale' num2str(scale) '.gif']);
%         if exist(gif_file,'file')
%             bdry_im_fname = fullfile(curr_bdry_dir,[im_name, '_' num2str(q_thresh),...
%                 '_mthresh', num2str(m_thresh), '_scale', num2str(scale), '_bdry.jpg']);
%             labels = imread(gif_file,1);
%             segs{j} = labels+1;
%             if ~exist(bdry_im_fname,'file')
%                 edge_map = seg2bdry(labels,'imageSize');
%                 % change thickness of edges
%                 edge_map = imdilate(edge_map, strel('disk',1));
%                 edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
%                 imwrite(edge_map_im,bdry_im_fname);
%             end
%         end
%     end
%     parsave(fullfile(curr_matfile_dir,[im_name '_scale' num2str(scale) '.mat']),segs);
%     t = toc(T); fprintf(' Done in %.2f seconds\n',t);
% end

%% check the print out of parameter and output
% matfile_dir = 'D:\Documents\Tiles_Norm\JSEG_results\segmented_images';
% bdry_dir = 'D:\Documents\Tiles_Norm\JSEG_results\test_bdry';
% if ~exist(bdry_dir,'dir'); mkdir(bdry_dir); end;
% 
% im_list = dir(fullfile(matfile_dir,'*.mat'));
% im_list = {im_list.name}';
% 
% parfor i =1:length(im_list)
%    im_name = im_list{i}(1:end-4);
%    if exist(fullfile(curr_dir,['ocmmhhrtzz5-' im_name '.jpg']),'file')
%        I = imread( fullfile(curr_dir,['ocmmhhrtzz5-' im_name '.jpg']));
%        tmp = load(fullfile(matfile_dir,[im_name '.mat']));
%        segs = tmp.data;
%        for j = 1:length(segs)
%            q_thresh = params_scale1(1,j); m_thresh = params_scale1(2,j); scale = 1;
%            bdry_im_name = fullfile(bdry_dir,[im_name, '_' num2str(q_thresh),...
%                '_mthresh', num2str(m_thresh), '_scale', num2str(scale), '_bdry.jpg']);
%            if ~exist(bdry_im_name,'file')
%                labels = segs{j};
%                edge_map =  seg2bdry(labels,'imageSize');
%                edge_map = imdilate(edge_map, strel('disk',1));
%                edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));
%                imwrite(edge_map_im,bdry_im_name);
%            end
%        end
%    end
% end

%% convert images to jpg
%tiff_dir = '/home/lun5/ColorNorm/Target15Norm';
%jpg_dir = '/home/lun5/ColorNorm/Target15Norm_jpg';

%method_names = {'Khan','Luong','Macenko','Reinhard','Macenko','Vahadane','VahadaneFast'};

%if ~exist(jpg_dir,'dir'); mkdir(jpg_dir); end;

%{
for mm = 1:length(method_names)
    method = method_names{mm};
    cd(fullfile(tiff_dir,method));
    if (~exist(fullfile(tiff_dir,method,'*.jpg'),'file')) & ...
	(~exist(fullfile(jpg_dir,method,'*.jpg'),'file'))
    fprintf('mogrify method %s\n',method);
    system('mogrify -format jpg *.tif'); 
    end
end

for mm = 1:length(method_names)
   method = method_names{mm};
   fprintf('method %s\n',method);
   if ~exist(fullfile(jpg_dir,method),'dir'); mkdir(fullfile(jpg_dir,method)); end
   %if exist(fullfile(tiff_dir,method,'*.jpg'),'file')
	fprintf('move file\n');
   	movefile(fullfile(tiff_dir,method,'*.jpg'),fullfile(jpg_dir,method));
   %end
end
%}
%im_dir = 'Z:\HEproject\data\normalization_512_jpg\';
%im_dir = 'Z:\HEproject\data\normalization_512_jpg_gray\';
%im_dir = 'Z:\Tiles_512_jpg\';
%output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','multi_scale');
%output_dir = fullfile('Z:\HEproject\evaluation_results\JSEG','one_scale');
%output_dir = 'Z:\HEproject\normalized_evaluation_results\JSEG_grayscale';
%output_dir = fullfile('Z:\HEproject\normalized_evaluation_results\JSEG','new_params');


%% run JSEG
im_dir = 'D:\Documents\Tiles_Norm\Target15Norm_jpg';
result_dir = 'D:\Documents\Tiles_Norm\JSEG_results';
bdry_dir = fullfile(result_dir,'bdry_files'); 
mat_dir = fullfile(result_dir,'mat_files');
gif_dir = fullfile(result_dir,'gif_files');

if ~exist(result_dir,'dir');mkdir(result_dir);end
if ~exist(mat_dir,'dir'); mkdir(mat_dir); end
if ~exist(gif_dir,'dir'); mkdir(gif_dir);end
if ~exist(bdry_dir,'dir'); mkdir(bdry_dir);end

%seg_dir = fullfile(result_dir,'seg_im');
%if ~exist(seg_dir,'dir');mkdir(seg_dir);end
method_names = {'Luong','Khan','Macenko','Reinhard','Vahadane','VahadaneFast'};
% q_thresh = 250;m_thresh = 0.4;scale = 1;
nrow = 512; ncol = 512;

for mm = length(method_names)
   method = method_names{mm};
   curr_dir = fullfile(im_dir,method);

   curr_matfile_dir = fullfile(mat_dir,method);
   if ~exist(curr_matfile_dir,'dir'); mkdir(curr_matfile_dir); end;

   curr_gif_dir = fullfile(gif_dir,method);
   if ~exist(curr_gif_dir,'dir'); mkdir(curr_gif_dir); end;
 
   curr_bdry_dir = fullfile(bdry_dir,method);
   if ~exist(curr_bdry_dir,'dir'); mkdir(curr_bdry_dir); end;

   im_list = dir(fullfile(curr_dir,'*.jpg'));
   im_list = {im_list.name}';  
   num_images = length(im_list); 
  
   tt = tic;
   parfor i = 1:num_images 
      im_name = im_list{i}(1:end-4);
      gif_file = fullfile(curr_gif_dir, [im_name '.gif']);
      expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 6 -r9 ',gif_file, ...
             ' -s ', num2str(nrow),' ', num2str(ncol),' -m 0.3',' -l 3' ];
%       % for color image
%       expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 6 -r9 ',gif_file, ...
%           ' -s ', num2str(nrow),' ', num2str(ncol),' -q ' num2str(q_thresh),...
%           ' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
%       % for gray image
%       expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 3 -s ',...
%           num2str(nrow),' ', num2str(ncol),' -r9 ', gif_file,... 
%           ' -q ' num2str(q_thresh),' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
      if ~exist(gif_file,'file')
          s = evalc_parfor(expr);
      end
   end
   parfor i = 1:num_images
      segs = cell(1,1);
      im_name = im_list{i}(1:end-4);
      %fprintf('imname is %s\n',im_name);
      try
        I =  imread( fullfile(curr_dir,[im_name '.jpg']));
      catch
        warning(sprintf('image %s cannot be read for method %s\n',im_name, method));
        I = zeros(nrow,ncol,3);
      end
      gif_file = fullfile(curr_gif_dir, [im_name '.gif']);
      bdry_im_fname = fullfile(curr_bdry_dir,[im_name, '_bdry.jpg']);
      if ~exist(gif_file,'file')
          fprintf('gif file does not exist for file %s with method %s\n',im_name,method); 
      end
      if exist(gif_file,'file') && (~exist(bdry_im_fname,'file'))
      	labels = imread(gif_file,1);
      	segs{1} = labels + 1;
        edge_map = seg2bdry(labels,'imageSize');
        edge_map = imdilate(edge_map, strel('disk',1));
        edge_map_im = I.*uint8(repmat(~edge_map,[1 1 3]));        
        imwrite(edge_map_im,bdry_im_fname);
        parsave(fullfile(curr_matfile_dir,[im_name '.mat']),segs);
      end
   end
   t = toc(tt); fprintf('Done in %.2f seconds\n',t);
end
