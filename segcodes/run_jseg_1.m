%% Run JSEG algorithm on the data
% This actually behaves quite well
% Luong Nguyen 10/05/2015

github_dir = 'C:\Users\luong_nguyen\Documents\GitHub\HE-tumor-object-segmentation';
addpath(genpath(github_dir));
jseg_dir = fullfile(github_dir,'otherMethods','JSEG');
cd(jseg_dir);

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


im_dir = 'D:\Documents\TilesNorm\Target15Norm';
result_dir = 'D:\Documents\TilesNorm\JSEG_results'
seg_dir = fullfile(result_dir,'seg_im');
mat_dir = fullfile(result_dir,'mat_files');
gif_dir = fullfile(result_dir,'gif_files');
bdry_dir = fullfile(result_dir,'bdry_files'); 

if ~exist(result_dir,'dir')
    mkdir(result_dir);
end

if ~exist(seg_dir)
    mkdir(seg_dir);
end

if ~exist(mat_dir)
    mkdir(mat_dir);
end

if ~exist(git_dir)
    mkdir(git_dir)
end

if ~exist(bdry_dir)
    mkdir(bdry_dir);
end

method_names = {'Luong','Khan','Macenko','Reinhard','Vahadane','VahadaneFast'};
q_thresh = 250;
m_thresh = 0.4;
scale = 1;
nrow = 512; ncol = 512;

for mm = 1:length(method_names)
   method = method_names{mm};
   curr_dir = fullfile(im_dir,method);

   curr_outdir = fullfile(seg_dir,method);
   if ~exist(curr_outdir,'dir'); mkdir(curr_outdir);end

   curr_matfile_dir = fullfile(mat_dir,method);
   if ~exist(curr_matfile_dir,'dir'); mkdir(curr_matfile_dir); end;

   curr_gif_dir = fullfile(gif_dir,method);
   if ~exist(curr_gif_dir,'dir'); mkdir(curr_gif_dir); end;
 
   im_list = dir(fullfile(curr_dir,'*.jpg'));
   im_list = {im_list.name}';  
   num_images = length(im_list); 
  
   tt = tic;
   parfor i = 1:num_images 
      imname = im_list{i}(1:end-4);
      gif_file = fullfile(curr_gif_dir, [im_name '.gif']);
      expr = ['segwin -i ', fullfile(curr_dir,im_list{i}), ' -t 3 -s ',...
	num2str(nrow),' ', num2str(ncol),' -r9 ', gif_file,... 
	' -q ' num2str(q_thresh),' -m ', num2str(m_thresh) ' -l ' num2str(scale)];
      if ~exist(gif_file,'file')
	s = evalc_parfor(expr);
      end
   end
   parfor i = 1:num_images
      gif_file = fullfile(curr_gif_dir, [im_name '.gif']);
      if exist(gif_file,'file')
      	labels = imread(gif_file,1);
      	segs{1} = labels + 1;
      end
      parsave(fullfile(curr_matfile_dir,[im_name '.mat']),segs);
   end
   t = toc(tt); fprintf('Done in %.2f seconds\n',t);
end
