%% script to run Felz-Hutt EGB method
% Luong Nguyen 09/21/2015

%github_dir = '/home/lun5/github/HE-tumor-object-segmentation';
%addpath(genpath(github_dir));

code_dir =  '/home/lun5/ColorNorm/codes';
addpath(genpath(code_dir));
egb_dir = '/home/lun5/Documents/segment';
cd(egb_dir);

%im_dir = '/home/lun5/HEproject/data/Tiles_512_ppm';
im_dir = '/home/lun5/ColorNorm/Target15Norm';
%result_dir = '/home/lun5/HEproject/evaluation_results/EGB';
result_dir = '/home/lun5/ColorNorm/EGB_results';
seg_result_dir = fullfile(result_dir,'segmentations_seism');
matfile_result_dir = fullfile(result_dir,'segmented_images_seism');  

if ~exist(result_dir,'dir')
    mkdir(result_dir)
end

if ~exist(seg_result_dir,'dir')
    mkdir(seg_result_dir);
end

if ~exist(matfile_result_dir,'dir')
    mkdir(matfile_result_dir);
end

method_names = {'Luong','Khan','Macenko','Reinhard','Vahadane','VahadaneFast'};
%tmp = load('params_seism.mat'); params = tmp.params;

sig = 0.5; k = 500; min_val = 50;
for mm = 1:length(method_names)
   method = method_names{mm};
   curr_dir = fullfile(im_dir,method);

   curr_outdir = fullfile(seg_result_dir,method);
   if ~exist(curr_outdir,'dir'); mkdir(curr_outdir);end

   curr_matfile_dir = fullfile(matfile_result_dir,method);
   if ~exist(curr_matfile_dir,'dir'); mkdir(curr_matfile_dir); end;

   im_list = dir(fullfile(curr_dir,'*.ppm'));
   im_list = {im_list.name}';
   
   tt = tic;
   parfor i = 1:length(im_list)
      	imname = im_list{i}(1:end-4);
	segs = cell(1);
 	cmm = ['./segment ' num2str(sig) ' ' num2str(k) ' ' ...
		num2str(min_val) ' ' fullfile(curr_dir,[im_name '.ppm']) ' ' ...
		 fullfile(curr_outdir,[im_name '.ppm'])]; 
	s = evalc_parfor(cmm); 
        I = imread(fullfile(curr_outdir,[im_name '.ppm']));
	segs{1} = rgb2label(I);
        outFile = fullfile(curr_matfile_dir,[im_name,'.mat']);  
        parsave(outFile, segs);
   end
   fprintf('Done with method %s in %.2f seconds\n',method, toc(tt));
end

disp('Done');


