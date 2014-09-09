%script for training 
% okay, first of all, find out how many images for 1 svs
addpath(pwd); close all;clear all; 
addpath(fullfile(matlabroot,'toolbox','nhist'))
addpath(fullfile(matlabroot,'toolbox','CircStat2012a'))

packagedir = '/Users/lun5/Research/github';
addpath(genpath(fullfile(packagedir,'matlabPyrTools'))); 
addpath(genpath(fullfile(packagedir,'textureSynth'))); 
addpath(genpath(fullfile(packagedir,'colorTextureSynth'))); 

workdir = '/Users/lun5/Research/color_deconvolution';
datadir = fullfile(workdir, 'TissueImages');

target_svs = 'tp10-611';
source_svs = 'tp10-867-1';
% %[ training_data_purple_source, training_data_pink_source] = wsi_get_training( workdir, source_svs);
% purple_source = load(fullfile(workdir,'results','140806',[source_svs 'training_purple.mat']));
% %purple_source = load(fullfile(workdir,'results','140908',[source_svs 'training_purple.mat']));
% training_data_purple_source =purple_source.training_data_purple;
% pink_source = load(fullfile(workdir,'results','140806',[source_svs 'training_pink.mat']));
% training_data_pink_source = pink_source.training_data_pink;
% %[ training_data_purple_target, training_data_pink_target] = wsi_get_training( workdir, target_svs);
% purple_target = load(fullfile(workdir,'results','140806',[target_svs 'training_purple.mat']));
% %purple_target = load(fullfile(workdir,'results','140908',[target_svs 'training_purple.mat']));
% training_data_purple_target =purple_target.training_data_purple;
% pink_target = load(fullfile(workdir,'results','140806',[target_svs 'training_pink.mat']));
% training_data_pink_target = pink_target.training_data_pink;
% 
% %% get the rotation matrix 
% % source image
% source_training_data = [training_data_purple_source(:,1:2000) training_data_pink_source(:,1:6000)];
% [U,~,~] = svd(source_training_data,0);
% source_rotation_matrix = [-U(:,1) U(:,2:3)]';
% % target image
% target_training_data = [training_data_purple_target(:,1:2000) training_data_pink_target(:,1:6000)];
% [U,~,~] = svd(target_training_data,0);
% target_rotation_matrix = [-U(:,1) U(:,2:3)]';
% 
% %% get the rotation matrix 
% % source image
% source_training_data = [training_data_purple_source(:,1:2000) training_data_pink_source(:,1:6000)];
% cov = (source_training_data - repmat(mean(source_training_data,2),1,size(source_training_data,2)))*...
%     (source_training_data - repmat(mean(source_training_data,2),1,size(source_training_data,2)))';
% [Us,Ss,~] = svd(cov,0);
% source_rotation_matrix = [-Us(:,1) Us(:,2:3)]';
% % target image
% target_training_data = [training_data_purple_target(:,1:2000) training_data_pink_target(:,1:6000)];
% cov = (target_training_data - repmat(mean(target_training_data,2),1,size(target_training_data,2)))*...
%     (target_training_data - repmat(mean(target_training_data,2),1,size(target_training_data,2)))';
% [Ut,St,~] = svd(cov,0);
% target_rotation_matrix = [-Ut(:,1) Ut(:,2:3)]';
%% Normalization or equalization
% R-G, Y-B
rotation_matrix = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3);...
    1/sqrt(6) 1/sqrt(6) -2/sqrt(6);...
    -1/sqrt(2) 1/sqrt(2) 0];

source_rotation_matrix = rotation_matrix;
target_rotation_matrix = rotation_matrix;

imnames = {'mws09-778a_20480_22528_2048_2048.tif','tp10-834-2_36864_12288_2048_2048.tif',...
    'tp10-762-1_28672_20480_2048_2048.tif','tp09-2031-1_40960_20480_2048_2048.tif',...
    'tp10-1105-1_8192_8192_2048_2048.tif','tp10-1546-1_8192_20480_2048_2048.tif',...
    'tp10-420-1_8192_18432_2048_2048.tif','tp09-1003-1_16384_24576_2048_2048.tif',...
    'tp10-354-1_8192_2048_2048_2048.tif'}; %,'tp09-16-1_26624_4096_2048_2048.tif'};

%target_imname = 'tp10-611_53248_12288_2048_2048.tif';
target_imname = 'tp10-611_22528_14336_2048_2048.tif';
target_im = imread(fullfile(datadir,target_imname));
source_imname = 'tp10-867-1_47104_22528_2048_2048.tif';
%source_imname = 'tp10-867-1_38912_20480_2048_2048.tif';
source_im = imread(fullfile(datadir,source_imname));

options = struct('matchMethod','stats');
source_eq_image = oppColNormalization( source_im, target_im, source_rotation_matrix, target_rotation_matrix, options );

%%
%% test the effect of white patches on normalization
% half_white_pink_imname = 'tp10-611_53248_18432_2048_2048.tif';
% target_im = imread(fullfile(datadir,half_white_pink_imname));
% source_imname = 'tp10-867-1_40960_18432_2048_2048.tif';
% %source_imname = 'tp10-867-1_38912_20480_2048_2048.tif';
% source_im = imread(fullfile(datadir,source_imname));
% 
% source_eq_image = oppColNormalization( source_im, target_im, source_rotation_matrix, target_rotation_matrix, options );
% %% get rotation matrices
% figure; imshow(target_im);
% rect = getrect; 
% target_training = imcrop(target_im,rect);
% target_training_rgb = raw2rgb(target_training);
% target_training_rgb(:,mean(target_training_rgb,1) > 230) = [];
% target_training_rgb(:,mean(target_training_rgb,1) < 20) = [];
% [U,~,~] = svd(target_training_rgb(:,randsample(size(target_training_rgb,2),10000)));
% target_rotation_matrix = [-U(:,1) U(:,2:3)]';
% 
% figure; imshow(source_im);
% rect = getrect; 
% source_training = imcrop(source_im,rect);
% source_training_rgb = raw2rgb(source_training);
% source_training_rgb(:,mean(source_training_rgb,1) > 230) = [];
% source_training_rgb(:,mean(source_training_rgb,1) < 20) = [];
% [U,~,~] = svd(source_training_rgb(:,randsample(size(source_training_rgb,2),10000)));
% source_rotation_matrix = [-U(:,1) U(:,2:3)]';
% 
% close all;