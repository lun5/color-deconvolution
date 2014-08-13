% script for training 
% okay, first of all, find out how many images for 1 svs
addpath(pwd); close all;

workdir = '/Users/lun5/Research/color_deconvolution';
datadir = fullfile(workdir, 'TissueImages');

source_svs = 'tp10-867-1';
purple_source = load(fullfile(workdir,'results','140806',[source_svs 'training_purple.mat']));
training_data_purple_source =purple_source.training_data_purple;
pink_source = load(fullfile(workdir,'results','140806',[source_svs 'training_pink.mat']));
training_data_pink_source = pink_source.training_data_pink;
 %[ training_data_purple_source, training_data_pink_source] = wsi_get_training( workdir, source_imname);
target_svs = 'tp10-611';
%[ training_data_purple_target, training_data_pink_target] = wsi_get_training( workdir, target_imname);
purple_source = load(fullfile(workdir,'results','140806',[target_svs 'training_purple.mat']));
training_data_purple_target =purple_source.training_data_purple;
pink_source = load(fullfile(workdir,'results','140806',[target_svs 'training_pink.mat']));
training_data_pink_target = pink_source.training_data_pink;

%% get the rotation matrix 
% source image
training_data = [training_data_purple_source(:,1:2000) training_data_pink_source(:,1:8000)];
[U,~,~] = svd(training_data,0);
source_rotation_matrix = [-U(:,1) U(:,2:3)]';
% target image
training_data = [training_data_purple_target(:,1:2000) training_data_pink_target(:,1:8000)];
[U,~,~] = svd(training_data,0);
target_rotation_matrix = [-U(:,1) U(:,2:3)]';

%% Normalization or equalization
%target_imname = 'tp10-611_53248_12288_2048_2048.tif';
target_imname = 'tp10-611_22528_14336_2048_2048.tif';
target_im = imread(fullfile(datadir,target_imname));
source_imname = 'tp10-867-1_47104_22528_2048_2048.tif';
%source_imname = 'tp10-867-1_38912_20480_2048_2048.tif';
source_im = imread(fullfile(datadir,source_imname));

source_eq_image = oppColNormalization( source_im, target_im, source_rotation_matrix, target_rotation_matrix );

