% script to test effect of different color spaces on segmentation
% Luong Nguyen
% July 7, 2014

% test the cv classification
% Options are: RGB, oppCol, HSV, Lab
close all;
num_folds = 10;
color_space = 'Lab';%'HSV'; %'RGB'; %'oppCol';
Acc = logReg_CV(num_folds, color_space, options);

