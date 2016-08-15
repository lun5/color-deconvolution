% color normalization manual
% Normalization script
% Luong Nguyen
% July 12, 2014

workdir = '/Users/lun5/Research/color_deconvolution'; 
datadir = fullfile(workdir, '20x_images','TIFF');
close all;
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));
if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

prompt = 'Please enter target image name ';
target_imname = input(prompt,'s');
while isempty(target_imname)
    target_imname = input(prompt,'s');
end

prompt = 'Please enter source image name ';
source_imname = input(prompt,'s');
while isempty(source_imname)
    source_imname = input(prompt,'s');
end

% get the names of the files
target_im = imread(fullfile(datadir,target_imname));
target_rgb = raw2rgb(target_im);figure;imshow(target_im);
source_im = imread(fullfile(datadir,source_imname));
source_rgb = raw2rgb(source_im);figure;imshow(source_im);

[ target_purple_manual_rgb,target_pink_manual_rgb, target_stain_mat_man,target_saturation_mat_man] = ...
    deconvolutionManual( target_imname,datadir,resultdir);% options ); 
% rescale the image in OD space then convert back to rgb space
% plot the intensity/saturation/amounts of purple
figure; hist(target_saturation_mat_man(1,:), 100);range(target_saturation_mat_man(1,:))
mean(target_saturation_mat_man(1,:))
% plot the intensity/saturation/amounts of purple
hist(target_saturation_mat_man(2,:), 100);range(target_saturation_mat_man(2,:))
mean(target_saturation_mat_man(2,:))

%% do the same for the source image
[ source_purple_manual_rgb,source_pink_manual_rgb, source_stain_mat_man,source_saturation_mat_man] = ...
    deconvolutionManual( source_imname,datadir,resultdir);% options ); 
% rescale the image in OD space then convert back to rgb space
% plot the intensity/saturation/amounts of purple
hist(source_saturation_mat_man(1,:), 100);range(source_saturation_mat_man(1,:))
mean(source_saturation_mat_man(1,:))
% plot the intensity/saturation/amounts of purple
figure; hist(source_saturation_mat_man(2,:), 100);range(source_saturation_mat_man(2,:))
mean(source_saturation_mat_man(2,:))

% calculate the OD values for source and target image
% stats(target) = func(stats(source))
% target_stats_purple = 
% target_stats_pink = 
% source_stats_purple = 
% source_stats_pink = 
% we can first try with mean
ratio_purple = mean(target_saturation_mat_man(1,target_saturation_mat_man(1,:)>0))...
    /mean(source_saturation_mat_man(1,source_saturation_mat_man(1,:)>0));
ratio_pink = mean(target_saturation_mat_man(2,target_saturation_mat_man(2,:)>0))...
    /mean(source_saturation_mat_man(2,source_saturation_mat_man(2,:)>0));

% convert source to target color space
source_in_target_space_od = target_stain_mat_man*source_saturation_mat_man;
source_in_target_space_rgb = od2rgb(source_in_target_space_od, size(source_im,1), size(source_im,2));
figure; imshow(source_in_target_space_rgb);

source_normalize_saturation = [source_saturation_mat_man(1,:)*ratio_purple;...
    source_saturation_mat_man(2,:)* ratio_pink];
source_normalize_od = target_stain_mat_man*source_normalize_saturation;
source_normalize_rgb = od2rgb(source_normalize_od, size(source_im,1), size(source_im,2));
figure; imshow(source_normalize_rgb)
% purple(source) => func(purple(source))
% pink(source) => func(pink(source))
% new source Sat matrix = [func(purple(source));func(pink(source))];
% normalize source OD = Sat * source stain matrix;
% reconstructed

