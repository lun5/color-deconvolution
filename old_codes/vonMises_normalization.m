% input the images
workdir = '/Users/lun5/Research/color_deconvolution'; 
datadir = fullfile(workdir, 'TissueImages');
target_imname = 'tp10-611_22528_14336_2048_2048.tif';
source_imname = 'tp10-867-1_47104_22528_2048_2048.tif';

% get the names of the files
target_im = imread(fullfile(datadir,target_imname));
target_rgb = raw2rgb(target_im);figure;imshow(target_im);
source_im = imread(fullfile(datadir,source_imname));
source_rgb = raw2rgb(source_im);figure;imshow(source_im);
% get the training data
color_space = 'oppCol';
[ training_data, labels, rotation_matrix ] = import_training_data( color_space );
% convert the image to oppCol space
options = struct('Normalize','off');
mu_s = 4; sigma_s = 2; % values for normalization
[ target_oppCol, target_first_component] = rgb2oppCol( target_rgb, mu_s, sigma_s, rotation_matrix, options); 
[ source_oppCol, source_first_component] = rgb2oppCol( source_rgb, mu_s, sigma_s, rotation_matrix, options); 

