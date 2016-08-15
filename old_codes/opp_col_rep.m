% Examine images projected onto opponent color space. 
% also should know why we use svd

% directories
addpath(pwd);
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/color_deconvolution'; 
%datadir = uigetdir('*.tiff', 'Please select the image folder');
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));
close all;
if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('Read the input')
imagepaths = getImInfo(datadir,'TIF');
% load the training data
color_space = 'oppCol';
[ training_data, labels, rotation_matrix ] = import_training_data( color_space);
options = struct('Normalize','on');
mu_s = 4; sigma_s = 2; % values for normalization
% loop through the image and convert to oppCol space
numImages = length(imagepaths);
for i = 1:numImages
    imname = imagepaths{i};
    raw_image = imread([datadir filesep imname]);
    h = figure;
    subplot(1,2,1); imshow(raw_image);
    [xsize, ysize] = size(raw_image(:,:,1));
    rgb_image = raw2rgb(raw_image);    
    imtool(raw_image); 
    rgb_image = raw2rgb(raw_image);    
    [ im_oppCol, first_component] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options);
    im_brightness = reshape(first_component,[xsize, ysize]);
    im_brightness(im_brightness<255) = 0;
    imtool(im_brightness)

%     [ img, ~] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
%     % plot the color space 
%     nstep = 100;
%     subplot(1,2,2)
%     scatter(img(1,1:nstep:end),img(2,1:nstep:end),2,rgb_image(:,1:nstep:end)'./255,'filled');
%     axis([-1 1 -1 1]);
%     axis equal
%     % save the plot
%     split_string = regexp(imname,'\.','split');
%     savename = fullfile(resultdir,split_string{1});
%     print(h,'-dpng', [savename '_oppCol_rep_size2.png']);
%     
%     close all;
end



    