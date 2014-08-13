% opponent color space for the 1x images
% also apply this to the expanded images from bruce 
% sample out 5000 pixels and plotted them on the plane
% Luong Nguyen
% July 29th

addpath(pwd);close all;

workdir = '/Users/lun5/Research/color_deconvolution'; 
datadir = fullfile(workdir, 'aperio_scans');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

imagepaths = getImInfo(datadir,'svs');
numImages = length(imagepaths);

% calculate the rotation matrix
% need to load new training data
trainingdir = fullfile(workdir, 'results', '140625');
% Read the input in RGB form
training_purple = load([trainingdir filesep 'training_purple.mat'],'training_data_purple');
training_pink = load([trainingdir filesep 'training_pink.mat'],'training_data_pink');
X_purple_rgb = training_purple.training_data_purple;
X_pink_rgb = training_pink.training_data_pink;
training_data = [X_purple_rgb(:,1:3000) X_pink_rgb(:,1:9000)];
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

% plot the sic coordinates
options = struct('Normalize','on');
mu_s = 4; sigma_s = 2;

%h = figure;
for i = 1:numImages
    imname = imagepaths{i}; imname
    raw_image = imread(fullfile(datadir,imname),'Index',2);
    imshow(raw_image); rect = getrect; raw_image = imcrop(raw_image,rect);
    xsize = size(raw_image,1); ysize = size(raw_image,2);
    close all
    imshow(raw_image);
    rgb_image = raw2rgb(raw_image);    
    [ im_oppCol, first_component] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options);
    first_component = uint8(first_component/max(first_component) * 255);
    im_brightness = reshape(first_component,[xsize, ysize]);
    %im_brightness(im_brightness<200) = 0;
    imtool(im_brightness)
    
    % plot the opponent color space
    figure;nstep = 1;
    scatter(im_oppCol(1,1:nstep:end),im_oppCol(2,1:nstep:end),2,rgb_image(:,1:nstep:end)'./255,'filled');
    axis([-1 1 -1 1]);
    axis equal
    
end

% Read in the images

% omit the pixels that are white

% sample 5000 pixels

% calculate 