% gather training data for images
% loop through the images
% ask weather it's white/not
% if not then start collecting the data
% we can put a patch of pink and a patch of purple
% until there are about 30 qualified images
% then we can stop
% let's write this as a function

function [ training_data_purple, training_data_pink, rotation_matrix] = wsi_get_training( workdir, svs_fname) 

datadir = fullfile(workdir, 'TissueImages');
fileNames = dir(fullfile(datadir,[svs_fname '*.' 'tif']));
imagepaths = {fileNames.name}';

resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

numImages = length(imagepaths);% 420
numImagesQualified = 0;

svs_image = imread(fullfile(workdir, 'aperio_scans',[svs_fname '.svs']),'Index',2);
imshow(svs_image);

%training_data = zeros(3,numImages*nrepeats*2);
training_data_purple = zeros(3, 1e6);
training_data_pink = zeros(3, 1e6);

count_purple = 1; % number of training examples
count_pink = 1;
for i = 1:numImages
    % double the size of training_data_purple if needed
    if count_purple + 1e3 > size(training_data_purple,2)
        temp = training_data_purple;
        training_data_purple = zeros(3,size(training_data_purple,2)*2);
        training_data_purple(:, size(temp,2)) = temp;
    end
    
    % double the size of training_data_pink if needed
    if count_pink + 1e3 > size(training_data_pink,2)
        temp = training_data_pink;
        training_data_pink = zeros(3,size(training_data_pink,2)*2);
        training_data_pink(:, size(temp,2)) = temp;
    end
    
    imname = imagepaths{i}; 
    raw_image = imread(fullfile(datadir,imname));
    imshow(raw_image);
    choice = menu('Image qualified for traning?','white','qualified');
    if choice == 2
        % get the training data
        % select purple
        disp('Please crop a purple region');
        zoom on;
        rect = getrect; 
        zoom on;
        rect = floor(rect);
        purple_image_crop = imcrop(raw_image,rect);
        purple_image_crop_rgb = raw2rgb(purple_image_crop);
        npixels = size(purple_image_crop_rgb,2);
        training_data_purple(:,count_purple:count_purple+npixels-1) = purple_image_crop_rgb;
        count_purple = count_purple + npixels;
    
        % select pink
        disp('Please crop a pink region');
        rect = getrect; rect = floor(rect);
        pink_image_crop = imcrop(raw_image,rect);
        pink_image_crop_rgb = raw2rgb(pink_image_crop);
        npixels = size(pink_image_crop_rgb,2);
        training_data_pink(:,count_pink:count_pink+npixels-1) = pink_image_crop_rgb;
        count_pink = count_pink + npixels;
        % update
        numImagesQualified = numImagesQualified + 1;                        
    end
    % only consider up to 20 qualified images
    if numImagesQualified > 20
        break;
    end
        
end

% remove excessive space in the pink/purple arrays
training_data_purple(:, count_purple+1:end) = [];
training_data_pink(:, count_pink+1:end) = [];
% save the training data for that image
save([resultdir filesep svs_fname 'training_purple.mat'],'training_data_purple');
save([resultdir filesep svs_fname 'training_pink.mat'],'training_data_pink');

% should I or should I not calculate the rotation matrix?
% then I should do the conversion to see as well, without normalization
training_data = [training_data_purple(:,1:min(2000,size(training_data_purple,2)))...
    training_data_pink(:,1:min(8000,size(training_data_pink,2)))];
[U,~,~] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

end