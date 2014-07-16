% get training data for the pink and purple stains

addpath(pwd);
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/color_deconvolution'; 
%datadir = uigetdir('*.tiff', 'Please select the image folder');
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('Read the imput')
imagepaths = getImInfo(datadir,'TIF');
numImages = length(imagepaths);
nrepeats = 1;
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
    % read the raw img
    raw_image = imread([datadir filesep imagepaths{i}]);
    imshow(raw_image);zoom('on'); zoom(2);
    % select purple
    for j = 1:nrepeats        
    disp('Please crop a purple region');
    %msgbox(['Please select ' num2str(nrepeats) ' regions of purple'],'Success');
    rect = getrect; rect = floor(rect);
    purple_image_crop = imcrop(raw_image,rect);
    purple_image_crop_rgb = raw2rgb(purple_image_crop);
    npixels = size(purple_image_crop_rgb,2);
    %purple_manual_rgb = mean(raw2rgb(purple_image_crop),2);
    training_data_purple(:,count_purple:count_purple+npixels-1) = purple_image_crop_rgb;
    count_purple = count_purple + npixels;
    end
    
    % select pink
    for j = 1:nrepeats
    disp('Please crop a pink region');
    %msgbox(['Please select ' num2str(nrepeats) ' regions of pink'],'Success');
    rect = getrect; rect = floor(rect);
    pink_image_crop = imcrop(raw_image,rect);
    pink_image_crop_rgb = raw2rgb(pink_image_crop);
    npixels = size(pink_image_crop_rgb,2);
    %pink_manual_rgb = mean(raw2rgb(pink_image_crop),2);
    training_data_pink(:,count_pink:count_pink+npixels-1) = pink_image_crop_rgb;
    count_pink = count_pink + npixels;
    end
end

training_data_purple(:, count_purple+1:end) = [];
training_data_pink(:, count_pink+1:end) = [];

save([resultdir filesep 'training_purple.mat'],'training_data_purple');
save([resultdir filesep 'training_pink.mat'],'training_data_pink');

training_data = [training_data_purple(:,1:4000) training_data_pink(:,1:12000)];
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

% training3 = load([workdir filesep 'results' filesep '140507' filesep  'training_pink_purple.mat'],'training_data');
% training_data = training3.training_data(:,1:15000);
% [U,D,V] = svd(training_data,0);
% rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

rgb_image = [training_data_purple training_data_pink];
[X_purple_oppCol,~] = rgb2oppCol(training_data_purple, mu_s, sigma_s, rotation_matrix, options);
[X_pink_oppCol,~] = rgb2oppCol(training_data_pink, mu_s, sigma_s, rotation_matrix, options);
X_oppCol = [X_purple_oppCol';X_pink_oppCol'];

% apply the classifier on images
nstep = 1;
h1 = figure;
%subplot(1,3,1);
scatter(X_oppCol(1:nstep:end,1),X_oppCol(1:nstep:end,2),20,rgb_image(:,1:nstep:end)'./255,'filled');
axis([-1 1 -1 1])

% how does pink and purple look like?
% purple
h2 = figure;
%subplot(1,3,2);
scatter(X_purple_oppCol(1,:), X_purple_oppCol(2,:),20,training_data_purple'./255,'filled');
axis([-1 1 -1 1])

% pink
h3 = figure;
%subplot(1,3,3);
scatter(X_pink_oppCol(1,:), X_pink_oppCol(2,:),20,training_data_pink'./255,'filled');  
axis([-1 1 -1 1])


