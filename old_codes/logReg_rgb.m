% Logistic regression model for classification of purple and pink pixels
% Luong Nguyen
% June 23, 2014

% Work directory and training data directory
close all;
workdir = '/Users/lun5/Research/color_deconvolution'; 
trainingdir = fullfile(workdir, 'results', '140625');
% Read in the labels and the input
training_purple = load([trainingdir filesep 'training_purple.mat'],'training_data_purple');
training_pink = load([trainingdir filesep 'training_pink.mat'],'training_data_pink');
X_purple_rgb = training_purple.training_data_purple;
X_pink_rgb = training_pink.training_data_pink;
X_rgb = [X_purple_rgb'; X_pink_rgb'];

label_purple = ones(size(X_purple_rgb,2),1);
label_pink = zeros(size(X_pink_rgb,2),1);
labels = [label_purple; label_pink];

%% input an image
prompt = 'Please enter image type, either aperio or 20x: ';
imtype = input(prompt,'s');
if strcmp(imtype,'20x')
    %% 20X
    prompt = 'Please enter the image name ';
    datadir = fullfile(workdir, '20x_images','TIFF');%,'setaside');
    imname = input(prompt,'s');
    while isempty(imname)
        imname = input(prompt,'s');
    end

    raw_image = imread([datadir filesep imname]);
    imshow(raw_image);
    [xsize, ysize] = size(raw_image(:,:,1));
    rgb_image = raw2rgb(raw_image);
elseif strcmp(imtype,'aperio')
    %% Aperio
    prompt = 'Please enter the image name ';
    datadir = fullfile(workdir, 'aperio_scans');
    imname = input(prompt,'s');
    while isempty(imname)
        imname = input(prompt,'s');
    end
    raw_image = imread(fullfile(datadir,imname),'Index',2); 
    imshow(raw_image);
    rect = getrect; 
    cropped_image = imcrop(raw_image,rect); 
    figure; imshow(cropped_image);
    [xsize, ysize] = size(cropped_image(:,:,1));
    rgb_image = raw2rgb(cropped_image);
else
    error('Wrong input for image type');
end

% options = struct('Normalize','on');
% [ oppCol_image, ~] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 

%% Train the classificer on the whole training set
[b,dev,stats] = glmfit(X_rgb,labels,'binomial','link','logit');        
% Use the calculated parameters to classify a new image
labels_prob = glmval(b,rgb_image','logit'); thres = 0.5;
labels_pred = labels_prob > thres;
% plot the classified 
% each class is 0.5
indx_bound = abs(labels_prob - thres) < 5*1e-2;
% plot the predicted classes of test data        
h2 = figure;
plot(rgb_image(1,labels(train_indx) == 0),rgb_image(3,labels(train_indx) == 0),'r.', 'MarkerSize',10);
hold on;
plot(rgb_image(1,labels(train_indx) == 1),rgb_image(3,labels(train_indx) == 1),'b.', 'MarkerSize',10);
plot(rgb_image(1,indx_bound),rgb_image(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
%axis([-1 1 -1 1]);
%axis equal

nstep = 100;
h3 = figure;
scatter3(rgb_image(1,1:nstep:end),rgb_image(2, 1:nstep:end),rgb_image(3,1:nstep:end),20,rgb_image(:, 1:nstep:end)'./255,'filled');
hold on;
plot3(rgb_image(1,indx_bound),rgb_image(2, indx_bound),rgb_image(3, indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
%axis([-1 1 -1 1]);
%axis equal
%% plot the image based on the probabilities and predictions
prob_map = reshape(labels_prob,[xsize, ysize]);
figure; imagesc(prob_map);
colorbar;

pred_map = reshape(labels_pred,[xsize, ysize]);
h4 = figure; %imagesc(pred_map);
imshow(pred_map);

segmented_purple = double(raw_image).*repmat(pred_map,1,1,3);
segmented_purple = uint8(segmented_purple);
h5 = figure; imshow(segmented_purple);

raw_cross_prob = double(raw_image).*repmat(prob_map,1,1,3);
raw_cross_prob = uint8(raw_cross_prob);
h6 = figure; imshow(raw_cross_prob);

subtract_img = raw_cross_prob - segmented_purple;
h7 = figure; imshow(subtract_img);
