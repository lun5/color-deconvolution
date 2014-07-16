% July 1st
% Get training data for images on which the algorithm performs poorly then
% retrain the model to see if it makes a different in segmentation 

close all;
workdir = '/Users/lun5/Research/color_deconvolution'; 
prompt = 'Please enter the image name ';
%datadir = fullfile(workdir, '20x_images','TIFF','setaside');
datadir = fullfile(workdir, '20x_images','bluebrown');

imname = input(prompt,'s');
while isempty(imname)
    imname = input(prompt,'s');
end

% read the raw img
raw_image = imread([datadir filesep imname]);
[xsize, ysize] = size(raw_image(:,:,1));
imshow(raw_image)
rgb_image = raw2rgb(raw_image);

imshow(raw_image);zoom('on'); zoom(1);
nrepeats = 10;
%training_data = zeros(3,numImages*nrepeats*2);
training_data_purple = zeros(3, 1e6);
training_data_pink = zeros(3, 1e6);

count_purple = 1; % number of training examples
count_pink = 1;

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
    
disp('GOING TO CHANGE TO PINK!!!!');
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

training_data_purple(:, count_purple+1:end) = [];
training_data_pink(:, count_pink+1:end) = [];
label_purple = ones(size(training_data_purple,2),1);
label_pink = zeros(size(training_data_pink,2),1);
labels = [label_purple; label_pink];
%Convert from RGB --> OppCol space
options = struct('Normalize','on');
training_data = [training_data_purple training_data_pink];
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one
mu_s = 4; sigma_s = 2; % values for normalization
%rgb_image = [X_purple_rgb X_pink_rgb];
%[ oppCol_image, indx_sat] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
[X_purple_oppCol,~] = rgb2oppCol(training_data_purple, mu_s, sigma_s, rotation_matrix, options);
[X_pink_oppCol,~] = rgb2oppCol(training_data_pink, mu_s, sigma_s, rotation_matrix, options);
X_oppCol = [X_purple_oppCol';X_pink_oppCol'];
%% Train the classificer on the whole training set
[b,dev,stats] = glmfit(X_oppCol,labels,'binomial','link','logit'); 

% train the model
options = struct('Normalize','on');
[ oppCol_image, ~] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
       
% Use the calculated parameters to classify a new image
labels_prob = glmval(b,oppCol_image','logit');
labels_pred = labels_prob > 0.5;
% plot the classified 
% each class is 0.5
indx_bound = abs(labels_prob - 0.5) < 5*1e-2;
% plot the predicted classes of test data        
figure;
plot(oppCol_image(1,labels_pred == 0),oppCol_image(2,labels_pred == 0),'r.', 'MarkerSize',10);
hold on;
plot(oppCol_image(1,labels_pred == 1),oppCol_image(2,labels_pred == 1),'b.', 'MarkerSize',10);
plot(oppCol_image(1,indx_bound),oppCol_image(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
axis([-1 1 -1 1]);
axis equal

nstep = 1;
figure;
scatter(oppCol_image(1,1:nstep:end),oppCol_image(2,1:nstep:end),20,rgb_image(:,1:nstep:end)'./255,'filled');
hold on;
plot(oppCol_image(1,indx_bound),oppCol_image(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
axis([-1 1 -1 1]);
axis equal

%% plot the image based on the probabilities and predictions
prob_map = reshape(labels_prob,[xsize, ysize]);
figure; imagesc(prob_map);
%colormap(hsv)
%colormap(spring)
% colormap(winter)
% colormap(gray)
% colormap(bone)
% colormap(pink)
colormap(lines)
colorbar;

pred_map = reshape(labels_pred,[xsize, ysize]);
figure; %imagesc(pred_map);
imshow(pred_map);

segmented_purple = double(raw_image).*repmat(pred_map,1,1,3);
segmented_purple = uint8(segmented_purple);
figure; imshow(segmented_purple);

raw_cross_prob = double(raw_image).*repmat(prob_map,1,1,3);
raw_cross_prob = uint8(raw_cross_prob);
figure; imshow(raw_cross_prob);

subtract_img = raw_cross_prob - segmented_purple;
figure; imshow(subtract_img);
