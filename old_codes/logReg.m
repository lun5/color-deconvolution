% Logistic regression model for classification of purple and pink pixels
% Luong Nguyen
% June 23, 2014
% modified July 7, 2014

function logReg(color_space) %(imname,datadir,resultdir, rotation_matrix, options, varargin)
close all;
workdir = '/Users/lun5/Research/color_deconvolution'; 
% defaultopt = struct('PlotResults','on',...
%     'filterOD',0.15,...
%     'extremeCutoff',1,...
%     'SavePlots','on'); % flag for plotting the results
% 
% if nargin < 5
%     options = [];
%     if nargin < 4
%         rotation_matrix = [];
%         if nargin <3
%           error('Need at least 3 inputs: image name, data directory, and result directory')
%         end
%     end
% end

if nargin < 1
    error('Need at least 1 input: color space')
end
    
[ training_data, labels, rotation_matrix ] = import_training_data( color_space );

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

% convert the input image to appropriate color space
% if rgb only then return
if strcmpi(color_space,'RGB')
    img = rgb_image;
elseif strcmpi(color_space,'oppCol') % if opponent color space
    options = struct('Normalize','on');
    mu_s = 4; sigma_s = 2; % values for normalization
    [ img, ~] = rgb2oppCol( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
elseif strcmpi(color_space, 'HSV') % if hsv space
    img = rgb2hsv(raw_image);
    img = raw2rgb(img);
elseif strcmpi(color_space, 'Lab') % if L*a*b* space
    colorTransform = makecform('srgb2lab');
    img = double(applycform(raw_image, colorTransform));
    img = raw2rgb(img);
else
    error('Wrong color space input');    
end

%% Train the classificer on the whole training set
[b,dev,stats] = glmfit(training_data,labels,'binomial','link','logit');        
% Use the calculated parameters to classify a new image
labels_prob = glmval(b,img','logit'); thres = 0.5;
labels_pred = labels_prob > thres;
% plot the classified 
% each class is 0.5
indx_bound = abs(labels_prob - thres) < 5*1e-2;
% plot the predicted classes of test data        
figure;
plot(img(1,labels_pred == 0),img(2,labels_pred == 0),'r.', 'MarkerSize',10);
hold on;
plot(img(1,labels_pred == 1),img(2,labels_pred == 1),'b.', 'MarkerSize',10);
%plot(img(1,indx_bound),img(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
%axis([-1 1 -1 1]);
%axis equal

nstep = 100;
figure;
scatter(img(1,1:nstep:end),img(2,1:nstep:end),20,rgb_image(:,1:nstep:end)'./255,'filled');
hold on;
%plot(img(1,indx_bound),img(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
hold off;
%axis([-1 1 -1 1]);
%axis equal

% h3 = figure;
% scatter3(rgb_image(1,1:nstep:end),rgb_image(2, 1:nstep:end),rgb_image(3,1:nstep:end),20,rgb_image(:, 1:nstep:end)'./255,'filled');
% hold on;
% plot3(rgb_image(1,indx_bound),rgb_image(2, indx_bound),rgb_image(3, indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
% hold off;

%% plot the image based on the probabilities and predictions
prob_map = reshape(labels_prob,[xsize, ysize]);
figure; imagesc(prob_map);
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

end

%% extra
% %% Plot the most purple and most pink pixels from the classification
% % extremeCutoff = 1;
% % extreme_values = prctile(labels_prob,[extremeCutoff 100-extremeCutoff]);
% % % find the stain vectors that are extreme on first coordinate of s1
% % % get index saturation indx_sat from calculate_SIC
% % indx_min = (abs(img(1,:) - extreme_values(1)) < 1e-2) & (~ indx_sat);
% % indx_max = (abs(img(1,:) - extreme_values(2)) < 1e-2) & (~ indx_sat);
% % most pink prob ~ 0
% %most_pink_pixels_indx = abs(labels_prob - extreme_values(1)) < 1e-2;
% most_pink_pixels_indx = abs(labels_prob - 0) < 1e-5;
% % most purple prob ~ 1
% %most_purple_pixels_indx = abs(labels_prob - extreme_values(2)) < 1e-2;
% most_purple_pixels_indx = abs(1 - labels_prob) < 1e-3;
% % plot these
% % figure;
% % scatter(img(1,most_pink_pixels_indx),img(2,most_pink_pixels_indx),20,rgb_image(:,most_pink_pixels_indx)'./255,'filled');
% % hold on;
% % scatter(img(1,most_purple_pixels_indx),img(2,most_purple_pixels_indx),20,rgb_image(:,most_purple_pixels_indx)'./255,'filled');
% % plot(img(1,indx_bound),img(2,indx_bound),'k-', 'MarkerSize',30,'LineWidth',5); % dec boundary
% % hold off;
% % axis([-1 1 -1 1]);
% % axis equal
% 
% % stain vectors:
% pink_oppCol = median(rgb_image(:,most_pink_pixels_indx),2);
% purple_oppCol = median(rgb_image(:,most_purple_pixels_indx),2);
% 
% stain_mat = rgb2od([pink_oppCol purple_oppCol]);
% 
% % get the extreme cutoff and filter optical density and plot flag
% filterOD = 0.15;
% extremeCutoff = 1;
% plotflag = 'on';
% 
% % calculate OD space for deconvolution 
% calculate_optical_density;
% saturation_mat = pinv(stain_mat)*opticalDensity;
% 
% stain1_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
% stain2_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
% remain_rgb = raw_image - stain1_rgb - stain2_rgb;
% 
% if strcmpi(plotflag,'on')
%     % deconvolved images
%     h = figure;
%     subplot(2,2,1)
%     imshow(raw_image)
%     subplot(2,2,2)
%     imshow(stain1_rgb)
%     subplot(2,2,3)
%     imshow(stain2_rgb)
%     subplot(2,2,4)
%     imshow(remain_rgb)
%     %if strcmpi(saveflag,'on')
%     %    print(h,'-dtiff', [resultdir filesep imname '_sic_deconv_' num2str(extremeCutoff) '.tiff']);
%     %end
% end
% 
