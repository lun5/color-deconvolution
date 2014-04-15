% get training data for the pink and purple stains

addpath(pwd);
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/color_deconvolution'; 
%datadir = uigetdir('*.tiff', 'Please select the image folder');
%datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

% disp('Read the imput')
% imagepaths = getImInfo(datadir);
% numImages = length(imagepaths);
% nrepeats = 1;
% training_data = zeros(3,numImages*nrepeats*2);
% count = 1; % number of training examples
% for i = 1:numImages
%     raw_image = imread([datadir filesep imagepaths{i}]);
%     imshow(raw_image);zoom('on'); zoom(2);
%     % select pink
%     for j = 1:nrepeats        
%     %msgbox(['Please select ' num2str(nrepeats) ' regions of purple'],'Success');
%     rect = getrect; rect = floor(rect);
%     purple_image_crop = imcrop(raw_image,rect);
%     purple_manual_rgb = mean(raw2rgb(purple_image_crop),2);
%     training_data(:,count) = purple_manual_rgb;
%     count = count + 1;
%     end 
%     % select purple
%     for j = 1:nrepeats
%     %msgbox(['Please select ' num2str(nrepeats) ' regions of pink'],'Success');
%     rect = getrect; rect = floor(rect);
%     pink_image_crop = imcrop(raw_image,rect);
%     pink_manual_rgb = mean(raw2rgb(pink_image_crop),2);
%     training_data(:,count) = pink_manual_rgb;
%     count = count + 1;
%     end
% end
% 
% save([resultdir filesep 'training_pink_purple.mat'],'training_data');
% 
% get new rotation matrix
[U,D,V] = svd(training_data);
rotation_matrix = [-U(:,1) U(:,2:3)];
datadir = uigetdir('*.tiff', 'Please select the image folder');
prompt = 'Please enter the image name ';
imname = input(prompt,'s');
while isempty(imname)
    imname = input(prompt,'s');
end

raw_image = imread([datadir filesep imname]);
rgb_image_whole = raw2rgb(raw_image);
[xsize, ysize] = size(raw_image(:,:,1));
nsamples = 100000;
indx_pixel = randperm(xsize*ysize,nsamples);
rgb_image = rgb_image_whole(:,indx_pixel);
mu_s = 4;
sigma_s = 2;
sic_image = rgb2sic( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
figure;scatter(sic_image(1,:),sic_image(2,:),20,rgb_image'./255,'filled');
hold on
xlabel('s1');ylabel('s2');
line([-1 1],[0 0]);line([0 0],[-1 1])
hold off
title(['mu_s = ' num2str(mu_s) ' sigma_s = ' num2str(sigma_s)]); axis([-1 1 -1 1]);

    
