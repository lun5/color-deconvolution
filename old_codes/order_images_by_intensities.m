% Luong Nguyen
% 5/21/14
% Order images by intensities  
% Inputs: small microscope images by Whitney
% Reason: can't read svs file. 33003 decode is not implemented

addpath(pwd);
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/color_deconvolution'; 
%datadir = uigetdir('*.tiff', 'Please select the image folder');
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('Read the input')
imagepaths = getImInfo(datadir);
disp('Set parameters for automated de-convolution')
extremeCutoff = 1;
filterOD = 0.15;
disp('Calculate the rotation matrix')
training1 = load([workdir filesep 'training_pink_purple.mat'],'training_data');
training2 = load([workdir filesep 'results' filesep '140428' filesep  'training_pink_purple.mat'],'training_data');
training3 = load([workdir filesep 'results' filesep '140507' filesep  'training_pink_purple.mat'],'training_data');
%training_data = [training1.training_data training2.training_data];
training_data = training3.training_data(:,1:2000);
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

disp('Analyze the image set');
numImages = length(imagepaths);
stain1_mean_saturation = zeros(1,numImages);
stain2_mean_saturation = zeros(1,numImages);

options = struct('PlotResults','off',...
    'filterOD',0.15,...
    'extremeCutoff',1,...
    'SavePlots','off');
for i = 1:numImages
    disp(['analyze image number ' num2str(i)]);
    imname = imagepaths{i};
    raw_image = imread([datadir filesep imname]);
    [ stain1_sic, stain2_sic, saturation_mat] = deconvolutionColOpp(imname,datadir,resultdir,rotation_matrix, options);
    stain1_mean_saturation(i) = mean(saturation_mat(1,:));
    stain2_mean_saturation(i) = mean(saturation_mat(2,:));
end

% sort the images by mean saturation
[sorted_mean1, indx_stain1] = sort(stain1_mean_saturation);
% 
% for j = indx_stain2
%     raw_image = imread([datadir filesep imagepaths{j}]);
%     figure; imshow(raw_image)
%     pause;
%     close
% end

% plot images by saturation
h1 = figure; plot(1:numImages, sorted_mean1, 'LineWidth',3);
axis([0 100 0 1]);
grid on;
ylabel('Mean saturation of purple stain','FontSize',15);
im = imread([datadir filesep imagepaths{indx_stain1(2)}]);
axes('position',[0.1,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain1(30)}]);
axes('position',[0.3,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain1(60)}]);
axes('position',[0.6,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain1(93)}]);
axes('position',[0.8,0, 0.1, 0.1]);
imshow(im);
print(h1,'-dtiff', [resultdir filesep 'sortPurpleSat' '.tiff']);

[sorted_mean2, indx_stain2] = sort(stain2_mean_saturation);
plot(1:numImages, sorted_mean2)

% for j = indx_stain2
%     raw_image = imread([datadir filesep imagepaths{j}]);
%     figure; imshow(raw_image)
%     pause;
%     close
% end

h2 = figure; plot(1:numImages, sorted_mean2, 'LineWidth',3);
axis([0 100 0 1]);
grid on
ylabel('Mean saturation of pink stain','FontSize',15);
im = imread([datadir filesep imagepaths{indx_stain2(2)}]);
axes('position',[0.1,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain2(30)}]);
axes('position',[0.3,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain2(60)}]);
axes('position',[0.6,0, 0.1, 0.1]);
imshow(im);
im = imread([datadir filesep imagepaths{indx_stain2(93)}]);
axes('position',[0.8,0, 0.1, 0.1]);
imshow(im);
print(h2,'-dtiff', [resultdir filesep 'sortPinkSat' '.tiff']);
