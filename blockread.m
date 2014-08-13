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

options = struct('PlotResults','off',...
    'filterOD',0.15,...
    'extremeCutoff',2,...
    'SavePlots','on');

prompt = 'Please enter the image name ';
imname = input(prompt,'s');
while isempty(imname)
    imname = input(prompt,'s');
end

%% get the image name for plotting
split_string = regexp(imname,'\.','split');
savename = fullfile(resultdir,split_string{1});

%% calculate the saturation levels
[ stain1_sic, stain2_sic, saturation_mat] = deconvolutionColOpp(imname,datadir,resultdir,rotation_matrix, options);
step = 100;
% Plot the saturation histograms for each of the stains. 
%% distribution of stain saturations (
h=figure;hist(saturation_mat(1,1:step:end),30);
%b1 = bar(hist(sic_image(1,:),30) ./ sum(sic_image(1,:)));
h3 = findobj(gca,'Type','patch');
set(h3,'FaceColor',[.8 .8 .8],'EdgeColor','w')
title('Distribution of saturation levels of stain 1','FontSize',15);
set(gca,'FontSize',15);
print(h,'-dtiff', [savename '_sat1_dist.tiff']);
