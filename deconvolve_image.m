% deconvolve_image
% called in masterscript 
% Analyze image whose name was input by user
% this is a script being called inside the master script

raw_image = imread([datadir filesep imname]);
options = struct('PlotResults','on',...
    'filterOD',0.15,...
    'extremeCutoff',1,...
    'SavePlots','on');

training3 = load([workdir filesep 'results' filesep '140507' filesep  'training_pink_purple.mat'],'training_data');
%training_data = [training1.training_data training2.training_data];
training_data = training3.training_data(:,1:4000);
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

%[ purple_manual_rgb,pink_manual_rgb, stain_mat_man] = deconvolutionManual( imname,datadir,resultdir);% options ); 
%[ stain1_svd, stain2_svd, saturation_mat_svd] = deconvolutionSVD(imname,datadir,resultdir);%,options);

%options = struct('PlotResults','on',...
%    'filterOD',0.15,'CropImage','off',...
%    'extremeCutoff',1,'ColorSpace','SIC');
%[ stain1_nnmf, stain2_nnmf, saturation_mat_nnmf] = deconvolutionNNMF(imname,datadir,resultdir, options);
%plot_color_cloud;
[ stain1_oppCol, stain2_oppCol, saturation_mat] = deconvolutionOppCol(imname,datadir,resultdir,rotation_matrix, options);
% SIC color representation
% calculate_SIC;