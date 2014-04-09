% deconvolve_image
% called in masterscript 
% Analyze image whose name was input by user
% this is a script being called inside the master script

raw_image = imread([datadir filesep imname]);
options = struct('PlotResults','on',...
    'filterOD',0.15,...
    'extremeCutoff',1);
[ purple_stain_man,pink_stain_man] = deconvolutionManual( imname,datadir,resultdir);% options ); % why do I need workdir hmm
[ stain_mat_svd, saturation_mat_svd, project2svds, U] = deconvolutionSVD(imname,datadir,resultdir);%,options);
[ stain_mat_nnmf, saturation_mat_nnmf] = deconvolutionNNMF(imname,datadir,resultdir);%, options);
%plot_color_cloud;