% deconvolve_image_sequence
% script name with _
% function name NO _
% called in masterscript 
% Analyze image sequence
% this is a script being called inside the master script
numImages = length(imagepaths);
min_stainvec_auto = zeros(3,numImages);
max_stainvec_auto = zeros(3,numImages);
purple_stainvec_man = zeros(3,numImages);
pink_stainvec_man = zeros(3,numImages);
for i = 1:numImages
    raw_image = imread([datadir filesep imagepaths{i}]);
    [ stain_mat_auto, saturation_mat_auto, project2svds, U] = deconvolutionSVD(imagepaths{i},workdir, filterOD, extremeCutoff);
    min_stainvec_auto(:,i) = stain_mat_auto(:,1);
    max_stainvec_auto(:,i) = stain_mat_auto(:,2);
    [ purple_stain_man,pink_stain_man] = deconvolutionManual( imagepaths{i},workdir, filterOD );
    purple_stainvec_man(:,i) = purple_stain_man;
    pink_stainvec_man(:,i) = pink_stain_man;
    plot_color_cloud;
end

plotResults;

