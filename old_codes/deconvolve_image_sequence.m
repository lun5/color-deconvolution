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
for i = 1:10 %numImages
    imname = imagepaths{i};
    deconvolve_image;
    close all;
end

%plotResults;

