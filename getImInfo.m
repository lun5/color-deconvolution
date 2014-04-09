function [ imagepaths ] = getImInfo( datadir )
%getImInfo gets the paths to images from data dir
% INPUT: datadir - where the images are located
% OUTPUT: imagepaths - path to image
% We will add more to the output once there are more images to analyze

fileNames = dir(fullfile(datadir,'*.TIF'));
imagepaths = {fileNames.name}';

end

