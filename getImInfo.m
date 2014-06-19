function [ imagepaths ] = getImInfo( datadir, filetype )
%getImInfo gets the paths to images from data dir
% INPUT: datadir - where the images are located, filetype is the format to
% read
% OUTPUT: imagepaths - path to image
% We will add more to the output once there are more images to analyze

%fileNames = dir(fullfile(datadir,'*.TIF'));
fileNames = dir(fullfile(datadir,['*.' filetype]));
imagepaths = {fileNames.name}';

end

