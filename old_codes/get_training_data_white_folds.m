
% Need to sample pixels from these images
% calculate svd from all of these
% what is the benefit
% there are 420 images
% we need to loop through them, 
% omit the one that are just white
% omit the one that are part white
% we probably need to train the model again
% how to sample our training data

% i can loop through them all
% close the white one
% or train also which one is white and which one isn't
% if it is then take some samples from it
% actually we need a way to collect this training data for not good piece

% loop through the images
% click on yes/no for enough info + no blemish 
workdir = '/Users/lun5/Research/color_deconvolution'; 
svs_fname = 'tp10-958-1'; % this will now be the input
resize_ratio = 1/4;
%[ stitched_image ] = stitchSVS( workdir, svs_fname, resize_ratio );

datadir = fullfile(workdir, 'TissueImages');
fileNames = dir(fullfile(datadir,[svs_fname '*.' 'tif']));
imagepaths = {fileNames.name}';

resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

numImages = length(imagepaths)% 420
% read the coordinate from the file names
% file name system is: svsname_x coord_y coord_2048_2048.tif 

% image types for classification
imlabels = zeros(500,1);
for i = 1:numImages
    imname = imagepaths{i}; 
    raw_image = imread(fullfile(datadir,imname));
    imshow(raw_image);
    image_types = {'white','half tissue','full tissue'};
    choice = menu('Choose an image type','white','half tissue','full tissue');
    imlabels(i) = choice;    
end

imlabels(imlabels==0) = [];
% need to save the results
save([resultdir filesep 'labels_tiles.mat'],'imlabels');


