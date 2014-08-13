% Need to stitch the images together and then do the analysis
% what do people do with such large image
% first read the image, 
% resize the image
% then stitch them together.

% okay, first of all, find out how many images for 1 svs
% addpath(pwd); %close all;
% 
% workdir = '/Users/lun5/Research/color_deconvolution'; 
% 

%svs_fname = 'tp10-958-1'; % this will now be the input

function [ stitched_image ] = stitchSVS( workdir, svs_fname, resize_ratio )

datadir = fullfile(workdir, 'TissueImages');
fileNames = dir(fullfile(datadir,[svs_fname '*.' 'tif']));
imagepaths = {fileNames.name}';

resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('number of images is '); 
numImages = length(imagepaths)% 420
xarray = []; yarray = [];
% read the coordinate from the file names
% file name system is: svsname_x coord_y coord_2048_2048.tif 

for i = 1:numImages
    imname = imagepaths{i}; 
    split_imname = regexp(imname,'_','split');
    ycoord = split_imname(2);
    xcoord = split_imname(3);
    xarray = [xarray str2double(xcoord)];        
    yarray = [yarray str2double(ycoord)];
end

svs_image = imread(fullfile(workdir, 'aperio_scans',[svs_fname '.svs']),'Index',2);
imshow(svs_image);

tile_size = 2048*resize_ratio;
stitched_image = zeros(tile_size*length(unique(xarray)),...
    tile_size*length(unique(yarray)),3);
for i = 1:numImages
    imname = imagepaths{i}; imname
    raw_image = imread(fullfile(datadir,imname));
    resize_image = imresize(raw_image,resize_ratio);
    stitched_image(xarray(i)*resize_ratio+1:xarray(i)*resize_ratio+tile_size,...
        yarray(i)*resize_ratio+1:yarray(i)*resize_ratio+tile_size,:) = resize_image;
end

rs_im = uint8(imresize(stitched_image,1));
figure; imshow(rs_im); 
split_string = regexp(svs_fname,'\.','split');
savename = fullfile(resultdir,split_string{1});
%imwrite(stitched_image,[savename '_stitch1to4.tiff'],'Compression','none');

end