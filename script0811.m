% filter out optically saturated (white) pixels
% chemically saturated (black) pixels
% need to load in an image with both white and pink 
addpath(pwd); close all;

workdir = '/Users/lun5/Research/color_deconvolution';
datadir = fullfile(workdir, 'TissueImages');

target_svs = 'tp10-611';
all_white_imname = 'tp10-611_55296_26624_2048_2048.tif';
all_pink_imname = 'tp10-611_53248_14336_2048_2048.tif';
%half_white_pink_imname = 'tp10-611_53248_18432_2048_2048.tif';
half_white_pink_imname = 'tp10-867-1_47104_22528_2048_2048.tif';
%half_white_pink_imname = 'tp10-611_24576_18432_2048_2048.tif';
all_white_im = double(imread(fullfile(datadir, all_white_imname)));
mean_intensities = (all_white_im(:,:,1)+all_white_im(:,:,2)+ all_white_im(:,:,3))./3;
indx = mean_intensities > 230;

sum(indx(:))

im = double(imread(fullfile(datadir, all_pink_imname)));
figure; imshow(uint8(im));
mean_intensities = (im(:,:,1)+im(:,:,2)+ im(:,:,3))./3;
indx = mean_intensities > 230;
indx = cat(3,indx,indx,indx);
im(indx) = 0;
im = uint8(im);
figure; imshow(im);

im = double(imread(fullfile(datadir, half_white_pink_imname)));
figure; imshow(uint8(im));
mean_intensities = (im(:,:,1)+im(:,:,2)+ im(:,:,3))./3;
indx = mean_intensities > 230;
indx = cat(3,indx,indx,indx);
im(indx) = 0;
im = uint8(im);
figure; imshow(im);

% so 200 seems to be good

%% Filter chemically saturated (black) pixels
source_imname = 'tp10-867-1';
tissue_fold_imname = 'tp10-867-1_26624_26624_2048_2048.tif';

im = imread(fullfile(datadir,tissue_fold_imname));
figure; imshow(im); 
mean_intensities = (double(im(:,:,1)) + double(im(:,:,2)) + double(im(:,:,3)))/3;
indx = mean_intensities < 20;
indx = cat(3,indx,indx,indx);
im(indx) = 255;
figure; imshow(im);

%% For an image that have both black & white
black_white_imname = 'tp10-867-1_38912_20480_2048_2048.tif';
im = imread(fullfile(datadir,black_white_imname));

figure; imshow(im); 
mean_intensities = (double(im(:,:,1)) + double(im(:,:,2)) + double(im(:,:,3)))/3;
indx = mean_intensities < 20;
indx_black = cat(3,indx,indx,indx);
im(indx_black) = 255;
indx = mean_intensities > 230;
indx_white = cat(3,indx,indx,indx);
im(indx_white) = 0;
figure; imshow(im);