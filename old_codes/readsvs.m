% readsvs.m
% Luong Nguyen
% Jun 18, 2014
% read in WSI from Aperio
% these files are from Adrian Lee's lab
% plot them side by side to show some variability in colors. 

addpath(pwd);close all;

workdir = '/Users/lun5/Research/color_deconvolution'; 
datadir = fullfile(workdir, 'aperio_scans');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

imagepaths = getImInfo(datadir,'svs');
numImages = length(imagepaths);

h = figure;
for i = 1:5
    for j = 1:6
        imname = imagepaths{(i-1)*6+j}; imname
        im = imread(fullfile(datadir,imname),'Index',2); 
        subplot(5,6,(i-1)*6+j); imshow(im)
    end
end

savefig(h,fullfile(resultdir,'allWSI.fig'));

% prompt = 'Please enter the image name ';
% imname = input(prompt,'s');
% while isempty(imname)
%     imname = input(prompt,'s');
% end

% im2 = imread(fullfile(datadir,imname),'Index',2); figure; imshow(im2);
% im3 = imread(fullfile(datadir,imname),'Index',3); figure; imshow(im3);
% % 
% im4 = imread(fullfile(datadir,imname),'Index',4); figure; imshow(im4);
% im5 = imread(fullfile(datadir,imname),'Index',5); figure; imshow(im5);
% im6 = imread(fullfile(datadir,imname),'Index',6); figure; imshow(im6);
