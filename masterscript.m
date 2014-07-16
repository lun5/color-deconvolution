% Luong Nguyen
% March 20, 2014
% Demo code for color-deconvolution of H+E stained slides
% Data: 30 slides of BC primary tumors from Dr. Lee
% Method: Marc Macenko, Nancy E. Thomas, IEEE 2009
% Method: NNMF
% Method: SIC
% New method: color opponency
% May 26, 2014
% LNguyen modified the method to show color opponency space instead of
% saturation independent coordinates. 
%%
%clearvars;
addpath(pwd);
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/color_deconvolution'; 
%datadir = uigetdir('*.tiff', 'Please select the image folder');
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));
close all;
if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('Read the input')
imagepaths = getImInfo(datadir,'TIF');
disp('Set parameters for automated de-convolution')
extremeCutoff = 1;
filterOD = 0.15;
disp('Analyze the image set');

%deconvolve_image_sequence;

prompt = 'Please enter the image name ';
imname = input(prompt,'s');
while isempty(imname)
    imname = input(prompt,'s');
end

deconvolve_image;
