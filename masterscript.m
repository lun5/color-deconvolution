% Luong Nguyen
% March 20, 2014
% Demo code for color-deconvolution of H+E stained slides
% Data: 30 slides of BC primary tumors from Dr. Lee
% Method: Marc Macenko, Nancy E. Thomas, IEEE 2009

%%
clearvars;
disp('Directories of inputs and results')
workdir = '/Users/lun5/Research/deconvolution'; 
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

disp('Read the imput')
imagepaths = getImInfo(datadir);
disp('Set parameters for automated de-convolution')
extremeCutoff = 1;
filterOD = 0.15;
disp('Analyze the image set');

% deconvolve_image_sequence;
% ask user for name of image, input that into function
% deconvole_image

% debug flag: on off
% plot flag: on off
% put in options struct: parameter-method, flags
% save image to powerpoint: 
% take an optimal pink and an optimal purple. project pixels on the two
% vectors, select the most purple and most pink based on the projection,
% then calculate the saturation. 

prompt = 'Please enter the image name ';
imname = input(prompt,'s');
while isempty(imname)
    imname = input(prompt,'s');
end

deconvolve_image;
