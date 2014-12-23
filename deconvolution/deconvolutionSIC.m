%colorDeconvolution: deconvole H+E images into blue + pink stained images
%   Algorithm: ElMa
% INPUTS: 
%         imname, - path to image + imagename
%         filterOD - threshold for low stain (lowOD) - default 0.15
%         extremeCutoff - to account for noise noise - default 1
% OUTPUTS: saturationMat and matrix of stain vectors stainVector
% Author: Luong Nguyen
% Modified 4/9/14
function [ stain1_sic, stain2_sic, saturation_mat] = deconvolutionSIC(imname,datadir,resultdir, options, varargin) 

defaultopt = struct('PlotResults','on',...
    'filterOD',0.15,...
    'extremeCutoff',1,...
    'SavePlots','on'); % flag for plotting the results

if nargin < 4
    options = [];
    if nargin < 3
          error('Need at least 3 inputs: image name, data directory, and result directory')
    end
end

% get the extreme cutoff and filter optical density and plot flag
filterOD = optimget(options,'filterOD',defaultopt,'fast');
extremeCutoff = optimget(options,'extremeCutoff',defaultopt,'fast');
plotflag = optimget(options,'PlotResults',defaultopt,'fast');
saveflag = optimget(options,'SavePlots',defaultopt,'fast');

raw_image = imread([datadir filesep imname]);
[xsize, ysize] = size(raw_image(:,:,1));
rgb_image = raw2rgb(raw_image);
% calculate OD space for deconvolution 
calculate_optical_density;
% calculate SIC coordinates
calculate_SIC;

% identify the mode for sic_image(1,:) and sic_image(2,:)
% [N,bin] = hist(sic_image(1,:),100);
% mode_value1 = mode(bin);
% indx_stain1 = abs(sic_image(1,:) - mode_value1) < 1e-2;
% stain1_sic = median(rgb_image(:,indx_stain1),2);
% 
% % identify the mode for sic_image(1,:) and sic_image(2,:)
% [N,bin] = hist(sic_image(2,:),100);
% mode_value2 = mode(bin);
% indx_stain2 = abs(sic_image(2,:) - mode_value2) < 1e-2;
% stain2_sic = median(rgb_image(:,indx_stain2),2);

extreme_values = prctile(sic_image(1,:),[extremeCutoff 100-extremeCutoff]);
% find the stain vectors that are extreme on first coordinate of s1
indx_min = abs(sic_image(1,:) - extreme_values(1)) < 1e-2;
indx_max = abs(sic_image(1,:) - extreme_values(2)) < 1e-2;
% stain vectors in rgb space
stain1_sic = median(rgb_image(:,indx_min),2);
stain2_sic = median(rgb_image(:,indx_max),2);

% calculate OD space for deconvolution 
% rgb_image = raw2rgb(raw_image);
% calculate_optical_density;

% OD = VS where V=stainVectors, S = saturationMat
stain_mat = rgb2od([stain1_sic stain2_sic]);
saturation_mat = pinv(stain_mat)*opticalDensity;

stain1_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
stain2_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
remain_rgb = raw_image - stain1_rgb - stain2_rgb;

if strcmpi(plotflag,'on')
    % deconvolved images
    h = figure;
    subplot(2,2,1)
    imshow(raw_image)
    subplot(2,2,2)
    imshow(stain1_rgb)
    subplot(2,2,3)
    imshow(stain2_rgb)
    subplot(2,2,4)
    imshow(remain_rgb)
    if strcmpi(saveflag,'on')
        print(h,'-dtiff', [resultdir filesep imname '_sic_deconv.tiff']);
    end
end

%close all
end


