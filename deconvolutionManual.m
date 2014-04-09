function [ purple_manual_rgb,pink_manual_rgb ] = deconvolutionManual( imname, datadir, resultdir, options,  varargin)
%deconvolution_man: manually deconvole colors by selecting regions that are
%most purple and pink respectively
defaultopt = struct('PlotResults','on',...
    'filterOD',0.15); % flag for plotting the results

if nargin < 4
    options = [];
    if nargin < 3
          error('Need at least 3 inputs: image name, data directory, and result directory')
    end
end

% get the extreme cutoff and filter optical density and plot flag
filterOD = optimget(options,'filterOD',defaultopt,'fast');
plotflag = optimget(options,'PlotResults',defaultopt,'fast');

raw_image = imread([datadir filesep imname]);
[xsize, ysize] = size(raw_image(:,:,1));

disp('manually select purple stain vector...');
figure; imshow(raw_image);zoom('on'); zoom(2);
rect = getrect; rect = floor(rect);
purple_image_crop = imcrop(raw_image,rect);
purple_od = rgb2od(purple_image_crop); % convert to OD
purple_manual_od = mean(purple_od,2);
purple_manual_rgb = mean(raw2rgb(purple_image_crop),2);
%
disp('manually select pink stain vector...');
rect = getrect; rect = floor(rect);
pink_image_crop = imcrop(raw_image,rect);
pink_od = rgb2od(pink_image_crop); % convert to OD
pink_manual_od = mean(pink_od,2);
pink_manual_rgb = mean(raw2rgb(pink_image_crop),2);
%
disp('calculate optical density tuple ...');
calculate_optical_density
%
disp('deconvole image by manually selected stain vectors...')
stain_mat_man = [purple_manual_od pink_manual_od]; % do i need to normalize this?
saturation_mat_man = pinv(stain_mat_man)*opticalDensity;
purple_stain_rgb = stainvec2rgb(stain_mat_man(:,1),saturation_mat_man(1,:),xsize,ysize);
pink_stain_rgb = stainvec2rgb(stain_mat_man(:,2),saturation_mat_man(2,:),xsize,ysize);
remain_rgb_man = raw_image - purple_stain_rgb - pink_stain_rgb;

if strcmpi(plotflag,'on')
    plot_manual_results;
end

end

