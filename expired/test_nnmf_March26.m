% Test script March 26
defaultopt = struct('PlotResults','on',...
    'filterOD',0.15); % flag for plotting the results
options = [];
workdir = '/Users/lun5/Research/deconvolution'; 
datadir = fullfile(workdir, '20x_images','TIFF');
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));
imname = 'Slide#1-2.TIF';
raw_image = imread([datadir filesep imname]);
[xsize, ysize] = size(raw_image(:,:,1));
% get the extreme cutoff and filter optical density and plot flag
filterOD = optimget(options,'filterOD',defaultopt,'fast');
plotflag = optimget(options,'PlotResults',defaultopt,'fast');
%calculate optical density
%opticalDensity = raw2rgb(raw_image)./255;
calculate_optical_density 
[stain_mat_nnmf,saturation_mat_nnmf] = nnmf(opticalDensity,2);
min_stain_rgb = stainvec2rgb(stain_mat_nnmf(:,1),saturation_mat_nnmf(1,:),xsize,ysize);
max_stain_rgb = stainvec2rgb(stain_mat_nnmf(:,2),saturation_mat_nnmf(2,:),xsize,ysize);
remain_rgb = raw_image - min_stain_rgb - max_stain_rgb;

if strcmpi(plotflag,'on')
    plot_nnmf_results;
end
