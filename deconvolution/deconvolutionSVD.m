%colorDeconvolution: deconvole H+E images into blue + pink stained images
%   Algorithm: Macenko et al. A method for normalizing histology slides for
%   quantitative analysis
% INPUTS: varargin should include; we actually only need stain_mat,
% saturation_mat
%         imname, - path to image + imagename
%         filterOD - threshold for low stain (lowOD) - default 0.15
%         extremeCutoff - to account for noise noise - default 1
% OUTPUTS: saturationMat and matrix of stain vectors stainVector
%[imname, filterOD, extremeCutoff] = parse_inputs(varargin{:});
% Author: Luong Nguyen
% Modified 4/9/14
function [ stain1_svd, stain2_svd, saturation_mat] = deconvolutionSVD(imname,datadir,resultdir, options, varargin) 

defaultopt = struct('PlotResults','on',...
    'filterOD',0.15,...
    'extremeCutoff',1); % flag for plotting the results

if nargin < 4
    options = [];
    if nargin < 3
          error('Need at least 3 inputs: image name, data directory, and result directory')
    end
end

%datadir = fullfile(workdir, '20x_images','TIFF'); % need to move this outside 
raw_image = imread([datadir filesep imname]);
[xsize, ysize] = size(raw_image(:,:,1));
rgb_image = raw2rgb(raw_image);

% get the extreme cutoff and filter optical density and plot flag
filterOD = optimget(options,'filterOD',defaultopt,'fast');
extremeCutoff = optimget(options,'extremeCutoff',defaultopt,'fast');
plotflag = optimget(options,'PlotResults',defaultopt,'fast');

%num_pixels = numel(raw_image(:,:,1));
calculate_optical_density;
[U, S, V] = svd(opticalDensity,'econ'); % calculate SVD from OD
% project onto the first 2 svd directions normalize to unit length
project2svds = S(1:2,1:2)*V(:,1:2)';
% cosine of angles to the first svd directions is first row of
% norm_project2svds, sine is the second row
tan_angles = project2svds(2,:)./project2svds(1,:);
angles_1st_svd = atan(tan_angles);
%%
% extreme values of angles 
extreme_values = prctile(angles_1st_svd,[extremeCutoff 100-extremeCutoff]);
indx_min = abs(angles_1st_svd - extreme_values(1)) < 1e-5;
indx_max = abs(angles_1st_svd - extreme_values(2)) < 1e-5;
proj_stain_mat = [median(project2svds(:,indx_min),2) median(project2svds(:,indx_max),2)];
% OD = VS where V=stainVectors, S = saturationMat
stain_mat = [median(opticalDensity(:,indx_min),2) median(opticalDensity(:,indx_max),2)];
% keep the stain vectors in RGB form
stain1_svd = od2rgb(stain_mat(:,1),1,1);
stain2_svd = od2rgb(stain_mat(:,2),1,1);
saturation_mat = pinv(stain_mat)*opticalDensity;
stain1_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
stain2_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
% stain_mat = normc([median(opticalDensity(:,indx_min),2) median(opticalDensity(:,indx_max),2)]);
% saturation_mat = pinv(stain_mat)*opticalDensity;
% min_stain_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
% max_stain_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
remain_rgb = raw_image - stain1_rgb - stain2_rgb;

if strcmpi(plotflag,'on')
    plot_svd_results;
end

end


