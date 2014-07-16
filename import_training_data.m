function [ training_data, labels, rotation_matrix ] = import_training_data( color_space, varargin )
%import_training_data Prepare training data in different color spaces
% Options are: RGB, oppCol, HSV, Lab
workdir = '/Users/lun5/Research/color_deconvolution'; 
trainingdir = fullfile(workdir, 'results', '140625');
% Read the input in RGB form
training_purple = load([trainingdir filesep 'training_purple.mat'],'training_data_purple');
training_pink = load([trainingdir filesep 'training_pink.mat'],'training_data_pink');
X_purple_rgb = training_purple.training_data_purple;
X_pink_rgb = training_pink.training_data_pink;
X_rgb = [X_purple_rgb'; X_pink_rgb'];
% Read in the labels
label_purple = ones(size(X_purple_rgb,2),1);
label_pink = zeros(size(X_pink_rgb,2),1);
labels = [label_purple; label_pink];
% parse the input argument
if nargin ~= 1
    error('Need exactly one input for the color space');
end

if ((nargout < 3) ||  strcmpi(color_space,'RGB') || strcmpi(color_space, 'HSV')...
        || strcmpi(color_space, 'Lab') )
   rotation_matrix = [];
end

rgb_3d = reshape(X_rgb,[length(X_rgb) 1 3]);
% if rgb only then return
if strcmpi(color_space,'RGB')
    training_data = X_rgb;
elseif strcmpi(color_space,'oppCol') % if opponent color space
    options = struct('Normalize','on');
    training_data = [X_purple_rgb(:,1:3000) X_pink_rgb(:,1:9000)];
    [U,D,V] = svd(training_data,0);
    rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one
    mu_s = 4; sigma_s = 2; % values for normalization
    [X_purple_oppCol,~] = rgb2oppCol(X_purple_rgb, mu_s, sigma_s, rotation_matrix, options);
    [X_pink_oppCol,~] = rgb2oppCol(X_pink_rgb, mu_s, sigma_s, rotation_matrix, options);
    training_data = [X_purple_oppCol';X_pink_oppCol'];

elseif strcmpi(color_space, 'HSV') % if hsv space
    training_data = rgb2hsv(rgb_3d);
    training_data = reshape(training_data, [size(training_data,1) 3]);
elseif strcmpi(color_space, 'Lab') % if L*a*b* space
    colorTransform = makecform('srgb2lab');
    training_data = double(applycform(uint8(rgb_3d), colorTransform));
    training_data = reshape(training_data, [size(training_data,1) 3]);
else
    error('Wrong color space input');    
end

end

% convert from RGB -> OD before -> oppCol 