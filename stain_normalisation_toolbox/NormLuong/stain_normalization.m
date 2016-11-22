
function  normalized_image = stain_normalization( orig_image, pink_purple_mask)

rotation_matrix = load('rotation_matrix_tp10-867-1.mat','rotation_matrix');
numClusters = 3; % only purple and pink this time
opts_mixture.noise = 1;

%% calculate features
which_features = {'hue opp', 'brightness opp','saturation opp'};
%tic;
im_rgb = double(orig_image)./255;
nrows = size(im_rgb,1); ncols = size(im_rgb,2);
X = reshape(im_rgb,[nrows*ncols,3]);
rotated_coordinates = rotation_matrix.rotation_matrix*X';
indx_purple_pink = pink_purple_mask(:);

%% cluster using von Mises
theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
sat = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
brightness = rotated_coordinates(1,:);
% von Mises
X_cart = [cos(theta); sin(theta)]';
X_cart_pp = X_cart(pink_purple_mask(:),:); % pp = purple + pink
% Call the function
%[ mu_hat_polar,~, kappa_hat,posterior_probs_pp, prior_probs_pp] =...
%    moVM_fixWhite(X_cart_pp,numClusters,opts_mixture);
[ mu_hat_polar,~, kappa_hat,posterior_probs_pp, prior_probs_pp] =...
    moVM(X_cart_pp,numClusters,opts_mixture);
num_pixels = length(theta);
% posterior_probs = zeros(num_pixels,4);
% posterior_probs(indx_purple_pink,[1:2 4]) = posterior_probs_pp;
% posterior_probs(~indx_purple_pink,3) = 1; %white, red, mask
% prior_probs_pp = prior_probs_pp/num_pixels *sum(indx_purple_pink);
% prior_probs = [prior_probs_pp(1:2), sum(~indx_purple_pink)/num_pixels,...
%     prior_probs_pp(3)];
posterior_probs = zeros(num_pixels,5);
posterior_probs(indx_purple_pink,[1:3 5]) = posterior_probs_pp;
posterior_probs(~indx_purple_pink,4) = 1; %white, red, mask
prior_probs_pp = prior_probs_pp/num_pixels *sum(indx_purple_pink);
prior_probs = [prior_probs_pp(1:3), sum(~indx_purple_pink)/num_pixels,...
    prior_probs_pp(4)];

opts_matching.source_stats = struct('mu_hat_polar',mu_hat_polar,'kappa_hat',kappa_hat,...
    'posterior_probs',posterior_probs,'prior_probs',prior_probs);
f_maps_source = {reshape(theta,[nrows,ncols]), reshape(brightness,[nrows,ncols]),...
    reshape(sat,[nrows,ncols])};
%fprintf('done with stats of source image in %.2f\n',toc);

tmp = load('target_stats_aug24');%tmp = load('target_stats_new');%tmp = load('target_stats_may19_lowres');
%tmp = load('target_stats_may19');
opts_matching.target_stats = tmp.data;
tmp = load('f_maps_target_aug24');%tmp = load('f_maps_target_new');%tmp = load('f_maps_target_may19_lowres');
%tmp = load('f_maps_target_may19');
f_maps_target = tmp.data;
f_maps_source_normalized = cell(1,3);
for feature_iter = 1:length(which_features)
    %tic;
    f_map_source_curr = f_maps_source{feature_iter};
    f_map_target_curr = f_maps_target{feature_iter};
    %% normalization
    f_map_normalized_curr = matchingMoments(f_map_source_curr, f_map_target_curr,which_features{feature_iter}, opts_matching);
    f_maps_source_normalized{feature_iter} = f_map_normalized_curr;
    %fprintf('done with matching %s in %.2f\n',which_features{feature_iter},toc);
end

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
source_rotated_eq = zeros(3,length(f_map_normalized_curr));
source_rotated_eq(1,:) = f_maps_source_normalized{2}; % brightness normalized/equalized
source_rotated_eq(2,:) = f_maps_source_normalized{3}.*cos(f_maps_source_normalized{1}); % c2
source_rotated_eq(3,:) = f_maps_source_normalized{3}.*sin(f_maps_source_normalized{1}); % c3
source_rgb_eq = rotation_matrix.rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
source_rgb_eq(source_rgb_eq < 0) = 0;
source_rgb_eq(source_rgb_eq > 1) = 1;
source_rgb_eq_uint8 = uint8(source_rgb_eq*255); 
normalized_image = reshape(source_rgb_eq_uint8', size(orig_image));
end





