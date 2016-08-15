function [ source_eq_image ] = oppColNormalization( source_im, target_im,...
     source_rotation_matrix, target_rotation_matrix, options, varargin)

defaultopt = struct('matchMethod','cdf'); 

if nargin < 5
    options = [];
    if nargin < 4
          error('Need at least 4 inputs: source image, target image, source rotation matrix, target rotation matrix')        
    end
end 

matchMethod = optimget(options,'matchMethod',defaultopt,'fast');
% Color normalization using opponent color space
%% input images
target_rgb = raw2rgb(target_im);
[target_xsize, target_ysize] = size(target_im(:,:,1));
figure;imshow(target_im);
source_rgb = raw2rgb(source_im);
[source_xsize, source_ysize] = size(source_im(:,:,1));
figure;imshow(source_im);
% take care of the white and black pixels
target_indx_white = mean(target_rgb,1) > 230;
target_indx_black = mean(target_rgb,1) < 30;
target_indx_saturated = target_indx_white | target_indx_black;
source_indx_white = mean(source_rgb,1) > 230;
source_indx_black = mean(source_rgb,1) < 30;
source_indx_saturated = source_indx_white | source_indx_black;

%% converting to opponent color spaces
[ target_oppCol, target_brightness, target_theta, target_sat] = rgb2oppCol( target_rgb, target_rotation_matrix); 
[ source_oppCol, source_brightness, source_theta, source_sat] = rgb2oppCol( source_rgb, source_rotation_matrix); 
%[ source_oppCol, source_brightness, source_theta, source_sat] = rgb2oppCol( source_rgb, target_rotation_matrix);
% plot the opponent color space
nstep = 100;
figure;
subplot(1,2,1);scatter(target_oppCol(1,1:nstep:end),target_oppCol(2,1:nstep:end),20,target_rgb(:,1:nstep:end)'./255,'filled');
axis([-1 1 -1 1]); axis square
subplot(1,2,2);scatter(source_oppCol(1,1:nstep:end),source_oppCol(2,1:nstep:end),20,source_rgb(:,1:nstep:end)'./255,'filled');
axis([-1 1 -1 1]); axis square
% 
figure;
subplot(1,2,1);circ_plot(target_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,2,2);circ_plot(source_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% 
% % 
% % **************************** temporary patched this here
% % ****************************%
% % options = struct('Normalize','on');
% % mu_s = 4; sigma_s = 2; % values for normalization
% % [ source_oppCol, ~] = rgb2oppCol_old( source_rgb, mu_s, sigma_s, source_rotation_matrix, options); 
% % 
%% matching hues, brightness, saturation
if strcmp(matchMethod,'cdf')
    matchHueCdf;
    matchBrightnessCdf;
    matchSatCdf;
elseif strcmp(matchMethod,'stats')
    matchHueStats
    matchBrightnessStats;
    matchSatStats
else
    error('wrong input for the matching method')
end

%% displaying brightness, theta, and saturation in the same size as the image
% display_brightness_theta_sat_im; % a script attached

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
%sat = sqrt(source_rotated(2,:).^2 + source_rotated(3,:).^2); % calculate r
source_rotated_eq = zeros(3,length(source_brightness));
source_rotated_eq(1,:) = source_brightness_eq; % brightness normalized/equalized
source_rotated_eq(2,:) = source_sat_eq.*cos(source_theta_eq); % c2
source_rotated_eq(3,:) = source_sat_eq.*sin(source_theta_eq); % c3
source_rgb_eq = target_rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
%source_rgb_eq = rotation_matrix_source\source_rotated_eq;
source_rgb_eq_uint8 = uint8(source_rgb_eq); 
source_rgb_eq_uint8(:,source_indx_saturated) = source_rgb(:,source_indx_saturated);

r = reshape(source_rgb_eq_uint8(1,:),[source_xsize, source_ysize]);
g = reshape(source_rgb_eq_uint8(2,:),[source_xsize, source_ysize]);
b = reshape(source_rgb_eq_uint8(3,:),[source_xsize, source_ysize]);

source_eq_image = cat(3,r,g,b);
figure; imshow(source_eq_image);

%close all;
[ source_oppCol_eq, source_brightness_eq, source_theta_eq, source_sat_eq] = rgb2oppCol( double(source_rgb_eq_uint8), target_rotation_matrix);

% edge detection for the normalized image. 
edge_detect; 
%% plot the distribution
figure;
subplot(1,3,1);circ_plot(target_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,3,2);circ_plot(source_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
subplot(1,3,3);circ_plot(source_theta_eq(:),'hist',[],40,true,true,'linewidth',2,'color','r');

%test090614;
end

