function [ source_eq_image ] = oppColNormalization( source_im, target_im, source_rotation_matrix, target_rotation_matrix )
% Color normalization using opponent color space
target_rgb = raw2rgb(target_im);
figure;imshow(target_im);
source_rgb = raw2rgb(source_im);
figure;imshow(source_im);
% take care of the white and black pixels
target_indx_white = mean(target_rgb,1) > 230;
target_indx_black = mean(target_rgb,1) < 30;
target_indx_saturated = target_indx_white | target_indx_black;
source_indx_white = mean(source_rgb,1) > 230;
source_indx_black = mean(source_rgb,1) < 30;
source_indx_saturated = source_indx_white | source_indx_black;

% converting to opponent color spaces
[ target_oppCol, target_brightness, target_theta, target_rotated] = rgb2oppCol( target_rgb, target_rotation_matrix); 
[ source_oppCol, source_brightness, source_theta, source_rotated] = rgb2oppCol( source_rgb, source_rotation_matrix); 

% plot the opponent color space
nstep = 100;
figure;
scatter(target_oppCol(1,1:nstep:end),target_oppCol(2,1:nstep:end),20,target_rgb(:,1:nstep:end)'./255,'filled');
axis([-1 1 -1 1]); axis square

figure;
scatter(source_oppCol(1,1:nstep:end),source_oppCol(2,1:nstep:end),20,source_rgb(:,1:nstep:end)'./255,'filled');
axis([-1 1 -1 1]); axis square

%% Histogram equalization (matching) 
% normalize the range of brightness to between 0 and 1
source_brightness_norm = source_brightness./range(source_brightness);
target_brightness_norm = target_brightness./range(target_brightness);
% count number of elements in each bin 
binranges = 0:0.01:1;
bincounts = histc(target_brightness_norm(~target_indx_saturated),binranges);
% equalize source's brightness --> target's brightness in [0 1]
source_brightness_norm_eq = histeq(source_brightness_norm(~ source_indx_saturated),bincounts);
% multiply by the range of target's brightness to complete normalization
source_brightness_eq = source_brightness;
source_brightness_eq(~ source_indx_saturated) = source_brightness_norm_eq * range(target_brightness);

figure;
subplot(1,3,1); hist(source_brightness);
subplot(1,3,2); hist(target_brightness);
subplot(1,3,3); hist(source_brightness_eq);

%% normalize the angle information
% convert theta's of source's and target's to [0 1]
target_theta_norm = target_theta./(2*pi) + 0.5;
source_theta_norm = source_theta./(2*pi) + 0.5;

% histogram equalization of source --> target
binranges = 0:0.01:1;
bincounts = histc(target_theta_norm(~ target_indx_saturated), binranges);
source_theta_norm_eq = histeq(source_theta_norm(~ source_indx_saturated),bincounts);
% convert back to between -pi and pi
source_theta_eq = source_theta;
source_theta_eq(~ source_indx_saturated) = (source_theta_norm_eq - 0.5).*(2*pi);

% figure; hist(source_theta_eq);
figure;
subplot(1,3,1); hist(source_theta, -pi:pi/6:pi);
subplot(1,3,2); hist(target_theta,-pi:pi/6:pi);
subplot(1,3,3); hist(source_theta_eq, -pi:pi/6:pi);

%% calculate c2, c3, recover rotated coordinate
% then calculate new rgb
radii = sqrt(source_rotated(2,:).^2 + source_rotated(3,:).^2); % calculate r
source_rotated_eq = zeros(3,length(source_brightness));
source_rotated_eq(1,:) = source_brightness_eq; % brightness normalized/equalized
source_rotated_eq(2,:) = radii.*cos(source_theta_eq); % c2
source_rotated_eq(3,:) = radii.*sin(source_theta_eq); % c3
source_rgb_eq = target_rotation_matrix\source_rotated_eq; % RGB = Rot_matrix \ rotated coordinates
%source_rgb_eq = rotation_matrix_source\source_rotated_eq;
source_rgb_eq_uint8 = uint8(source_rgb_eq); 
source_rgb_eq_uint8(:,source_indx_saturated) = source_rgb(:,source_indx_saturated);

[xsize,ysize] = size(source_im(:,:,1));
r = reshape(source_rgb_eq_uint8(1,:),[xsize, ysize]);
g = reshape(source_rgb_eq_uint8(2,:),[xsize, ysize]);
b = reshape(source_rgb_eq_uint8(3,:),[xsize, ysize]);

source_eq_image = cat(3,r,g,b);
figure; imshow(source_eq_image);

end

