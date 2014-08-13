function calibrate_deconvolution
%calibrate_deconvolution 
% identify the stain vectors
%%
close all; clear all;
clc
disp('read a test image that the algorithm does not perform well on ...'); 
fname = '/Users/lun5/Research/deconvolution/20x_images/TIFF/Slide#2-3.TIF';
%fname = '/Users/lun5/Research/deconvolution/20x_images/TIFF/Slide#1-1.TIF';
raw_image = imread(fname,'tif');
disp('range of intensity')
range(raw_image(:)); filterOD = 0.15;
%[ purple_stain_rgb ,pink_stain_rgb,remain_rgb_man ] = deconvolution_man( raw_image, filterOD );

optical_density_clouds(raw_image,200);
for ntimes = 1: 10 % 
    disp('select a region to examine the distribution of color pixels');
    imshow(raw_image);
    rect = getrect; rect = floor(rect);
    optical_density_clouds(imcrop(raw_image,rect),5);
end

    
% disp('deconvole the image into R, G, B');
% figure;
% subplot(2,2,1)
% imshow(raw_image)
% subplot(2,2,2)
% imshow(raw_image(:,:,1))
% subplot(2,2,3)
% imshow(raw_image(:,:,2))
% subplot(2,2,4)
% imshow(raw_image(:,:,3))
%%
[xsize, ysize] = size(raw_image(:,:,1));
num_pixels = numel(raw_image(:,:,1));
disp('calculate optical density tuple ...');
opticalDensity = rgb2od(raw_image); % convert image to optical density values
filterOD = 0.15; % threshold low stain OD
opticalDensity(opticalDensity <= filterOD) = 0; % threshold low stain OD
opticalDensity(isnan(opticalDensity)) = 0; % omit NAN
opticalDensity(isinf(opticalDensity)) = 0; % omit infinity
%%
disp('manually select purple stain vector...');
figure; imshow(raw_image);
rect = getrect; rect = floor(rect);
purple_image_crop = imcrop(raw_image,rect);
purple_od = rgb2od(purple_image_crop); % convert to OD
purple_manual = mean(purple_od,2);
%
disp('manually select pink stain vector...');
rect = getrect; rect = floor(rect);
pink_image_crop = imcrop(raw_image,rect);
pink_od = rgb2od(pink_image_crop); % convert to OD
pink_manual = mean(pink_od,2);
%
disp('deconvole image by manually selected stain vectors...')
stain_mat_man = [purple_manual pink_manual];
saturation_mat_man = pinv(stain_mat_man)*opticalDensity;
purple_stain_rgb = stainvec2rgb(stain_mat_man(:,1),saturation_mat_man(1,:),xsize,ysize);
pink_stain_rgb = stainvec2rgb(stain_mat_man(:,2),saturation_mat_man(2,:),xsize,ysize);
remain_rgb_man = raw_image - purple_stain_rgb - pink_stain_rgb;
    fig_manual = figure;
    subplot(2,2,1)
    imshow(raw_image)
    subplot(2,2,2)
    imshow(purple_stain_rgb)
    subplot(2,2,3)
    imshow(pink_stain_rgb)
    subplot(2,2,4)
    imshow(remain_rgb_man)
%%
disp('calculate SVD on the optical density tuple ...')
[U, S, V] = svd(opticalDensity,'econ'); % calculate SVD from OD
% project onto the first 2 svd directions normalize to unit length
% cosine of angles to the first svd directions is first row of
% norm_project2svds, sine is the second row
%norm_project2svds = normc(S(1:2,1:2)*V(:,1:2)');
project2svds = S(1:2,1:2)*V(:,1:2)';
% cos_angles = norm_project2svds(1,:); sin_angles = norm_project2svds(2,:);
% disp('histogram of cosine and sine of angles to 1st svd ...')
% %
% figure;hist(cos_angles,101); title('cosine of angles');
% [n,s] = hist(cos_angles,101);
% figure;semilogy(s,n/sum(n))
% ax = axis;
% axis([-1 1 0 1])
% %
% figure;hist(sin_angles,101); title('sine of angles');
% [n,s] = hist(sin_angles,101);
% figure;semilogy(s,n/sum(n)); title('sine of angles in semilogy');
% ax = axis;
% axis([-1 1 0 1])
% % calculate tangents of the angles
% tan_angles = norm_project2svds(2,:)./norm_project2svds(1,:);
% figure;hist(tan_angles,101)
% [n,s] = hist(tan_angles,101);title('tan of angles');
% figure;semilogy(s,n/sum(n));title('tan of angles in semilogy');
% ax = axis;
% axis([-1 1 0 1])
% % calculate angles from tangents
% angles = atan(tan_angles);
% figure;hist(angles,101);title('distribution of angles');
% [n,s] = hist(angles,101);
% figure;semilogy(s,n/sum(n))
% 
% % deconvole the image into 2 colors
% % extreme values of cos_angles 
% extremeCutoff = 1;
% extreme_values = prctile(angles,[extremeCutoff 100-extremeCutoff]);
% indx_min = abs(angles - extreme_values(1)) < 1e-5;
% indx_max = abs(angles - extreme_values(2)) < 1e-5;
% % OD = VS where V=stainVectors, S = saturationMat%stain_mat = normc([mean(opticalDensity(:,indx_min),2) mean(opticalDensity(:,indx_max),2)]);
% stain_mat = normc([median(opticalDensity(:,indx_min),2) median(opticalDensity(:,indx_max),2)]);
% saturation_mat = pinv(stain_mat)*opticalDensity;
% 
% min_stain_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
% max_stain_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
% remain_rgb = raw_image - min_stain_rgb - max_stain_rgb;
%     figure;
%     subplot(2,2,1)
%     imshow(raw_image)
%     subplot(2,2,2)
%     imshow(min_stain_rgb)
%     subplot(2,2,3)
%     imshow(max_stain_rgb)
%     subplot(2,2,4)
%     imshow(remain_rgb)
end

