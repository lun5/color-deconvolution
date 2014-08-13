% March 21, 2014
% This script

%Purple 	128 	0 	128 
%Pink 	255 	192 	203
close all;
standard_purple_rgb = [255;1;128];
standard_pink_rgb = [255;192;203];
standard_purple_rgb = [255;1;255];
standard_pink_rgb = [1;255;1];
% convert to OD space
standard_purple_od = -log10(standard_purple_rgb./255);
standard_pink_od = -log10(standard_pink_rgb./255);
% read in an image
workdir = '/Users/lun5/Research/deconvolution'; 
datadir = fullfile(workdir, '20x_images','TIFF');
fname = fullfile(datadir,'Slide#2-3.TIF');

% manually deconvolved the image
[ purple_stain_man,pink_stain_man] = deconvolutionManual( imname,datadir,resultdir);% options ); % why do I need workdir hmm

% deconvolved using the two standard vectors
raw_image = imread(fname);
[xsize, ysize] = size(raw_image(:,:,1));
% convert to OD space
calculate_optical_density
% project the image onto the 2 stain vectors
standard_purple_proj = normc(standard_purple_od)'*opticalDensity;
standard_pink_proj = normc(standard_pink_od)'*opticalDensity;
RGB = raw2rgb(raw_image)./255; % RGB color for the image
%standard_purple_proj = normc(standard_purple_rgb)'*RGB;
%standard_pink_proj = normc(standard_pink_rgb)'*RGB;
figure;hist(standard_purple_proj,100);
figure;hist(standard_pink_proj,100);
% plot the optical density space and also the two vectors
step = 200;
figure;
scatter(standard_purple_proj(1:step:end), standard_pink_proj(1:step:end),20,RGB(:,1:step:end)', 'filled');
hold on
h1 = plot(purple_stain_man'*normc(standard_purple_od),purple_stain_man'*normc(standard_pink_od),...
    'bs','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',3);
    %mean(raw2rgb(purple_image_crop)./255,2),'MarkerEdgeColor','k','LineWidth',3);
h2 = plot(pink_stain_man'*normc(standard_purple_od),pink_stain_man'*normc(standard_pink_od),...
    'ro','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',3);
hold off
legend([h1, h2],'manual purple stain', 'manual pink stain');
xlabel('[1 0 1]'); ylabel('[0 1 0]');
%xlabel('standard purple'); ylabel('standard pink');
% pick out the most similar
% find out which vector has the largest, 90, 95, etc. 
% this doesn't work because the most pink is also the most purple
%%
% cutoff = 90;
% most_purple_val = prctile(standard_purple_proj,cutoff);
% most_pink_val =  prctile(standard_pink_proj,cutoff);
% indx_purple = abs(standard_purple_proj - most_purple_val) < 1e-5;
% indx_pink = abs(standard_pink_proj - most_pink_val) < 1e-5;
% stain_mat = [median(opticalDensity(:,indx_purple),2) median(opticalDensity(:,indx_pink),2)];
% %stain_mat = [standard_purple_od standard_pink_od];
% saturation_mat = pinv(stain_mat)*opticalDensity;
% purple_stain_rgb = stainvec2rgb(stain_mat(:,1),saturation_mat(1,:),xsize,ysize);
% pink_stain_rgb = stainvec2rgb(stain_mat(:,2),saturation_mat(2,:),xsize,ysize);
% remain_rgb = raw_image - purple_stain_rgb - pink_stain_rgb;
% % look at the pixels that are most similar to standard purple vector
% purple_region_od = opticalDensity(:,indx_purple);
% purple_region_rgb = od2rgb(purple_region_od,size(purple_region_od,2),1);
% %purple_region_rgb = reshape(purple_region_rgb(1:440,:,:),44,10,3);
% figure;imshow(purple_region_rgb);
% % look at pixels that are similar to standard pink vector
% pink_region_od = opticalDensity(:,indx_pink);
% pink_region_rgb = od2rgb(pink_region_od,size(pink_region_od,2),1);
% %pink_region_rgb = reshape(pink_region_rgb(1:66,:),6,11,3);
% figure;imshow(pink_region_rgb);
% deconvolved images
h = figure;
subplot(2,2,1)
imshow(raw_image)
subplot(2,2,2)
imshow(purple_stain_rgb)
subplot(2,2,3)
imshow(pink_stain_rgb)
subplot(2,2,4)
imshow(remain_rgb)

% Try vector [1 0 1] and [0 1 0]
% then try nnmf

