%% test script for normalization
% load('source_image_data_all.mat');
% load('wsi_test_normalization.mat');
%load('wsi_new_mask.mat');
%source_im = tmpim3;
%load('newtargetWithPinkPurpleMask_highres');source_im = tmpim2;
%load('newtargetWithPinkPurpleMask_lowres');source_im = tmpim2;
%pink_purple_mask = pink_purple_mask;

%load '../normalizationImageWithIssue.mat';
%rect = [1.6139    2.3686    2.0752    1.5954]*1e3;
%rect = [5428, 2435, 1439, 1375];
%source_im = imcrop(tmpim3,rect);pink_purple_mask = imcrop(pink_purple_mask,rect);

%im_dir = 'C:\Users\luong_nguyen\Box Sync\ADH\segmentation trials\adh03-2a_bb';
%imname = 'adh03-2a_38595_39427_41186_42078';
%im_dir = 'C:\Users\luong_nguyen\Box Sync\ADH\segmentation trials\adh15-1a_bb';
%imname = 'adh15-1a_10933_11232_20928_21200';
%im_dir = 'C:\Users\luong_nguyen\Box Sync\ADH\segmentation trials\adh47-1a_bb';
%out_dir = 'C:\Users\luong_nguyen\Box Sync\ADH\segmentation trials\adh47-1a_bb_normalized';

%im_dir = 'Y:\MITOS\training_aperio\A03\frames\x20';
%out_dir = 'Y:\MITOS\training_aperio\A03\frames\x20\normalized';
im_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Images';
out_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\NormalizedImages\LuongNewTarget';
imlist = dir(fullfile(im_dir,'*.tif'));
imlist = {imlist.name}';

if ~exist(out_dir,'dir'); mkdir(out_dir);end;
for i = 1:length(imlist)
imname = imlist{i};
source_im = double(imread(fullfile(im_dir,imname)));
target_im = double(imread('1bhjqxecct.tif'));
pink_purple_mask = ones(size(source_im,1),size(source_im,2))>0;
% source_im = tmpim3;
% load 'source_image_crop2'
tic;
[ normalized_image] = target_features(source_im, target_im, pink_purple_mask);
imwrite(normalized_image,fullfile(out_dir,imname));

% [ source_eq_image] = col_normalization( source_im, mask_source);
toc
figure; 
subplot(1,3,1);imshow(source_im/255);title('Source');
subplot(1,3,2);imshow(target_im/255);title('Target');
% subplot(1,3,1);imshow(tmpim2./255);
subplot(1,3,3);imshow(normalized_image);title('Normalized Image');
% figure; imshow(normalized_image);
% imwrite(normalized_image,fullfile(tiles_dir,['new_eq_' source_im_name]));
end
