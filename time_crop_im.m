im_dir = 'D:\Documents\Tiles_Norm\time_test_2';
im_list = dir(fullfile(im_dir,'original_images','*.tif'));
im_list = {im_list.name}';
im_sizes = [512, 2048, 5000, 7500, 10000, 15000];

for i = 1:length(im_list)
   im_name = im_list{i}(1:end-4);
   im = imread(fullfile(im_dir,'original_images',im_list{i}));
   figure; imshow(im);
   for s = 5:length(im_sizes)
       [x,y] = getpts;
       cropped_im = imcrop(im,[x,y,im_sizes(s) - 1,im_sizes(s) - 1]);
       imwrite(cropped_im, fullfile(im_dir,'cropped_im',[im_name '_' num2str(im_sizes(s)) '.tif']));
   end
   close all;
end