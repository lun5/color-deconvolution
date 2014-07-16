% plot_manual_results;
% called by deconvolutionManual
split_string = regexp(imname,'\.','split');
savename = fullfile(resultdir,split_string{1});
h=figure;
subplot(2,2,1)
imshow(raw_image)
subplot(2,2,2)
imshow(purple_stain_rgb)
subplot(2,2,3)
imshow(pink_stain_rgb)
subplot(2,2,4)
imshow(remain_rgb_man)
print(h,'-dtiff', [savename '_man_deconvolution.tiff']);
