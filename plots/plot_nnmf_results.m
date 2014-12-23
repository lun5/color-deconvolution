% plot_nnmf_results
% called from deconvolution NNMF
% deconvolved images
RGB = raw2rgb(raw_image)./255; step = 200;
savename = [resultdir filesep imname];
h = figure;
subplot(2,2,1)
imshow(raw_image)
subplot(2,2,2)
imshow(stain1_rgb)
subplot(2,2,3)
imshow(stain2_rgb)
subplot(2,2,4)
imshow(remain_rgb)
print(h,'-dtiff', [savename '_nnmf_deconvolution.tiff']);

% plot of saturation amount
h = figure;
scatter(saturation_mat_nnmf(1,1:step:end),saturation_mat_nnmf(2,1:step:end),20,RGB(:,1:step:end)','filled');
ax=axis; hold on; plot(min(ax):0.1:max(ax),min(ax):0.1:max(ax));
line([0 0],[min(ax) max(ax)]);line([min(ax) max(ax)],[0 0]);
title(['Saturation level of ' imname]);
