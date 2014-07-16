% called by colorDeconvolutionSVD
%plot_svd_results
RGB = raw2rgb(raw_image)./255; step = 200;
savename = [resultdir filesep imname];
% distribution of angles normal and log scale
h = figure;
subplot(2,1,1);hist(angles_1st_svd,101);title('distribution of angles');
[n,s] = hist(angles_1st_svd,101);
subplot(2,1,2);semilogy(s,n/sum(n))
%print(h,'-dtiff', [savename '_auto_angleDistribution.tiff']);

% deconvolved images
h = figure;
subplot(2,2,1)
imshow(raw_image)
subplot(2,2,2)
imshow(stain1_rgb)
subplot(2,2,3)
imshow(stain2_rgb)
subplot(2,2,4)
imshow(remain_rgb)
%print(h,'-dtiff', [savename '_auto_deconv.tiff']);

% plot of saturation amount
h = figure;
scatter(saturation_mat(1,1:step:end),saturation_mat(2,1:step:end),20,RGB(:,1:step:end)','filled');
ax=axis; hold on; plot(min(ax):0.1:max(ax),min(ax):0.1:max(ax));
line([0 0],[min(ax) max(ax)]);line([min(ax) max(ax)],[0 0]);
title(['Saturation level of ' imname]);
%print(h,'-dtiff', [savename '_auto_saturation.tiff']);
