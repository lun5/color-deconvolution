function optical_density_clouds( raw_image, step )
%optical_density_clouds generates scatter plot of pixels in optical density
%plane defined by the first 2 svds
RGB = raw2rgb(raw_image)./255; % RGB color for the image
disp('calculate optical density tuple ...');
opticalDensity = rgb2od(raw_image); % convert image to optical density values
filterOD = 0.15; % threshold low stain OD
opticalDensity(opticalDensity <= filterOD) = 0; % threshold low stain OD
opticalDensity(isnan(opticalDensity)) = 0; % omit NAN
opticalDensity(isinf(opticalDensity)) = 0; % omit infinity
disp('calculate SVD on the optical density tuple ...')
[U, S, V] = svd(opticalDensity,'econ');
project2svds = S(1:2,1:2)*V(:,1:2)';
%%
disp('manually select purple stain vector...');
figure; imshow(raw_image);zoom(2)
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
disp('compute projections of purple and pink vectors on the svd plane...');
purple_proj = purple_manual - purple_manual'*U(:,3)*U(:,3);
pink_proj = pink_manual - pink_manual'*U(:,3)*U(:,3);

%%
disp('calculate automated stain vectors')
norm_project2svds = normc(project2svds);
tan_angles = norm_project2svds(2,:)./norm_project2svds(1,:);
angles_1st_svd = atan(tan_angles);
extremeCutoff = 1; 
%angles_1st_svd = atan2(project2svds(2,:), project2svds(1,:));
disp('distribution of angles to the 1st svd')
figure;hist(angles_1st_svd,101);title('distribution of angles');
[n,s] = hist(angles_1st_svd,101);
figure;semilogy(s,n/sum(n))

% extreme values of cos_angles 
extreme_values = prctile(angles_1st_svd,[extremeCutoff 100-extremeCutoff]);
indx_min = abs(angles_1st_svd - extreme_values(1)) < 1e-5;
indx_max = abs(angles_1st_svd - extreme_values(2)) < 1e-5;
% OD = VS where V=stainVectors, S = saturationMat
stain_mat = normc([median(opticalDensity(:,indx_min),2) median(opticalDensity(:,indx_max),2)]);
disp('projections of automated vectors on the svd plane')
proj_stain_mat = [median(project2svds(:,indx_min),2) median(project2svds(:,indx_max),2)];
%%
disp('pixel color distribution on the plane of the first 2 svds');
figure;scatter(project2svds(1,1:step:end), project2svds(2,1:step:end),20,RGB(:,1:step:end)', 'filled')
hold on
h1 = plot(purple_proj'*U(:,1),purple_proj'*U(:,2),'bs','MarkerSize',10,'MarkerFaceColor',...
    mean(raw2rgb(purple_image_crop)./255,2),'MarkerEdgeColor','k','LineWidth',3);
h2 = plot(pink_proj'*U(:,1),pink_proj'*U(:,2),'ro','MarkerSize',10,'MarkerFaceColor',...
    mean(raw2rgb(pink_image_crop)./255,2),'MarkerEdgeColor','k','LineWidth',3);
h3 = plot(proj_stain_mat(1,1),proj_stain_mat(2,1),'ks','MarkerSize',10,'MarkerFaceColor',...
    'w','MarkerEdgeColor','k','LineWidth',3);
h4 = plot(proj_stain_mat(1,2),proj_stain_mat(2,2),'ko','MarkerSize',10,'MarkerFaceColor',...
    'w','MarkerEdgeColor','k','LineWidth',3);
hold off
legend([h1, h2, h3, h4],'manual purple stain', 'manual pink stain','automated min stain','automated max stain');
xlabel('first svd'); ylabel('second svd');
title('pixel colors for slides 1 #1');
ax = axis;
%axis([-3 .5 -1 1]);
line([0 0], [ax(3) ax(4)]);line([ax(1) ax(2)],[0 0]);

end

% angles_1st_svd = atan2(project2svds(2,:), project2svds(1,:));

%Purple 	128 	0 	128 
%Pink 	255 	192 	203