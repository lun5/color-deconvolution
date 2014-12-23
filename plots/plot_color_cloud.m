% plot_color_cloud
% called by deconvolve_image_sequence
savename = [resultdir filesep imagepaths{i}];
RGB = raw2rgb(raw_image)./255;step = 200;
%disp('compute projections of purple and pink vectors on the svd plane...');
purple_proj = purple_stain_man - purple_stain_man'*U(:,3)*U(:,3);
pink_proj = pink_stain_man - pink_stain_man'*U(:,3)*U(:,3);
%disp('compute projections of automated stain vectors on the svd plane
proj_stain_mat = stain_mat_auto - repmat((stain_mat_auto'*U(:,3))',3,1).*repmat(U(:,3),1,2);
%disp('pixel color distribution on the plane of the first 2 svds');
h = figure;scatter(project2svds(1,1:step:end), project2svds(2,1:step:end),20,RGB(:,1:step:end)', 'filled')
hold on
h1 = plot(purple_proj'*U(:,1),purple_proj'*U(:,2),'bs','MarkerSize',10,'MarkerFaceColor',...
    'b','MarkerEdgeColor','k','LineWidth',3);
    %mean(raw2rgb(purple_image_crop)./255,2),'MarkerEdgeColor','k','LineWidth',3);
h2 = plot(pink_proj'*U(:,1),pink_proj'*U(:,2),'ro','MarkerSize',10,'MarkerFaceColor',...
    'b','MarkerEdgeColor','k','LineWidth',3);
    %mean(raw2rgb(pink_image_crop)./255,2),'MarkerEdgeColor','k','LineWidth',3);
h3 = plot(proj_stain_mat(:,1)'*U(:,1),proj_stain_mat(:,1)'*U(:,2),'ks','MarkerSize',10,'MarkerFaceColor',...
    'w','MarkerEdgeColor','k','LineWidth',3);
h4 = plot(proj_stain_mat(:,2)'*U(:,1),proj_stain_mat(:,2)'*U(:,2),'ko','MarkerSize',10,'MarkerFaceColor',...
    'w','MarkerEdgeColor','k','LineWidth',3);
hold off
legend([h1, h2, h3, h4],'manual purple stain', 'manual pink stain','automated min stain','automated max stain');
xlabel('first svd'); ylabel('second svd');
title(['pixel colors for ' imagepaths{i}]);
ax = axis;
line([0 0], [ax(3) ax(4)]);line([ax(1) ax(2)],[0 0]);
%print(h,'-dtiff', [savename '_color_cloud.tiff']);
%close all;
