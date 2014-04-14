% calculate_SIC
% Thesis by Elmaraghi 2003
rotation_matrix = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
    1/sqrt(6) 1/sqrt(6) -2/sqrt(6); ...
    -1/sqrt(2) 1/sqrt(2) 0];
rotated_coordinates = rotation_matrix*rgb_image;

% mu_s = 4; sigma_s = 4; % transition control parameters
% F = @(x, mu, sigma) 1-exp(- max(x - mu,0)/2/sigma^2);
% Fmatrix = F(sqrt(sum(rotated_coordinates(2:3,:).^2,1)),mu_s,sigma_s); 
% sic_image = repmat(Fmatrix,2,1).* rotated_coordinates(2:3,:)./...
%     repmat(1 + sqrt(sum(rotated_coordinates(2:3,:).^2,1)),2,1);

% do not normalize
sic_image = rotated_coordinates(2:3,:)./...
     repmat(1 + sqrt(sum(rotated_coordinates(2:3,:).^2,1)),2,1);

purple_stain_man_sic = rgb2sic(purple_manual_rgb);
pink_stain_man_sic = rgb2sic(pink_manual_rgb);
figure; scatter(sic_image(1,:),sic_image(2,:),20,rgb_image'./255,'filled');
hold on
h1=plot(pink_stain_man_sic(1),pink_stain_man_sic(2),'bs','MarkerSize',20,'MarkerFaceColor',...
    pink_manual_rgb./255,'MarkerEdgeColor','k','LineWidth',3);

h2=plot(purple_stain_man_sic(1),purple_stain_man_sic(2),'bo','MarkerSize',20,'MarkerFaceColor',...
    purple_manual_rgb./255,'MarkerEdgeColor','k','LineWidth',3);
hold off
xlabel('s1');
ylabel('s2');
legend([h1, h2],'manual pink stain', 'manual purple stain');
title('SIC 2D representation of hue');
axis([-1 1 -1 1])
