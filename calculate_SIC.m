% calculate_SIC
% Thesis by Elmaraghi 2003

raw_image = imread([datadir filesep imname]);
rgb_image_whole = raw2rgb(raw_image);
[xsize, ysize] = size(raw_image(:,:,1));
nsamples = 10000;
indx_pixel = randperm(xsize*ysize,nsamples);
rgb_image = rgb_image_whole(:,indx_pixel);

options = struct('Normalize','on');

% testing effect of different values of mu and sigma
mu_s = 4;
sigma_s = 4;
sigma_s_vector = 1:1:10;
for i = 1:length(sigma_s_vector)
    sigma_s = sigma_s_vector(i);
    sic_image = rgb2sic( rgb_image, mu_s, sigma_s, options); 
    figure; scatter(sic_image(1,:),sic_image(2,:),20,rgb_image'./255,'filled');
    hold on
    xlabel('s1');ylabel('s2');
    line([-1 1],[0 0]);line([0 0],[-1 1])
    hold off
    title('SIC 2D representation of hue'); axis([-1 1 -1 1]);
end
       
% purple_stain_man_sic = rgb2sic(purple_manual_rgb);
% pink_stain_man_sic = rgb2sic(pink_manual_rgb);
% figure; scatter(sic_image(1,:),sic_image(2,:),20,rgb_image'./255,'filled');
% hold on
% h1=plot(pink_stain_man_sic(1),pink_stain_man_sic(2),'bs','MarkerSize',20,'MarkerFaceColor',...
%     pink_manual_rgb./255,'MarkerEdgeColor','k','LineWidth',3);
% 
% h2=plot(purple_stain_man_sic(1),purple_stain_man_sic(2),'bo','MarkerSize',20,'MarkerFaceColor',...
%     purple_manual_rgb./255,'MarkerEdgeColor','k','LineWidth',3);
% legend([h1, h2],'manual pink stain', 'manual purple stain');
% xlabel('s1');
% ylabel('s2');
% line([-1 1],[0 0])
% line([0 0],[-1 1])
% hold off
% title('SIC 2D representation of hue');
% axis([-1 1 -1 1])
