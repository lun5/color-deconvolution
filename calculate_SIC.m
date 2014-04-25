% calculate_SIC
% Thesis by Elmaraghi 2003

% sample pixels from the raw image
% close all;
% raw_image = imread([datadir filesep imname]);
% rgb_image_whole = raw2rgb(raw_image);
% [xsize, ysize] = size(raw_image(:,:,1));
% nsamples = 10000;
% indx_pixel = randperm(xsize*ysize,nsamples);
% rgb_image = rgb_image_whole(:,indx_pixel);
% rgb_image = raw2rgb(raw_image);
step = 1000; % for plotting 

% set parameters
% prompt = 'Please enter the mu value ';
% mu_s = input(prompt);
% while isempty(imname)
%     mu_s = input(prompt);
% end
% 
% prompt = 'Please enter the sigma value ';
% sigma_s = input(prompt);
% while isempty(imname)
%     sigma_s = input(prompt);
% end

options = struct('Normalize','on');
%options = struct('Normalize','off');
mu_s = 4; sigma_s = 2;
% obtain the color opponency matrix
workdir = '/Users/lun5/Research/color_deconvolution'; 
load([workdir filesep 'training_pink_purple.mat'],'training_data')
[U,D,V] = svd(training_data);
rotation_matrix = []; %[-U(:,1) U(:,2:3)];
%rotation_matrix = [-U(:,1) U(:,2:3)]';
%rotation_matrix = [-U(:,1) U(:,2:3)];
% calculate the sic coordinates
sic_image = rgb2sic( rgb_image, mu_s, sigma_s, rotation_matrix, options); 

% This is stupid move by me
[ purple_manual_rgb,pink_manual_rgb, stain_mat_man] = deconvolutionManual( imname,datadir,resultdir);% options ); 
%
purple_stain_man_sic = rgb2sic(purple_manual_rgb,mu_s, sigma_s, rotation_matrix, options); 
pink_stain_man_sic = rgb2sic(pink_manual_rgb,mu_s, sigma_s, rotation_matrix, options); 

%% Plot the SIC space
h=figure; scatter(sic_image(1,1:step:end),sic_image(2,1:step:end),20,rgb_image(:,1:step:end)'./255,'filled');
hold on
h1=plot(pink_stain_man_sic(1),pink_stain_man_sic(2),'bs','MarkerSize',10,...%'MarkerFaceColor',pink_manual_rgb./255,
    'MarkerEdgeColor','k','LineWidth',3);
h2=plot(purple_stain_man_sic(1),purple_stain_man_sic(2),'bo','MarkerSize',10,...%'MarkerFaceColor',purple_manual_rgb./255,
    'MarkerEdgeColor','r','LineWidth',3);
legend([h1, h2],'manual pink stain', 'manual purple stain');
xlabel('s1','FontSize',15);ylabel('s2','FontSize',15);
line([-1 1],[0 0])
line([0 0],[-1 1])
hold off
title('SIC 2D representation of hue','FontSize',15);
axis([-1 1 -1 1])
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [resultdir filesep imname '_sic_rep.tiff']);
end
%% distribution of coordinates s1
h=figure;hist(sic_image(1,1:step:end),30);
%b1 = bar(hist(sic_image(1,:),30) ./ sum(sic_image(1,:)));
ax = axis;
hold on
h1 = line([pink_stain_man_sic(1) pink_stain_man_sic(1)],[ax(3) ax(4)],'LineWidth',3, 'Color',pink_manual_rgb./255);
h2 = line([purple_stain_man_sic(1) purple_stain_man_sic(1)],[ax(3) ax(4)],'LineWidth',3, 'Color',purple_manual_rgb./255);
legend([h1, h2],'manual pink stain', 'manual purple stain');
hold off
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','w')
title('Distribution of s1','FontSize',15);
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [resultdir filesep imname '_s1_dist.tiff']);
end
% save this to the result

%% distribution of coordinates s2
h=figure;hist(sic_image(2,1:step:end),30);
ax = axis;
hold on
h1 = line([pink_stain_man_sic(2) pink_stain_man_sic(2)],[ax(3) ax(4)],'LineWidth',3, 'Color',pink_manual_rgb./255);
h2 = line([purple_stain_man_sic(2) purple_stain_man_sic(2)],[ax(3) ax(4)],'LineWidth',3, 'Color',purple_manual_rgb./255);
legend([h1, h2],'manual pink stain', 'manual purple stain');
hold off
title('Distribution of s2','FontSize',15);
set(gca,'FontSize',15);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','w')
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [resultdir filesep imname '_s2_sic.tiff']);
end

%% Contribution of the first element
% histogram of the first component
rotated_coordinates = rotation_matrix*rgb_image;
figure;hist(rotated_coordinates(1,:),100);
set(gca,'FontSize',15);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .8 .8],'EdgeColor','w')
xlabel('Distribution of the 1st svd','FontSize',15);

% make an image where each pixel value is the corresponding first coeff
firstsvd_image = reshape(rotated_coordinates(1,:),[xsize, ysize]);
figure;imshow(firstsvd_image,[]);
set(gca,'FontSize',15);
title('Pixel is first coeff','FontSize',15);
%remove the contribution of the first SVD vector from each pixel and 
%show how the residual image looks like? 

% OD = VS where V=stainVectors, S = saturationMat
saturation_mat = pinv(rotation_matrix(:,1))*opticalDensity;
%saturation_mat = pinv(rotation_matrix(1,:))'*opticalDensity;
firstsvd_rgb = stainvec2rgb(rotation_matrix(:,1),saturation_mat,xsize,ysize);
remain_1stsvd_rgb = raw_image - firstsvd_rgb;
figure;imshow(remain_1stsvd_rgb,[])
title('Remove contribution of 1st svd','FontSize',15);

% 
%% show an image of the saturated amplitude
sat_levels = sqrt(sic_image(1,:).^2 + sic_image(2,:).^2);
sat_image = reshape(sat_levels ,[xsize, ysize]);
figure; imshow(sat_image,[]);
title('Saturated level','FontSize',15);

%% show an image of the hue angle
% bin the hue angle into say 8 categories and give them separate colors

% testing effect of different values of mu and sigma
% sigma_s_vector = 1:1:10;
% mu_s_vector = 1:1:10;
% for i = 1:length(mu_s_vector) %length(sigma_s_vector)
%     sigma_s = sigma_s_vector(i);
%     %mu_s = mu_s_vector(i);
%     sic_image = rgb2sic( rgb_image, mu_s, sigma_s, [], options); 
%     subplot(2,length(mu_s_vector)/2,i);
%     scatter(sic_image(1,:),sic_image(2,:),20,rgb_image'./255,'filled');
%     hold on
%     xlabel('s1');ylabel('s2');
%     line([-1 1],[0 0]);line([0 0],[-1 1])
%     hold off
%     title(['mu_s = ' num2str(mu_s) ' sigma_s = ' num2str(sigma_s)]); axis([-1 1 -1 1]);
% end
