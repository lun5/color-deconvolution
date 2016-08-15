% calculate_SIC
% Inspired by Thesis by Elmaraghi 2003

% sample pixels from the raw image
% close all;
% raw_image = imread([datadir filesep imname]);
% rgb_image_whole = raw2rgb(raw_image);
% [xsize, ysize] = size(raw_image(:,:,1));
% nsamples = 10000;
% indx_pixel = randperm(xsize*ysize,nsamples);
% rgb_image = rgb_image_whole(:,indx_pixel);
% rgb_image = raw2rgb(raw_image);

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
%%
step = 1000; % for plotting 
options = struct('Normalize','on');
%options = struct('Normalize','off');
mu_s = 4; sigma_s = 2;
% obtain the color opponency matrix
% workdir = '/Users/lun5/Research/color_deconvolution'; 
% training1 = load([workdir filesep 'training_pink_purple.mat'],'training_data');
% training2 = load([workdir filesep 'results' filesep '140428' filesep  'training_pink_purple.mat'],'training_data');
% %training3 = load([workdir filesep 'results' filesep '140507' filesep  'training_pink_purple.mat'],'training_data');
% training_data = [training1.training_data training2.training_data];
% %training_data = training3.training_data(:,1:4000);
% [U,D,V] = svd(training_data,0);
% rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one
% calculate the sic coordinates
sic_image = rgb2sic( rgb_image, mu_s, sigma_s, rotation_matrix, options); 
%% Eliminate chemical saturation (black) and optical saturation (white)
rotated_coordinates = rotation_matrix*rgb_image;
first_component = uint8(rotated_coordinates(1,:));
indx_chemical_sat = first_component < 5;
indx_optical_sat = first_component > 250;
indx_sat = indx_chemical_sat + indx_optical_sat;
%%
% manual deconvolution
if strcmpi(plotflag,'on')
[ purple_manual_rgb,pink_manual_rgb, stain_mat_man] = deconvolutionManual( imname,datadir,resultdir);% options ); 
%
purple_stain_man_sic = rgb2sic(purple_manual_rgb,mu_s, sigma_s, rotation_matrix, options); 
pink_stain_man_sic = rgb2sic(pink_manual_rgb,mu_s, sigma_s, rotation_matrix, options); 

%% get the image name for plotting
split_string = regexp(imname,'\.','split');
savename = fullfile(resultdir,split_string{1});
%% Plot the color opponency space (cite SIC)
%if strcmpi(plotflag,'on')
h=figure; scatter(sic_image(1,1:step:end),sic_image(2,1:step:end),20,rgb_image(:,1:step:end)'./255,'filled');
hold on
h1=plot(pink_stain_man_sic(1),pink_stain_man_sic(2),'bs','MarkerSize',10,...%'MarkerFaceColor',pink_manual_rgb./255,
    'MarkerEdgeColor','k','LineWidth',3);
h2=plot(purple_stain_man_sic(1),purple_stain_man_sic(2),'bo','MarkerSize',10,...%'MarkerFaceColor',purple_manual_rgb./255,
    'MarkerEdgeColor','r','LineWidth',3);
%legend([h1, h2],'manual pink stain', 'manual purple stain');
xlabel('s1','FontSize',15);ylabel('s2','FontSize',15);
line([-1 1],[0 0],'Color','k')
line([0 0],[-1 1],'Color','k')
hold off
title('Color opponency 2D representation of hue','FontSize',15);
axis([-1 1 -1 1])
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_sic_rep.tiff']);
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
h3 = findobj(gca,'Type','patch');
set(h3,'FaceColor',[.8 .8 .8],'EdgeColor','w')
title('Distribution of s1','FontSize',15);
set(gca,'FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_s1_dist.tiff']);
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
h3 = findobj(gca,'Type','patch');
set(h3,'FaceColor',[.8 .8 .8],'EdgeColor','w')
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_s2_sic.tiff']);
end

%% Contribution of the first element
% histogram of the first component
rotated_coordinates = rotation_matrix*rgb_image;
first_component = uint8(rotated_coordinates(1,:));
h = figure;hist(double(first_component),50);
set(gca,'FontSize',15);
h3 = findobj(gca,'Type','patch');
set(h3,'FaceColor',[.8 .8 .8],'EdgeColor','w')
xlabel('Distribution of the 1st svd','FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_dist_first_comp.tiff']);
end

% make an image where each pixel value is the corresponding first coeff
firstsvd_image = reshape(rotated_coordinates(1,:),[xsize, ysize]);
h = figure;imshow(firstsvd_image,[]);
set(gca,'FontSize',15);
title('Pixel is first coeff','FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_first_svd_image.tiff']);
end

%remove the contribution of the first SVD vector from each pixel and 
%show how the residual image looks like? 
residual = rgb_image - rotation_matrix(1,:)'*(rotation_matrix(1,:)*rgb_image);
RGB = abs(residual);
r = reshape(RGB(1,:),[xsize, ysize]);
g = reshape(RGB(2,:),[xsize, ysize]);
b = reshape(RGB(3,:),[xsize, ysize]);
residual_image = uint8(cat(3,r,g,b));
%imtool(residual_image);
% OD = VS where V=stainVectors, S = saturationMat
% saturation_mat = pinv(rotation_matrix(1,:))'*opticalDensity;
% firstsvd_rgb = stainvec2rgb(rotation_matrix(1,:)',saturation_mat,xsize,ysize);
% remain_1stsvd_rgb = raw_image - firstsvd_rgb;
% figure;imshow(remain_1stsvd_rgb,[])
h = figure; imshow(residual_image);
title('Remove contribution of 1st svd','FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_remove_brightness.tiff']);
end
% 
%% show an image of the saturated amplitude
sat_levels = sqrt(rotated_coordinates(2,:).^2 + rotated_coordinates(3,:).^2);
sat_image = reshape(sat_levels ,[xsize, ysize]);
h = figure; imshow(sat_image,[]);
title('Saturated level','FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_saturated_values.tiff']);
end
%% show an image of the hue angle
% bin the hue angle into say 8 categories and give them separate colors
hue_vector = angle(rotated_coordinates(2,:) + rotated_coordinates(3,:)*1i);
pink_stain_man_rotated = rotation_matrix*pink_manual_rgb;
pink_stain_hue = angle(pink_stain_man_rotated(2)+pink_stain_man_rotated(3)*1i);
purple_stain_man_rotated = rotation_matrix*purple_manual_rgb;
purple_stain_hue = angle(purple_stain_man_rotated(2)+purple_stain_man_rotated(3)*1i);
%
h = figure;hist(hue_vector,8);
ax = axis;
hold on
h1 = line([pink_stain_hue pink_stain_hue],[ax(3) ax(4)],'LineWidth',3, 'Color',pink_manual_rgb./255);
h2 = line([purple_stain_hue purple_stain_hue],[ax(3) ax(4)],'LineWidth',3, 'Color',purple_manual_rgb./255);
legend([h1, h2],'manual pink stain', 'manual purple stain');
hold off
title('Histogram of Hues','FontSize',15);
if strcmpi(saveflag,'on')
    print(h,'-dtiff', [savename '_hue_distribution.tiff']);
end
end
% collect more training data. Include the whole patch, not just the average
% of the rgb vectors
% eliminate chemical and optical saturation vectors by projecting them onto
% the first component (brightness). Eliminate the two low and two high ones

%%
%h = findobj(gca,'Type','patch');
%cm=colormap(jet(8));
%set(h,'FaceColor',cm,'EdgeColor','w')
%set(gca,'color',cm);
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
