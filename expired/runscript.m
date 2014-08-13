% test run colorDeconvolution
clear all
%close all
clc
workdir = '/Users/lun5/Research/deconvolution';
cd(workdir)
imagedir = fullfile(workdir, '20x_images','TIFF');
%imagedir = fullfile(workdir,'0.7x_imagesP_TIFF','0.7x_images_TIFF');
fileNames = dir(fullfile(imagedir,'*.TIF'));
%fileNames = dir(fullfile(imagedir,'*.tif'));
fileNames = {fileNames.name}';

%%
% first image
% raw_image = imread([imagedir filesep fileNames{1}]);
% figure;
% subplot(2,2,1)
% imshow(raw_image)
% subplot(2,2,2)
% imshow(raw_image(:,:,1))
% subplot(2,2,3)
% imshow(raw_image(:,:,2))
% subplot(2,2,4)
% imshow(raw_image(:,:,3))

%%
% recreate figure 2
% num_pixels = numel(raw_image(:,:,1));
% indx = randsample(num_pixels,10000);
% RGB = raw2rgb(raw_image)./255; % R,G,B vector
% od = rgb2od(raw_image);
% color_vec = RGB(:,indx)'; 
% subplot(1,2,1); scatter3(RGB(1,indx), RGB(2,indx), RGB(3,indx), 2, color_vec);
% subplot(1,2,2); scatter3(od(1,indx), od(2,indx), od(3,indx), 2, color_vec);
% figure;
% scatter(od(1,indx), od(2,indx), 2, color_vec);
% xlabel('R'); ylabel('B'); title('OD space');
% figure;
% scatter(RGB(1,indx), RGB(2,indx), 2, color_vec);
% xlabel('R'); ylabel('B'); title('RGB space');

%%
% store the stain vector
%stain_vectors = [];
extremeCutoff = 1;
filterOD = 0.15;
for i = 5 %: length(fileNames)
    raw_image = imread([imagedir filesep fileNames{i}]);
    [ min_stain_rgb, max_stain_rgb, remain_rgb,stain_mat,~] = colorDeconvolution([imagedir filesep fileNames{i}], filterOD, extremeCutoff);
    title(['automated deconvolution of ' num2str(i)])
    [ purple_stain_rgb ,pink_stain_rgb,remain_rgb_man ] = deconvolution_man( raw_image, filterOD );
    title(['manual deconvolution of ' num2str(i)])
%     figure;
%     subplot(2,2,1)
%     imshow(raw_image)
%     subplot(2,2,2)
%     imshow(min_stain_rgb)
%     subplot(2,2,3)
%     imshow(max_stain_rgb)
%     subplot(2,2,4)
%     imshow(remain_rgb)
    %stain_vectors = [stain_vectors; stain_mat'];
end
% number 4, 5 not work so well
ind = 2:2:length(stain_vectors);
X = stain_vectors(ind,1);
Y = stain_vectors(ind,2);
Z = stain_vectors(ind,3);

subplot(2,2,1)
scatter(X,Y,50,stain_vectors(ind,:),'filled');
xlabel('R'); ylabel('G');
subplot(2,2,2)
scatter(X,Z,50,stain_vectors(ind,:),'filled');
xlabel('R'); ylabel('B');
subplot(2,2,3)
scatter(Y,Z,50,stain_vectors(ind,:),'filled');
xlabel('G'); ylabel('B');
subplot(2,2,4)
scatter3(X,Y,Z,50,stain_vectors(ind,:),'filled');
xlabel('R'); ylabel('G');zlabel('B')


%% test SVD
A = rand(3,1000);
figure; scatter(A(1,:), A(2,:));
[U,S,V] = svd(A,'econ');
B = normc(S*V');
figure; scatter(B(1,:), B(2,:)); 
hold on
plot([0 0], [-1 1],'-k')
plot([-1 1],[0 0], '-k')
hold off

C = S*V';
figure; scatter(C(1,:), C(2,:)); 
hold on
plot([0 0], [-1 1],'-k')
plot([-1 1],[0 0], '-k')
hold off

hist(norm_project2svds(2,:),100)
[n,s] = hist(norm_project2svds(2,:),51);
semilogy(s,n/sum(n))
ax = axis;
axis([-1 1 0 1])


