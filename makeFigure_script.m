%% Figure 2: H&E color space
%im_rgb = double(imread('D:\Documents\Tiles_Norm\Tiles_2k\95f7k8loesyevi.tif'))./255;
im_rgb = double(imread('D:\Documents\Tiles_Norm\Tiles_2k\2ale5ngryfnpo.tif'))./255;

rotation_matrix = load('rotation_matrix_tp10-867-1.mat','rotation_matrix');

nrows = size(im_rgb,1); ncols = size(im_rgb,2);
X = reshape(im_rgb,[nrows*ncols,3]);
%rotated_coordinates = rotation_matrix.rotation_matrix*X';
tmp = load('../stain_vectors_tissue_images_tp10-867-1.mat');
stains_vec = tmp.vahadane_stains;
stains_vec = cat(1,stains_vec{:});
[U,~,~] = svd(stains_vec',0);
rotmat =  [-U(:,1) U(:,2:3)]';

rotated_coordinates = rotmat*X';

theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));

im_theta = reshape(theta,[nrows, ncols]);
figure; imagesc(im_theta);
cmap = colormap(hsv); colorbar('southoutside');
axis equal; axis off; axis tight;set(gca,'FontSize',20);
cmap1 = cmap([8:end, 1:9],:);
colormap(cmap1);

numClusters = 3;
opts_mixture.noise = 1;
X_cart = [cos(theta); sin(theta)]';
%% Call the function
[ mu_hat_polar,~, kappa_hat,posterior_probs, prior_probs] =...
    moVM(X_cart,numClusters,opts_mixture);

x = -pi:0.1:pi;
%c = ['r','g','b'];
c = [ 128 0 128; 205 145 158; 0 0 0]./255;
figure;
%histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'BinWidth',0.2);
%histogram(theta,'Normalization','pdf','FaceColor',[0.8 0.8 0.8],'NumBins',50);
drawWheel(theta,50,[0.8 0.8 0.8]);
set(gca,'XTickLabel',[]) 
set(gca,'YTickLabel',[]) 
set(gca,'FontSize',20);

hold on;
for cl=1:numClusters
    yk = prior_probs(cl)*circ_vmpdf(x, mu_hat_polar(cl), kappa_hat(cl));
    %plot(x, yk,'Color',c(cl),'LineStyle','-','LineWidth',2); hold on;
    circ_line(x,yk,c(cl,:));
end

%         if opts_mixture.noise
%             yk = prior_probs(numClusters + opts_mixture.noise)./(2*pi);
%             %plot(x, yk,'Color','k','LineStyle','-','LineWidth',2); hold on;
%             circ_line(x,yk,'k');
%         end

hold off; %xlim([-pi pi]);
set(gcf,'color','white') % White background for the figure.
