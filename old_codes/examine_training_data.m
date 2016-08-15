% Luong Nguyen
% 5/16/14
% Examine training data 
% evaluate whether the 2nd and 3rd could serve as color opponency while the
% 1st svd is interpreted as brightness
% 

close all; clear all;
% read the training data, different options
workdir = '/Users/lun5/Research/color_deconvolution'; 
resultdir = fullfile(workdir, 'results',datestr(now,'yymmdd'));

if ~exist(resultdir,'dir')
    mkdir(resultdir);
end

training1 = load([workdir filesep 'training_pink_purple.mat'],'training_data');
training2 = load([workdir filesep 'results' filesep '140428' filesep  'training_pink_purple.mat'],'training_data');
training3 = load([workdir filesep 'results' filesep '140507' filesep  'training_pink_purple.mat'],'training_data');
%training_data = [training1.training_data training2.training_data];
training_data = training3.training_data(:,1:4000);
[U,D,V] = svd(training_data,0);
rotation_matrix = [-U(:,1) U(:,2:3)]'; % this is correct one

% pca of the training data
normalized_training = zscore(training_data'); % zcore to center at 0, sd = 1, rows = observations, cols = variables (RGB)
[pc, score, latent] = princomp(normalized_training); %pc = 3x3 matrix, each column is a pc
% latent, a vector containing the eigenvalues of the covariance matrix of
% normalized_training

disp('contribution of pcs to explaining the variance of data');
cumsum(latent)./sum(latent)
%biplot(pc,'Scores', score, 'Color',training_data'./255)

% plot data in PC1 PC2
h1 = figure; scatter(score(:,1), score(:,2), 20, training_data'./255,'filled');
title('{\bf PCA} using princomp','FontSize',15); xlabel('PC 1','FontSize',15); ylabel('PC 2','FontSize',15);
set(gca,'FontSize',15);
grid on
%print(h1,'-dtiff', [resultdir filesep 'PC1PC2scatter' datestr(now,'hhMM') '.tiff']);
%
h2 = figure; scatter(score(:,2), score(:,3), 20, training_data'./255,'filled');
title('{\bf PCA} using princomp','FontSize',15); xlabel('PC 2','FontSize',15); ylabel('PC 3','FontSize',15)
set(gca,'FontSize',15);grid on
%print(h2,'-dtiff', [resultdir filesep 'PC2PC3scatter' datestr(now,'hhMM') '.tiff']);

%
h3 = figure; scatter(score(:,1), score(:,3), 20, training_data'./255,'filled');
title('{\bf PCA} using princomp','FontSize',15); xlabel('PC 1','FontSize',15); ylabel('PC 3','FontSize',15)
set(gca,'FontSize',15);grid on
%print(h3,'-dtiff', [resultdir filesep 'PC1PC3scatter' datestr(now,'hhMM') '.tiff']);

% 3D scatter plot with PCs
%h = figure; scatter3(training_data(1,:)./255,training_data(2,:)./255,...
%    training_data(3,:)./255,20,training_data'./255,'filled');
h4 = figure; scatter3(normalized_training(:,1),normalized_training(:,2),...
   normalized_training(:,3),20,training_data'./255,'filled');
hold on
biplot(pc);
hold off
xlabel('R'); ylabel('G'); zlabel('B');
set(gca,'FontSize',15);
%saveas(h4, [resultdir filesep 'training3Dscatter' datestr(now,'hhMM')], 'fig');
%title('Training data for r','FontSize',15);

% plot the sic coordinates
options = struct('Normalize','on');
mu_s = 4; sigma_s = 2;
training_sic = rgb2sic(training_data,mu_s, sigma_s, rotation_matrix, options); 
% scatter plot of training data
h5 =figure; scatter(training_sic(1,:),training_sic(2,:),20,training_data'./255,'filled');
xlabel('s1','FontSize',15);ylabel('s2','FontSize',15);
line([-1 1],[0 0],'Color','k')
line([0 0],[-1 1],'Color','k')
hold off
title('Color opponency 2D representation','FontSize',15);
axis([-1 1 -1 1])
set(gca,'FontSize',15);grid on
print(h5,'-dtiff', [resultdir filesep 'trainingCOspace' datestr(now,'hhMM') '.tiff']);

% examine data if suitable for PCA
mean_training = normc(mean(training_data,2));
mean_contrib = mean_training'*training_data;
errors_mean = training_data - mean_training*mean_contrib;
sd_errors_mean = std(errors_mean); % this is similar to the residuals of training data
resid_mean = sqrt(sum(errors_mean.^2,1));

h6 = figure; 
scatter(mean_contrib, resid_mean,20,training_data'./255,'filled');
%scatter(mean_contrib, sd_errors_mean,20,training_data'./255,'filled');
xlabel('Fitted','FontSize',15);
ylabel('Residuals','FontSize',15);grid on
print(h6,'-dtiff', [resultdir filesep 'equalVarTest_mean' datestr(now,'hhMM') '.tiff']);

% normal plot of the residual
h7 = figure; qqplot(resid_mean);
set(gca,'FontSize',15);
print(h7,'-dtiff', [resultdir filesep 'qqplot_mean' datestr(now,'hhMM') '.tiff']);

% now project data onto the first pc
firstpc_contrib = pc(:,1)'*training_data;
errors_firstpc = training_data - pc(:,1)*firstpc_contrib;
sd_resid_firstpc = std(errors_firstpc);
resid_firstpc = sqrt(sum(errors_firstpc.^2,1));

h8 = figure;
scatter(firstpc_contrib, resid_firstpc,20,training_data'./255,'filled');
%scatter(firstpc_contrib, sd_resid_firstpc,20,training_data'./255,'filled');
xlabel('Fitted','FontSize',15);
ylabel('Residuals','FontSize',15);grid on
print(h8,'-dtiff', [resultdir filesep 'equalVarTest_firstpc' datestr(now,'hhMM') '.tiff']);

% normal plot of the residual
h9 = figure; qqplot(resid_firstpc);
set(gca,'FontSize',15);
print(h9,'-dtiff', [resultdir filesep 'qqplot_pc' datestr(now,'hhMM') '.tiff']);





