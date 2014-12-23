%% match hue by cdf
% non saturated pixels
addpath(fullfile(matlabroot,'toolbox','CircStat2012a'))
target_theta_nonsat = target_theta(~ target_indx_saturated);
source_theta_nonsat = source_theta(~ source_indx_saturated);

target_theta_deg = radtodeg(target_theta);
binsize = 9;
edgevec = -180:binsize:180;
[~,target_theta_hist] = histc(target_theta_deg,edgevec);
target_theta_im = reshape(target_theta_hist,[target_xsize, target_ysize]);
figure; imagesc(target_theta_im);% should I convert to uint8?
colormap(hsv);colorbar;
axis square

%% plot the distribution
% figure(1);
% subplot(1,2,1);circ_plot(target_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,2,2);circ_plot(source_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% figure(2);
% subplot(1,2,1);circ_plot(target_theta_nonsat(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,2,2);circ_plot(source_theta_nonsat(:),'hist',[],40,true,true,'linewidth',2,'color','r');

%% matching
%figure; rose2(target_theta_nonsat(:),30);
%figure; rose2(source_theta_nonsat(:),30);
% convert theta's of source's and target's to [0 1]
target_theta_norm = target_theta_nonsat./(2*pi) + 0.5;
source_theta_norm = source_theta_nonsat./(2*pi) + 0.5;

% histogram equalization of source --> target
binranges = 0:0.01:1;
bincounts = histc(target_theta_norm, binranges);
source_theta_norm_eq = histeq(source_theta_norm,bincounts);
% convert back to between -pi and pi
source_theta_eq = source_theta;
source_theta_eq_nonsat = (source_theta_norm_eq - 0.5).*(2*pi);
source_theta_eq(~ source_indx_saturated) = source_theta_eq_nonsat;

source_theta_eq_deg = radtodeg(source_theta_eq);
binsize = 9;
edgevec = -180:binsize:180;
[~,source_theta_eq_hist] = histc(source_theta_eq_deg,edgevec);
source_theta_eq_im = reshape(source_theta_eq_hist,[source_xsize, source_ysize]);
figure; imagesc(source_theta_eq_im);% should I convert to uint8?
colormap(hsv);colorbar;
axis square

%% statistics! 
sprintf('Hue Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f  %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_cdf %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    circ_mean(source_theta_nonsat,[],2), circ_std(source_theta_nonsat,[],[],2), min(source_theta),...
    max(source_theta), circ_skewness(source_theta_nonsat,[],2), circ_kurtosis(source_theta_nonsat,[],2), ...
    circ_mean(target_theta_nonsat,[],2), circ_std(target_theta_nonsat,[],[],2), min(target_theta),...
    max(target_theta), circ_skewness(target_theta_nonsat,[],2), circ_kurtosis(target_theta_nonsat,[],2),...
    circ_mean(source_theta_eq_nonsat,[],2), circ_std(source_theta_eq_nonsat,[],[],2),...
    min(source_theta_eq), max(source_theta_eq),...
    circ_skewness(source_theta_eq_nonsat,[],2), circ_kurtosis(source_theta_eq_nonsat,[],2))

%% histogram
bincol = [0.8,0.8,0.8];
legendstr = {'source','target','matched source'};
A = {source_theta,target_theta,source_theta_eq};
figure; nhist(A,'legend',legendstr,'color',bincol,...
    'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','pdf');

% B = {source_theta_norm(~ source_indx_saturated),target_theta_norm(~target_indx_saturated),source_theta_norm_eq(~ source_indx_saturated)};
% %B = {source_theta_norm,target_theta_norm,source_theta_norm_eq};
% figure; nhist(B, 'legend',legendstr,'color',bincol,...
%     'separate','samebins','noerror','smooth','binfactor',0.2,'smooth','minx',0,'maxx',1,'pdf');