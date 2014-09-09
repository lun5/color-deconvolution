% match statistics of Hue
% Hue is a circular value
% use CircStat2012a for this purpose
addpath(fullfile(matlabroot,'toolbox','CircStat2012a'))
target_theta_nonsat = target_theta(~target_indx_saturated);
source_theta_nonsat = source_theta(~source_indx_saturated);

% statistics of target
target_mean = circ_mean(target_theta_nonsat,[],2);
target_std = circ_std(target_theta_nonsat,[],[],2);
target_skewness = circ_skewness(target_theta_nonsat,[],2);
target_kurtosis = circ_kurtosis(target_theta_nonsat,[],2);

% statistics of source
source_mean = circ_mean(source_theta_nonsat,[],2);
source_std = circ_std(source_theta_nonsat,[],[],2);
source_skewness = circ_skewness(source_theta_nonsat,[],2);
source_kurtosis = circ_kurtosis(source_theta_nonsat,[],2);
% %% plot the distribution
% figure(1);
% subplot(1,2,1);circ_plot(target_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,2,2);circ_plot(source_theta(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% figure(2);
% subplot(1,2,1);circ_plot(target_theta_nonsat(:),'hist',[],40,true,true,'linewidth',2,'color','r');
% subplot(1,2,2);circ_plot(source_theta_nonsat(:),'hist',[],40,true,true,'linewidth',2,'color','r');

%figure; rose2(target_theta_nonsat(:),30);
%figure; rose2(source_theta_nonsat(:),30);

%% center to 0
source_eq = source_theta_nonsat - source_mean;
wrap_theta;

%% spread source to have same std as target
source_eq = source_eq * target_std/circ_std(source_eq,[],[],2);
wrap_theta;

%% fix mean to be equal target mean
source_eq = source_eq + target_mean;
wrap_theta;
%% fix skewness
[source_eq, ~] = modskew(source_eq,target_skewness);
wrap_theta;

%% fix kurtosis
[source_eq, ~] = modkurt(source_eq,target_skewness);
wrap_theta;

%% put it back together
source_theta_eq = source_theta;
source_theta_eq(~source_indx_saturated) = source_eq;

% cmap = [1 0 0; 0 0 1; 1 1 1];% ; 0 1 1];
% source_theta_eq_deg = radtodeg(source_theta_eq);
% binsize = 30;
% edgevec = -150:binsize:-60;
% [~,source_theta_eq_hist] = histc(source_theta_eq_deg,edgevec);
% source_theta_eq_im = reshape(source_theta_eq_hist,[source_xsize, source_ysize]);
% figure; imagesc(source_theta_eq_im);% should I convert to uint8?
% colormap(cmap);%colorbar;
% axis square
% source_theta_deg = radtodeg(source_theta);
% binsize = 30;
% edgevec = -180:binsize:-90;
% [~,source_theta_hist] = histc(source_theta_deg,edgevec);
% source_theta_im = reshape(source_theta_hist,[source_xsize, source_ysize]);
% figure; imagesc(source_theta_im);% should I convert to uint8?
% colormap(cmap);%colorbar;
% axis square

% close all;
% statistics
sprintf('Hue Mean std min max skewness kurtosis\nSource %.4f %.4f %.4f %.4f %.4f %.4f\nTarget %.4f %.4f %.4f %.4f %.4f %.4f\nSource_statistics %.4f %.4f %.4f %.4f %.4f %.4f\n',...
    source_mean, source_std, min(source_theta), max(source_theta), source_skewness, source_kurtosis,...
    target_mean, target_std, min(target_theta_nonsat), max(target_theta_nonsat), target_skewness, target_kurtosis,...
    circ_mean(source_eq,[],2), circ_std(source_eq,[],[],2),min(source_eq,[],2), max(source_eq,[],2),...
    circ_skewness(source_eq,[],2), circ_kurtosis(source_eq,[],2))