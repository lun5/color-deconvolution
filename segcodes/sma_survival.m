% load(fullfile(matlabroot,'examples','stats','readmissiontimes.mat'));
% female = [ReadmissionTime(Sex==1),Censored(Sex==1)];
% male = [ReadmissionTime(Sex==0),Censored(Sex==0)];
% 
% figure()
% ecdf(gca,female(:,1),'Censoring',female(:,2));
% hold on
% [f,x] = ecdf(male(:,1),'Censoring',male(:,2));
% stairs(x,f,'--r')
% hold off
% legend('female','male','Location','SouthEast')
% 
% figure()
% ax1 = gca;
% ecdf(ax1,female(:,1),'Censoring',female(:,2),'function','survivor');
% hold on
% [f,x] = ecdf(male(:,1),'Censoring',male(:,2),'function','survivor');
% stairs(x,f,'--r')
% legend('female','male')

%% script to examine SMA associating to survival 
% Assumption is that epithelial cells are from the tumor areas
% April 17, 2017
close all;
data_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data';
coordinate_dir = fullfile(data_dir,'coordinates_all_bm');
clinical_data = readtable(fullfile(data_dir,'clincal_data_june27.csv'),...
    'Delimiter',',');

output_dir = 'D:\Documents\multiplex\Immuno_bm';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% find the index of SMA
filelist = clinical_data.spot_name;
tmp = load(fullfile(coordinate_dir,[filelist{1} '.mat']));
bm_names = cellstr(tmp.bm_names);

curr_bm_names = {'Median.Cell.CD8'};
indx = ismember(bm_names,curr_bm_names);

sma_expr = zeros(length(filelist),3); % col1: epi, col2: stroma;

for i = 1:length(filelist)
    imname = filelist{i};
    
    if ~exist(fullfile(coordinate_dir, [imname '.mat']),'file')
       sma_expr(i,:) = -1;
       continue
    end
    
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    %x = max(1,tmp.x);
    %y = max(1,tmp.y);
    
    %indx_cells = sub2ind(size(seg),y,x);
    epi_stroma = tmp.epithelial == 1;
    %area = tmp.area;
    % stroma, val = 2
    bm_data = tmp.bm_data(~epi_stroma,indx);
    %num_cells = size(tmp.bm_data,1);
    sma_expr(i,1) = mean(tmp.bm_data(epi_stroma, indx));
    sma_expr(i,2) = mean(tmp.bm_data(~epi_stroma, indx));
    sma_expr(i,3) = mean(tmp.bm_data(:, indx));
end
output_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data';
exist_indx = sma_expr(:,1) > -1;
figure; histogram(sma_expr(exist_indx,1),'Normalization','probability');
xlim([0 10]); ylim([0 0.5]); xlabel('Epithelial SMA');
figure; histogram(sma_expr(exist_indx,2),'Normalization','probability');
xlim([0 10]); ylim([0 0.5]); xlabel('Stroma SMA');

sma_level = sma_expr(:,2) > prctile(sma_expr(:,2),50);
%sma_level = sma_expr(:,2) > median(sma_expr(:,2));
high_level_indx = exist_indx & sma_level;
figure()
ax1 = gca;
ecdf(ax1,clinical_data.recurtime5yrs_days(high_level_indx),...
    'function','survivor');
hold on
[f,x] = ecdf(clinical_data.recurtime5yrs_days(~high_level_indx),...
    'function','survivor');
stairs(x,f,'--r')
xlabel('Recurrent time');title(curr_bm_names);legend('high','low');
title(curr_bm_names)
print(fullfile(output_dir,[curr_bm_names{1} '_rec']) ,'-dpng');

figure()
ax1 = gca;
ecdf(ax1,clinical_data.survtime5yrs_days(high_level_indx),...
    'function','survivor');
hold on
[f,x] = ecdf(clinical_data.survtime5yrs_days(~high_level_indx),...
    'function','survivor');
stairs(x,f,'--r')
xlabel('Survival time');title(curr_bm_names);legend('high','low');
legend('high','low')
title(curr_bm_names)
print(fullfile(output_dir,[curr_bm_names{1} '_surv']) ,'-dpng');
