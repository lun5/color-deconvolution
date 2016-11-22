% load('total_time 232 7-18-16.mat', 'time_matrix');
% time_matrix=time_matrix(~cellfun('isempty', time_matrix));
% for k=1:6
%     for i=1:1160
%         time_values{k, i} = time_matrix{i,1}(1, k);
%     end
% end
% time_values = cell2mat(time_values);
% time_total = sum(time_values, 2);
% time_avg = time_total/1160;
% luong_avg = time_avg(1);       
% macenko_avg = time_avg(2);
% reinhard_avg = time_avg(3);
% khan_avg = time_avg(4);
% vahadane_avg = time_avg(5);
% vahadane_fast_avg = time_avg(6);
% 
% save('time avg 232 7-18-16.mat', 'time_avg', 'luong_avg', 'macenko_avg', 'reinhard_avg', 'khan_avg', 'vahadane_avg', 'vahadane_fast_avg');

%     time_numbers = cell2mat(time_values); 
%     time_methods{k} = time_numbers;
%     time_std(k) = std(time_numbers);
% end
% luong_std = time_std(1);       
% macenko_std = time_std(2);
% reinhard_std = time_std(3);
% khan_std = time_std(4);
% vahadane_std = time_std(5);
% vahadane_fast_std = time_std(6);
% 
% save('time std 232 7-18-16.mat', 'time_std', 'luong_std', 'macenko_std', 'reinhard_std', 'khan_std', 'vahadane_std', 'vahadane_fast_std');
            
% load('norm emd 232 7-19-16.mat', 'metric_values');
% metric_values=metric_values(~cellfun('isempty', metric_values));
% for k=1:6
%     metric_matrix = cell(1155, 3);
%     for i=1:1155
%         for j=1:3
%             metric_matrix{i, j} = metric_values{i, 1}(k,j);
%         end
%     end
%     metric_matrix = cell2mat(metric_matrix); 
%     for j=1:3
%         metric_methods{k} = metric_matrix;
%         metric_std(k, j) = std(metric_matrix(:, j));
%     end
% end
% luong_std = metric_std(1, :);       
% macenko_std = metric_std(2, :);
% reinhard_std = metric_std(3, :);
% khan_std = metric_std(4, :);
% vahadane_std = metric_std(5, :);
% vahadane_fast_std = metric_std(6, :)
% 
% save('emd std 232 7-19-16.mat', 'metric_std', 'luong_std', 'macenko_std', 'reinhard_std', 'khan_std', 'vahadane_std', 'vahadane_fast_std');
   

% for k=1:6
%     metric_matrix = cell2mat(metric_values{k}); 
%     metric_std(k,1) = std2(metric_matrix);
% end
% luong_std = metric_std(1);       
% macenko_std = metric_std(2);
% reinhard_std = metric_std(3);
% khan_std = metric_std(4);
% vahadane_std = metric_std(5);
% vahadane_fast_std = metric_std(6);
% 
% save('chisq_std_232_7-18-16.mat', 'qssim_std', 'luong_std', 'macenko_std', 'reinhard_std', 'khan_std',...
%     'vahadane_std', 'vahadane_fast_std');

load('norm qssim 232 7-25-16.mat', 'qssim_values');

qssim_values=qssim_values(~cellfun('isempty', qssim_values));
qssim_matrix = cell(6,1);
for k=1:6
    qssim_matrix = cell2mat(reshape(qssim_values{k, 1}, [1,1155]));
    qssim_methods{k} = qssim_matrix;
    qssim_std(k) = std(qssim_matrix);
end
luong_std = qssim_std(1);       
macenko_std = qssim_std(2);
reinhard_std = qssim_std(3);
khan_std = qssim_std(4);
vahadane_std = qssim_std(5);
vahadane_fast_std = qssim_std(6);

save('norm std qssim 232 7-25-16.mat', 'qssim_std', 'luong_std', 'macenko_std', 'reinhard_std', 'khan_std',...
    'vahadane_std', 'vahadane_fast_std');

% load('emd 7-13-16 5 target.mat', 'metric_values');
% metric_values=metric_values(~cellfun('isempty', metric_values));
% for k=1:8
%     metric_matrix = cell(390, 3);
%     for i=1:390
%         for j=1:3
%             metric_matrix{i, j} = metric_values{i, 1}(k,j);
%         end
%     end
%     metric_matrix = cell2mat(metric_matrix); 
%     for j=1:3
%         metric_methods{k} = metric_matrix;
%         metric_std(k, j) = std(metric_matrix(:, j));
%     end
% end
% luong_std = metric_std(1, :);       
% macenko_std = metric_std(2, :);
% reinhard_std = metric_std(3, :);
% khan_std = metric_std(4, :);
% khan_hist_std = metric_std(5, :);
% khan_leeds_std = metric_std(6, :);
% vahadane_std = metric_std(7, :);
% vahadane_fast_std = metric_std(8, :);
% 
% save('emd_std_7-13-16 5 target.mat', 'metric_std', 'luong_std', 'macenko_std', 'reinhard_std', 'khan_std', 'khan_hist_std',...
%  'khan_leeds_std', 'vahadane_std', 'vahadane_fast_std');
%     
