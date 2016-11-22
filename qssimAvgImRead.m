function [qssim_avg, qssim_matrix, qssim_total]  = qssimAvgImRead()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Quaternion Structural Similarity Index Averages              %
% Displays the QSSIM comparison matrix                                    %
%                                                                         %
% Uses QSSIM Index, Version 1.2                                           %
% Copyright(c) 2011 Amir Kolaman                                          %
% All Rights Reserved.                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
qssim_total = zeros(6, 1);
qssim_matrix = cell(5,231);
qssim_values = cell(6, 1);
    
for k=1:6
    list_methods = {'Luong', 'Macenko', 'Reinhard', 'Khan', 'Vahadane', 'Vahadane Fast'};
    im_dir = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\', list_methods{k});
    targ_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\Target\';
    source_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\';
    renorm_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Renorm Images\';
    imlist = dir(fullfile(im_dir,'*.tif')); imlist = {imlist.name}';
    targlist = dir(fullfile(targ_dir,'*.tif')); targlist = {targlist.name}';
    
    for m=1:length(targlist)
        targetname = targlist{m};
        target = imread(fullfile(targ_dir,targetname));
        for n=1:231 
            index =(n+(m-1)*231);
            imname = imlist{index};
            norm_source = imread(fullfile(im_dir,imname));
            im1 = target; im2 = norm_source;
%             im1 = source;
            %% Calculate QSSIM values and store for each method
            [mqssim, ~] = qssim(im1, im2);
            qssim_matrix{m,n} = mqssim;
            
            % Store QSSIM totals for each method for all images
            qssim_total(k) = qssim_total(k) + qssim_matrix{m,n};
        end
    end
    qssim_values{k, 1} = qssim_matrix;
end
%% Calculate QSSIM averages
num_values = sum(sum(~cellfun(@isempty,qssim_matrix),2));
qssim_avg = qssim_total ./num_values;
luong_avg = qssim_avg(1, 1);
macenko_avg = qssim_avg(2, 1);
reinhard_avg = qssim_avg(3, 1);
khan_leeds_avg = qssim_avg(4, 1);
vahadane_avg = qssim_avg(5, 1);
vahadane_fast_avg = qssim_avg(6, 1);

save('norm qssim 232 7-25-16.mat', 'qssim_avg', 'qssim_values', 'qssim_total',...
    'luong_avg', 'macenko_avg', 'reinhard_avg', 'khan_leeds_avg', 'vahadane_avg', 'vahadane_fast_avg');

end