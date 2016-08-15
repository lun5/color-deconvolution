function [metric_values] = histValidationImRead(metric)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric = 'emd'                                                          %
%   Earth Mover's Distance (EMD) between positive vectors (histograms).   %
% metric = 'chisq'                                                        %
%   The chi-squared distance between two vectors is defined as:           %
%    d(x,y) = sum( (xi-yi)^2 / (xi+yi) ) / 2;                             %
% metric = 'kl'                                                           %
%   Kullback-Leibler Divergence                                           %
%                                                                         %
% Uses Piotr's Computer Vision Matlab Toolbox      Version 2.52           %
% Copyright 2014 Piotr Dollar.  [pdollar-at-gmail.com]                    %
% Licensed under the Simplified BSD License [see external/bsd.txt]        %
% Uses Implementation of the Kullback-Leibler Divergence                  %                           
% @author: Nima Razavi @date: 2008                                        %
% @url: https://www.mathworks.com/matlabcentral/fileexchange/20688        %
%-kullback-leibler-divergence                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
rand('seed',0);
metric_total = zeros(6, 3);
metric_matrix = zeros(6, 3);
metric_values = cell(1155, 1);

for m=1:1155
    clf;
    for k=1:6
        method_names = {'Luong', 'Macenko', 'Reinhard', 'Khan', 'Vahadane', 'Vahadane Fast'};
        im_dir = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Renorm Images\', method_names{k});
        targ_dir = 'C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512\';
        imlist = dir(fullfile(im_dir,'*.tif')); imlist = {imlist.name}';
        targlist = dir(fullfile(targ_dir,'*.tif')); targlist = {targlist.name}';
        
        imname = imlist{m};
        filename = strsplit(imname, '-');
        %targetname = strcat(filename{1}, '.tif');
        targetname = filename{2};   %renorm
        target = imread(fullfile(targ_dir,targetname));
        norm_source = imread(fullfile(im_dir,imname));
        
            %% Hue, saturation, value histograms for normalized source image
            HSV1 = rgb2hsv(target);
            H1 = HSV1(:,:,1); S1 = HSV1(:,:,2); V1 = HSV1(:,:,3);
            HSV2 = rgb2hsv(norm_source);
            H2 = HSV2(:,:,1); S2 = HSV2(:,:,2); V2 = HSV2(:,:,3);
            
            %% Normalize histogram so area under curve = 1
            [hHist1, h1] = imhist(H1, 256); [sHist1, s1] = imhist(S1, 256); [vHist1, v1] = imhist(V1, 256);
            [hHist2, h2] = imhist(H2, 256); [sHist2, s2] = imhist(S2, 256); [vHist2, v2] = imhist(V2, 256);
            hHist1 = hHist1/sum(hHist1); sHist1 = sHist1/sum(sHist1); vHist1 = vHist1/sum(vHist1);
            hHist2 = hHist2/sum(hHist2); sHist2 = sHist2/sum(sHist2); vHist2 = vHist2/sum(vHist2);
            %% Calculate metric values 
            hmetric = pdist2(hHist1', hHist2', metric);
            smetric = pdist2(sHist1', sHist2', metric);
            vmetric = pdist2(vHist1', vHist2', metric);
            metric_total(k, 1) = metric_total(k, 1) + hmetric;
            metric_total(k, 2) = metric_total(k, 2) + smetric;
            metric_total(k, 3) = metric_total(k, 3) + vmetric;
            %% Store hue, saturation, and value totals for each method
            metric_matrix(k, 1) = hmetric;
            metric_matrix(k, 2) = smetric;
            metric_matrix(k, 3) = vmetric;
            
            %% Display histogram montage
%             h_figure = figure(1); h_figure.Position = [150, 100, 1500, 800];
%             subplot_tight(2, 3, k, .07);
%             plot(h1, hHist1);        %bar(h1, hHist1, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.3);
%             ylim([0 0.04]); title(method_names{k}); hold on;
%             plot(h2, hHist2);
%             legend('Target', 'Normalized Source','Location','northwest');
%             %legend('Original Source', 'Renormalized Source','Location','northwest'); %Renorm
%             
%             s_figure = figure(2);  s_figure.Position = [150, 100, 1500, 800];
%             subplot_tight(2, 3, k, .07);
%             plot(s1, sHist1);
%             ylim([0 0.04]); title(method_names{k}); hold on;
%             plot(s2, sHist2);
%             legend('Target', 'Normalized Source','Location','northwest');
%             %legend('Original Source', 'Renormalized Source','Location','northwest'); %Renorm
%             
%             v_figure = figure(3); v_figure.Position = [150, 100, 1500, 800];
%             subplot_tight(2, 3, k, .07);
%             plot(v1, vHist1);
%             ylim([0 0.04]); title(method_names{k}); hold on;
%             plot(v2, vHist2);
%             legend('Target', 'Normalized Source','Location','northwest');
%             %legend('Original Source', 'Renormalized Source','Location','northwest'); %Renorm
        end
        
%         figure(1); suptitle('Hue Histogram');
%         targetfile = strrep(targetname, '.tif', '-');
%         filename = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\Chi-Square\', strrep(imname, '.tif', '-Hue.tif'));
%         print(h_figure, '-dtiff', filename);
%         figure(2); suptitle('Saturation Histogram');
%         targetfile = strrep(targetname, '.tif', '-');
%         filename = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\Chi-Square\', strrep(imname, '.tif', '-Sat.tif'));
%         print(s_figure, '-dtiff', filename);
%         figure(3); suptitle('Value Histogram');
%         targetfile = strrep(targetname, '.tif', '-');
%         filename = strcat('C:\Users\tam128\Documents\MATLAB\Norm\Tiles_512 Images\Chi-Square\', strrep(imname, '.tif', '-Val.tif'));
%         print(v_figure, '-dtiff', filename);
        
        % Store hue, saturation, value for all methods for each image
        metric_values{m} = metric_matrix;
end

%% Calculate average metric values 
num_values = sum(sum(~cellfun(@isempty,metric_values), 1));
metric_avg = metric_total./num_values;
luong_avg = metric_avg(1, :);       
macenko_avg = metric_avg(2, :);
reinhard_avg = metric_avg(3, :);
khan_avg = metric_avg(4, :);
vahadane_avg = metric_avg(5, :);
vahadane_fast_avg = metric_avg(6, :);

save('renorm EMD 232 7-25-16.mat', 'metric_avg', 'metric_values', 'metric_total',...
    'luong_avg', 'macenko_avg', 'reinhard_avg', 'khan_avg', 'vahadane_avg', 'vahadane_fast_avg');

end
