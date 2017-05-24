%% script for run time
% im_dir = 'D:\Documents\Tiles_Norm\time_test';
% im_list = dir(fullfile(im_dir,'*.tif'));
% im_list = {im_list.name}';
% 5k:  rect = [36100 10000 4999 4999];
% 7.5k rect = [12000 23000 7499 7499];
% 10k rect = [10400 4000 9999 9999];
% 15k rect = [28000 5000 14999 14999];

% im_list = {'2ale5ngryfnpo.tif','jbakl4tseqt.tif','5k.tif','75k.tif','10k.tif','15k.tif'};
method_names = {'MK','RH','KH','VH','VHF','SCAN','FSCAN'};
target = imread(fullfile('D:\Documents\Tiles_Norm\Target_512', 'ocmmhhrtzz5.tif'));
% num_run = 10;
% %run_times = cell(num_run,1); %(zeros(length(im_list), length(method_names));
% im_sizes = zeros(size(im_list));

im_dir = 'D:\Documents\Tiles_Norm\time_test_2';
if ~exist(fullfile(im_dir,'normalized_im'),'dir');
    mkdir(fullfile(im_dir,'normalized_im'),'dir');
end
%im_list = dir(fullfile(im_dir,'original_images','*.tif'));
im_list = dir(fullfile(im_dir, '*.tif'));
im_list = {im_list.name}';
im_sizes = [512, 2048, 5000, 7500, 10000, 15000];
%run_times = zeros(length(im_sizes), length(method_names), length(im_list));
for ii = 2:length(im_list)
    im_name = im_list{ii}(1:end-4);
    for jj = 6%1:(length(im_sizes)-1)
        source_im = imread(fullfile(im_dir,'cropped_im',[im_name '_' num2str(im_sizes(jj)) '.tif']));
        if ~exist(fullfile(im_dir,'normalized_im',['Macenko_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with Macenko size %d\n', im_sizes(jj))
            start = tic; normMacenko = NormMacenko(source_im, target); run_times(jj, 1, ii) = toc(start);
            imwrite(normMacenko,fullfile(im_dir,'normalized_im',['Macenko_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear normMacenko
            fprintf('Done with Macenko size %d\n', im_sizes(jj))
        else
            fprintf('Already done with Macenko size %d\n', im_sizes(jj))
        end
        
        if ~exist(fullfile(im_dir,'normalized_im',['Reinhard_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with Reinhard size %d\n', im_sizes(jj))
            start = tic; normReinhard = im2uint8(NormReinhard(source_im, target));run_times(jj,2,ii) = toc(start);
            imwrite(normReinhard,fullfile(im_dir,'normalized_im',['Reinhard_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear normReinhard
            fprintf('Done with Reinhard size %d\n', im_sizes(jj))
        else
            fprintf('Already with Reinhard size %d\n', im_sizes(jj))
        end
          
        if ~exist(fullfile(im_dir,'normalized_im',['Khan_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with Khan size %d', im_sizes(jj));          
            start = tic; normSCDLeeds = NormSCDLeeds(source_im, target);run_times(jj,3,ii) = toc(start);
            imwrite(normSCDLeeds,fullfile(im_dir,'normalized_im',['Khan_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear normSCDLeeds;  
            fprintf('Done with Khan size %d', im_sizes(jj))   
        end
        
%         if ~exist(fullfile(im_dir,'normalized_im',['SCPN_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
%             fprintf('Start with VH size %d\n', im_sizes(jj))
%             start = tic; normVahadane = SNMFnorm(source_im, target);run_times(jj,4,ii) = toc(start);
%             imwrite(normVahadane,fullfile(im_dir,'normalized_im',['SCPN_' im_name '_' num2str(im_sizes(jj)) '.tif']));
%             clear normVahadane
%             fprintf('Done with VH size %d\n', im_sizes(jj))
%         else
%             fprintf('Already done with VH size %d\n', im_sizes(jj))
%         end
        
        if ~exist(fullfile(im_dir,'normalized_im',['FSCPN_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with VHF size %d', im_sizes(jj))
            start = tic; normVahadaneFast = Demo_colornorm(source_im, target);run_times(jj,5,ii) = toc(start);
            imwrite(normVahadaneFast,fullfile(im_dir,'normalized_im',['FSCPN_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear normVahadaneFast
            fprintf('Done with VHF size %d', im_sizes(jj))
        end
        
        if ~exist(fullfile(im_dir,'normalized_im',['SCAN_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with SCAN size %d\n', im_sizes(jj))
            start = tic; normLuong = NormLuong(source_im, target);run_times(jj,6,ii) = toc(start);
            imwrite(normLuong,fullfile(im_dir,'normalized_im',['SCAN_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear normLuong
            fprintf('Done with SCAN size %d\n', im_sizes(jj))        
        else
            fprintf('Already one with SCAN size %d\n', im_sizes(jj))   
        end
        
        if ~exist(fullfile(im_dir,'normalized_im',['FSCAN_' im_name '_' num2str(im_sizes(jj)) '.tif']),'file')
            fprintf('Start with FSCAN size %d', im_sizes(jj));
            start = tic; normLuongFast = NormLuong_Fast(source_im, target);run_times(jj,7,ii) = toc(start);
            imwrite(normLuongFast,fullfile(im_dir,'normalized_im',['FSCAN_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            clear nomLuongFast
            fprintf('Done with FSCAN size %d\n', im_sizes(jj));
        else
            fprintf('Already done with FSCAN size %d\n', im_sizes(jj))   
        end
        
        clear source_im
        fprintf('Done with image %s size of %d\n\n',im_name, im_sizes(jj));
        save(fullfile(im_dir, 'run_times_dec7.mat'),'run_times');
    end
    fprintf('Done with image %s\n',im_name);
    %save('run_times_sept22.mat','run_times');
end

%save('run_times_sept13.mat','run_times');
    %source_im = imread(fullfile(im_dir,im_list{ii}));
    %im_sizes(ii) = max(size(source_im));
% 
%     for rr = 1:num_run
%         if ii == 1; run_times{rr} = zeros(length(im_list), length(method_names)); end
%         start = tic; normMacenko = NormMacenko(source_im, target); run_times{rr}(ii,1) = toc(start);
%         start = tic; normReinhard = im2uint8(NormReinhard(source_im, target));run_times{rr}(ii,2) = toc(start);
%         
%         if ii < 6
%             start = tic; normSCDLeeds = NormSCDLeeds(source_im, target);run_times{rr}(ii,3) = toc(start);
%             start = tic; normVahadane = SNMFnorm(source_im, target);run_times{rr}(ii,4) = toc(start);
%             start = tic; normLuong = NormLuong(source_im, target);run_times{rr}(ii,6) = toc(start);
%             %imwrite(normVahadane,fullfile(im_dir,['Vahadane_' im_list{ii}]));
%         else
%             run_times{rr}(ii,3) = 900;
%             run_times{rr}(ii,4) = 900;
%             run_times{rr}(ii,6) = 900;
%         end
%         start = tic; normVahadaneFast = Demo_colornorm(source_im, target);run_times{rr}(ii,5) = toc(start);
%         start = tic; normLuongFast = NormLuong_Fast(source_im, target);run_times{rr}(ii,7) = toc(start);
%         fprintf('Done with %d iteration\n',rr);
%     end
%     fprintf('Done with image %s\n', im_list{ii});
  %end

% load 'run_times_sept21.mat'
med_run_times = median(run_times,3);
save(fullfile(im_dir, 'run_times_dec9.mat'),'run_times_2');  
%   for i = 1:size(med_run_times,1)
%       fprintf('%d ',im_sizes(i));
%       for j = 1:size(med_run_times,2)
%           fprintf('&%.1f ',med_run_times(i,j));
%       end
%       fprintf('\\\\ \n')
%   end
          
% %% try with larger images
% addpath(genpath('C:\Users\luong_nguyen\Documents\GitHub\color-deconvolution'));
% source_im = imread('D:\Documents\PMI_classify\sample_adh05\adh07-1c.tif');
% target = imread(fullfile('D:\Documents\Tiles_Norm\Target_512', 'ocmmhhrtzz5.tif'));
% 
% run_times = zeros(4,1);
% start = tic; normMacenko = NormMacenko(source_im, target); run_times(1) = toc(start);
% start = tic; normReinhard = im2uint8(NormReinhard(source_im, target));run_times(2) = toc(start);
% start = tic; normVahadaneFast = Demo_colornorm(source_im, target);run_times(3) = toc(start);
% start = tic; normLuongFast = NormLuong_Fast(source_im, target);run_times(4) = toc(start);
% 
%% input images dir
input_dir = 'D:\Documents\Tiles_Norm\time_test_2\normalized_im';
output_dir = 'D:\Documents\Tiles_Norm\time_test_2\resize_images';
%output_dir = 'D:\Documents\Tiles_Norm\time_test_2\images';
%imlist = dir(fullfile(input_dir,'*.tif'));
imlist = dir(fullfile(input_dir,'*.tif'));
imlist = {imlist.name}';

parfor i = 1:length(imlist)
   if exist(fullfile(output_dir, [imlist{i}(1:end-3) 'tif']),'file')
       continue;
   end
   im = imread(fullfile(input_dir,imlist{i}));
   im_rz = imresize(im,[512 512]);
   imwrite(im_rz,fullfile(output_dir, [imlist{i}(1:end-3) 'tif']));
   %imwrite(im_rz,fullfile(output_dir, [imlist{i}(1:end-3) 'png']));
end

source_im_dir = 'D:\Documents\Tiles_Norm\time_test_2\cropped_im';
output_dir = 'D:\Documents\Tiles_Norm\time_test_2\rz_cropped_im';
imlist = dir(fullfile(source_im_dir,'*.tif'));
imlist = {imlist.name}';

parfor i = 1:length(imlist)
   imname = imlist{i}(1:end-4);
   if exist(fullfile(output_dir, [imname '.tif']),'file')
       continue;
   end
   im = imread(fullfile(source_im_dir, [imname '.tif']));
   im_rz = imresize(im,[512 512]);
   imwrite(im_rz, fullfile(output_dir, [imname '.tif']));
end


