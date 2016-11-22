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
im_list = dir(fullfile(im_dir,'*.tif'));
im_list = {im_list.name}';
im_sizes = [512, 2048, 5000, 7500, 10000, 15000];
run_times = zeros(length(im_sizes), length(method_names), length(im_list));
for ii = 1: length(im_list)
    im_name = im_list{ii}(1:end-4);
    for jj = 1:(length(im_sizes)-1)
        source_im = imread(fullfile(im_dir,'cropped_im',[im_name '_' num2str(im_sizes(jj)) '.tif']));
        %start = tic; normMacenko = NormMacenko(source_im, target); run_times(jj, 1, ii) = toc(start);
        %start = tic; normReinhard = im2uint8(NormReinhard(source_im, target));run_times(jj,2,ii) = toc(start);
        %start = tic; normVahadaneFast = Demo_colornorm(source_im, target);run_times(jj,5,ii) = toc(start);
        start = tic; normLuongFast = NormLuong_Fast(source_im, target);run_times(jj,7,ii) = toc(start);
        if jj < 6
            start = tic; normSCDLeeds = NormSCDLeeds(source_im, target);run_times(jj,3,ii) = toc(start);
            start = tic; normVahadane = SNMFnorm(source_im, target);run_times(jj,4,ii) = toc(start);
            start = tic; normLuong = NormLuong(source_im, target);run_times(jj,6,ii) = toc(start);
            %imwrite(normSCDLeeds,fullfile(im_dir,'normalized_im',['Khan_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            %imwrite(normLuong,fullfile(im_dir,'normalized_im',['Luong_' im_name '_' num2str(im_sizes(jj)) '.tif']));
            %imwrite(normVahadane,fullfile(im_dir,'normalized_im',['Vahadane_' im_name '_' num2str(im_sizes(jj)) '.tif']));
        else
            run_times(jj,3,ii) = 800;
            run_times(jj,4,ii) = 800;
            run_times(jj,6,ii) = 800;
        end
        %imwrite(normMacenko,fullfile(im_dir,'normalized_im',['Macenko_' im_name '_' num2str(im_sizes(jj)) '.tif']));
        %imwrite(normReinhard,fullfile(im_dir,'normalized_im',['Reihard_' im_name '_' num2str(im_sizes(jj)) '.tif']));
        %imwrite(normVahadaneFast,fullfile(im_dir,'normalized_im',['VahadaneFast_' im_name '_' num2str(im_sizes(jj)) '.tif']));
        %imwrite(normLuongFast,fullfile(im_dir,'normalized_im',['LuongFast_' im_name '_' num2str(im_sizes(jj)) '.tif']));
        fprintf('Done with image size of %d\n',im_sizes(jj));
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
  
  for i = 1:size(med_run_times,1)
      fprintf('%d ',im_sizes(i));
      for j = 1:size(med_run_times,2)
          fprintf('&%.1f ',med_run_times(i,j));
      end
      fprintf('\\\\ \n')
  end
          
      




